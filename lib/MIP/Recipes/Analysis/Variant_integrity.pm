package MIP::Recipes::Analysis::Variant_integrity;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile devnull splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_variant_integrity };

}

## Constants
Readonly my $SPACE      => q{ };
Readonly my $NEWLINE    => qq{\n};
Readonly my $DOT        => q{.};
Readonly my $UNDERSCORE => q{_};
Readonly my $ASTERISK   => q{*};

sub analysis_variant_integrity {

## Function : Tests sample for correct relatives (only performed for samples with relatives defined in pedigree file) performed on sequence data.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => The variant call type
##          : $family_id               => Family id
##          : $file_info_href          => The file_info hash {REF}
##          : $infamily_directory      => In family directory
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outfamily_directory     => Out family directory
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infamily_directory;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $outfamily_directory;
    my $parameter_href;
    my $program_name;
    my $sample_info_href;

    ## Default(s)
    my $call_type;
    my $family_id;
    my $outaligner_dir;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        call_type =>
          { default => q{BOTH}, store => \$call_type, strict_type => 1, },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        infamily_directory => {
            defined     => 1,
            required    => 1,
            store       => \$infamily_directory,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir,
            strict_type => 1,
        },
        outfamily_directory => {
            defined     => 1,
            required    => 1,
            store       => \$outfamily_directory,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Pedigree qw{ create_fam_file };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_family_dead_end };
    use MIP::Program::Variantcalling::Variant_integrity
      qw{ variant_integrity_mendel variant_integrity_father };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my ( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            call_type             => $call_type,
            core_number           => $core_number,
            directory_id          => $family_id,
            FILEHANDLE            => $FILEHANDLE,
            job_id_href           => $job_id_href,
            process_time          => $time,
            program_directory =>
              catfile( $outaligner_dir, q{casecheck}, $program_name ),
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
        }
    );

    # Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $program_info_file ) =
      splitpath($program_info_path);

    # To enable submission to &sample_info_qc later
    my $stderr_file = $program_info_file . $DOT . q{stderr.txt};

    # To enable submission to &sa
    my $stdout_file = $program_info_file . $DOT . q{stdout.txt};

    ## Assign Directories
    my $outfamily_file_directory =
      catfile( $active_parameter_href->{outdata_dir}, $family_id );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$family_id}{pgatk_combinevariantcallsets}{file_tag};
    my $infile_prefix = $family_id . $infile_tag . $call_type;
    my $file_path_prefix = catfile( $temp_directory, $infile_prefix );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain_vcf_data
    my $infile_suffix = get_file_suffix(
        {
            jobid_chain =>
              $parameter_href->{pgatk_combinevariantcallsets}{chain},
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

    my $infile_path =
      catfile( $infamily_directory,
        $infile_prefix . $infile_suffix . $ASTERISK );

    my $family_file =
      catfile( $outfamily_file_directory, $family_id . $DOT . q{fam} );
    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            active_parameter_href => $active_parameter_href,
            fam_file_path         => $family_file,
            FILEHANDLE            => $FILEHANDLE,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $infile_path,
            outfile_path => $temp_directory
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Variant_integrity
    if ( scalar @{ $active_parameter_href->{sample_ids} } > 1 ) {
        ## Only perform if more than 1 sample

        if ( $parameter_href->{dynamic_parameter}{trio} ) {

            variant_integrity_mendel(
                {
                    family_file => $family_file,
                    family_type =>
                      $active_parameter_href->{genmod_models_family_type},
                    FILEHANDLE   => $FILEHANDLE,
                    infile_path  => $file_path_prefix . $infile_suffix,
                    outfile_path => catfile(
                        $outfamily_directory,
                        $family_id . $UNDERSCORE . q{mendel.txt}
                    ),
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            if ( $mip_program_mode == 1 ) {

                ## Collect QC metadata info for later use
                add_program_outfile_to_sample_info(
                    {
                        path => catfile(
                            $outfamily_directory,
                            $family_id . $UNDERSCORE . q{mendel.txt}
                        ),
                        program_name     => q{variant_integrity_mendel},
                        sample_info_href => $sample_info_href,
                    }
                );
            }
        }

    }

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        # Alias
        my $father_info =
          $sample_info_href->{sample}{$sample_id}{father};

        # Father is included in analysis
        if ($father_info) {

            variant_integrity_father(
                {
                    family_file => $family_file,
                    family_type =>
                      $active_parameter_href->{genmod_models_family_type},
                    FILEHANDLE   => $FILEHANDLE,
                    infile_path  => $file_path_prefix . $infile_suffix,
                    outfile_path => catfile(
                        $outfamily_directory,
                        $family_id . $UNDERSCORE . q{father.txt}
                    ),
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            if ( $mip_program_mode == 1 ) {

                ## Collect QC metadata info for later use
                add_program_outfile_to_sample_info(
                    {
                        path => catfile(
                            $outfamily_directory,
                            $family_id . $UNDERSCORE . q{father.txt}
                        ),
                        program_name     => q{variant_integrity_father},
                        sample_info_href => $sample_info_href,
                    }
                );
            }
        }
    }

    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {

        slurm_submit_job_sample_id_dependency_family_dead_end(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                sbatch_file_name => $file_path,
            }
        );
    }
    return;
}

1;
