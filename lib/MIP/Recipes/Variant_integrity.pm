package MIP::Recipes::Variant_integrity;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use autodie qw{ :all };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Spec::Functions qw{ catdir catfile devnull };

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_variant_integrity };

}

##Constants
Readonly my $SPACE   => q{ };
Readonly my $NEWLINE => qq{\n};

sub analysis_variant_integrity {

## analysis_variant_integrity

## Function : Tests sample for correct relatives (only performed for samples with relatives defined in pedigree file) performed on sequence data.
## Returns  : ""
## Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, family_id, $temp_directory, $outaligner_dir, $call_type
##          : $parameter_href             => Parameter hash {REF}
##          : $active_parameter_href      => Active parameters for this analysis hash {REF}
##          : $sample_info_href           => Info on samples and family hash {REF}
##          : $file_info_href             => The file_info hash {REF}
##          : $infile_lane_prefix_href    => Infile(s) without the ".ending" {REF}
##          : $job_id_href                => Job id hash {REF}
##          : $program_name               => Program name
##          : $family_id                  => Family id
##          : $temp_directory             => Temporary directory
##          : $outaligner_dir             => Outaligner_dir used in the analysis
##          : $call_type                  => The variant call type

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id;
    my $temp_directory;
    my $outaligner_dir;
    my $call_type;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Program::Variantcalling::Variant_integrity
      qw{ variant_integrity_mendel variant_integrity_father };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info};
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_family_dead_end };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $family_id,
            program_name          => $program_name,
            program_directory =>
              catfile( lc $outaligner_dir, q{casecheck}, lc $program_name ),
            call_type    => $call_type,
            core_number  => $core_number,
            process_time => $time,
        }
    );

    my ( $volume, $directory, $program_info_file ) =
      File::Spec->splitpath($program_info_path)
      ;    #Split to enable submission to &sample_info_qc later

    #To enable submission to &sample_info_qc later
    my $stderr_file = $program_info_file . q{.stderr.txt};

    #To enable submission to &sa
    my $stdout_file = $program_info_file . q{.stdout.txt};

    ## Assign Directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir );
    my $outfamily_directory = catfile( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir, q{casecheck}, lc $program_name );
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
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
            jobid_chain =>
              $parameter_href->{pgatk_combinevariantcallsets}{chain},
        }
    );

    my $infile_path =
      catfile( $infamily_directory, $infile_prefix . $infile_suffix . q{*} );

    my $family_file =
      catfile( $outfamily_file_directory, $family_id . q{.fam} );
    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            parameter_href        => $parameter_href,
            active_parameter_href => $active_parameter_href,
            sample_info_href      => $sample_info_href,
            FILEHANDLE            => $FILEHANDLE,
            fam_file_path         => $family_file,
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
    if ( scalar( @{ $active_parameter_href->{sample_ids} } ) > 1 )
    {    #Only perform if more than 1 sample

        if ( $parameter_href->{dynamic_parameter}{trio} ) {

            variant_integrity_mendel(
                {
                    infile_path  => $file_path_prefix . $infile_suffix,
                    outfile_path => catfile(
                        $outfamily_directory, $family_id . q{_mendel.txt}
                    ),
                    family_file => $family_file,
                    family_type =>
                      $active_parameter_href->{genmod_models_family_type},
                    FILEHANDLE => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            if ( $mip_program_mode == 1 ) {

                ## Collect QC metadata info for later use
                add_program_outfile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        program_name     => q{variant_integrity_mendel},
                        outdirectory     => $outfamily_directory,
                        outfile          => $family_id . q{_mendel.txt},
                    }
                );
            }
        }

    }

    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        #Alias
        my $father_info =
          $sample_info_href->{sample}{$sample_id}{father};

        #Father is included in analysis
        if ( $father_info ne q{0} ) {

            variant_integrity_father(
                {
                    infile_path  => $file_path_prefix . $infile_suffix,
                    outfile_path => catfile(
                        $outfamily_directory, $family_id . q{_father.txt}
                    ),
                    family_file => $family_file,
                    family_type =>
                      $active_parameter_href->{genmod_models_family_type},
                    FILEHANDLE => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            if ( $mip_program_mode == 1 ) {

                ## Collect QC metadata info for later use
                add_program_outfile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        program_name     => q{variant_integrity_father},
                        outdirectory     => $outfamily_directory,
                        outfile          => $family_id . q{_father.txt},
                    }
                );
            }
        }
    }

    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {

        slurm_submit_job_sample_id_dependency_family_dead_end(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $family_id,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }

    return;

}
1;
