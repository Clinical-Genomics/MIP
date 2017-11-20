package MIP::Recipes::Analysis::Peddy;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile splitpath };
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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_peddy };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_peddy {

## Function : Compares familial-relationships and sexes.
## Returns  :
## Arguments: $parameter_href          => Parameter hash {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $infamily_directory      => In family directory
##          : $outfamily_directory     => Out family directory
##          : $temp_directory          => Temporary directory
##          : $program_name            => Program name
##          : $family_id               => Family id
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $call_type               => Variant call type

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $infamily_directory;
    my $outfamily_directory;
    my $program_name;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
    my $call_type;
    my $temp_directory;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
        },
        temp_directory => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory,
        },
        infamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infamily_directory,
        },
        outfamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfamily_directory,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Pedigree qw{ create_fam_file };
    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_view_and_index_vcf };
    use MIP::Program::Variantcalling::Peddy qw{ peddy };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::QC::Record
      qw{ add_program_outfile_to_sample_info add_program_metafile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    my $program_directory =
      catfile( $outaligner_dir, q{casecheck}, $program_name );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            job_id_href                     => $job_id_href,
            FILEHANDLE                      => $FILEHANDLE,
            directory_id                    => $family_id,
            program_name                    => $program_name,
            program_directory               => $program_directory,
            call_type                       => $call_type,
            core_number                     => $core_number,
            process_time                    => $time,
            source_environment_commands_ref => [],
        }
    );

    # Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $program_info_file ) =
      splitpath($program_info_path);

    # To enable submission to &sample_info_qc later
    my $stderr_file = $program_info_file . $DOT . q{stderr.txt};

    #To enable submission to &sample_info_qc later
    my $stdout_file = $program_info_file . $DOT . q{stdout.txt};

    ## Paths
    my $outfamily_file_directory =
      catfile( $active_parameter_href->{outdata_dir}, $family_id );
    my $outfile_path_prefix = catfile( $outfamily_directory, $family_id );

    ## Files
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

    my $suffix = $DOT . q{vcf.gz};

    my $family_file =
      catfile( $outfamily_file_directory, $family_id . $DOT . q{fam} );

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

    my $infile_path = catfile( $infamily_directory,
        $infile_prefix . $infile_suffix . $ASTERISK );
    migrate_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $infile_path,
            outfile_path => $temp_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Reformat variant calling file and index
    bcftools_view_and_index_vcf(
        {
            infile_path         => $file_path_prefix . $infile_suffix,
            outfile_path_prefix => $file_path_prefix,
            output_type         => q{z},
            index               => 1,
            index_type          => q{tbi},
            FILEHANDLE          => $FILEHANDLE,
        }
    );

    ## Peddy
    peddy(
        {
            infile_path         => $file_path_prefix . $suffix,
            outfile_prefix_path => $outfile_path_prefix,
            family_file_path    => $family_file,
            FILEHANDLE          => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    if ( $mip_program_mode == 1 ) {

        my %peddy_output = (
            ped_check => q{csv},
            sex_check => q{csv},
            peddy     => q{ped},
        );

      PEDDY_OUTPUT_FILES:
        while ( my ( $file_key, $temp_suffix ) = each %peddy_output ) {

            my $outfile_suffix = $DOT . $file_key . $DOT . $temp_suffix;

            ## Collect QC metadata info for later use
            add_program_metafile_to_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    program_name     => $program_name,
                    metafile_tag     => $file_key,
                    path             => $outfile_path_prefix . $outfile_suffix,
                }
            );
        }

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

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    return;
}

1;
