package MIP::Recipes::Analysis::Picardtools_collectmultiplemetrics;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile devnull };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
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
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_picardtools_collectmultiplemetrics };

}

## Constants
Readonly my $ASTERIX    => q{*};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_picardtools_collectmultiplemetrics {

## Function : Calculates coverage and alignment metrics on BAM files.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $insample_directory      => In sample directory
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outsample_directory     => Out sample directory
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $insample_directory;
    my $job_id_href;
    my $outsample_directory;
    my $parameter_href;
    my $program_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id
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
        insample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$insample_directory,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir
        },
        outsample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outsample_directory,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Language::Java qw{ java_core };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_dead_end };
    use MIP::Program::Alignment::Picardtools
      qw{ picardtools_collectmultiplemetrics };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $DOT => q{.};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $program_mode = $active_parameter_href->{$program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$program_name}{chain};
    my ( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name      => $program_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$sample_id}{gatk_baserecalibration}{file_tag};
    my $outfile_tag =
      $file_info_href->{$sample_id}{gatk_baserecalibration}{file_tag};

    ## Files
    my $infile_prefix  = $merged_infile_prefix . $infile_tag;
    my $outfile_prefix = $merged_infile_prefix . $outfile_tag;

    ## Paths
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ## Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            jobid_chain    => $parameter_href->{gatk_baserecalibration}{chain},
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
        }
    );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $sample_id,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            job_id_href           => $job_id_href,
            log                   => $log,
            process_time          => $time,
            program_directory => catfile( $outaligner_dir, q{coveragereport} ),
            program_name      => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => catfile(
                $insample_directory,
                $infile_prefix . substr( $infile_suffix, 0, 2 ) . $ASTERIX
            ),
            outfile_path => $temp_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## CollectMultipleMetrics
    say {$FILEHANDLE} q{## Collecting multiple metrics on alignment};

    picardtools_collectmultiplemetrics(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $file_path_prefix . $infile_suffix,
            java_jar    => catfile(
                $active_parameter_href->{picardtools_path},
                q{picard.jar}
            ),
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            memory_allocation => q{Xmx4g},
            outfile_path      => $outfile_path_prefix,
            referencefile_path =>
              $active_parameter_href->{human_genome_reference},
            temp_directory => $temp_directory,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    my @outfiles = (
        $outfile_prefix . $DOT . q{alignment_summary_metrics},
        $outfile_prefix . $DOT . q{quality} . $ASTERIX,
        $outfile_prefix . $DOT . q{insert} . $ASTERIX,
    );

  OUTFILE:
    foreach my $outfile (@outfiles) {

        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => catfile( $temp_directory, $outfile ),
                outfile_path => $outsample_directory,
            }
        );
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    if ( $program_mode == 1 ) {

        ## Collect QC metadata info for later use
        my $metrics_outfile_path = catfile( $outsample_directory,
            $outfile_prefix . $DOT . q{alignment_summary_metrics} );
        add_program_outfile_to_sample_info(
            {
                infile           => $merged_infile_prefix,
                path             => $metrics_outfile_path,
                program_name     => q{collectmultiplemetrics},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );
        my $insertsize_outfile_path = catfile( $outsample_directory,
            $outfile_prefix . $DOT . q{insert_size_metrics} );
        add_program_outfile_to_sample_info(
            {
                infile           => $merged_infile_prefix,
                path             => $insertsize_outfile_path,
                program_name     => q{collectmultiplemetricsinsertsize},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );
    }
    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        slurm_submit_job_sample_id_dependency_dead_end(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_id               => $sample_id,
                sbatch_file_name        => $file_path,
            }
        );
    }
    return;
}

1;
