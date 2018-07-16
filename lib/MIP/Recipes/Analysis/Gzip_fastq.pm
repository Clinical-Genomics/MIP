package MIP::Recipes::Analysis::Gzip_fastq;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
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
    our @EXPORT_OK = qw{ analysis_gzip_fastq };

}

## Constants
Readonly my $DOT     => q{.};
Readonly my $NEWLINE => qq{\n};

sub analysis_gzip_fastq {

## Function : Gzips fastq files
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $family_id;

    my $tmpl = {
        active_parameter_href => {
            defined     => 1,
            default     => {},
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
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
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Cluster qw(check_max_core_number);
    use MIP::Cluster qw(update_core_number_to_seq_mode);
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_no_dependency_add_to_sample);
    use MIP::Program::Compression::Gzip qw(gzip);
    use MIP::Script::Setup_script qw(setup_script);

    ## No uncompressed fastq infiles
    return if ( not $file_info_href->{is_file_uncompressed}{$sample_id} );

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $program_mode = $active_parameter_href->{$program_name};

    ## Alias
    my @infiles      = @{ $file_info_href->{$sample_id}{mip_infiles} };
    my $job_id_chain = $parameter_href->{$program_name}{chain};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Adjust according to number of infiles to process
    # One full lane on Hiseq takes approx. 2 h for gzip to process
    $time = $time * scalar @infiles;

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

  INFILE_LANE:
    foreach my $infile ( @{ $infile_lane_prefix_href->{$sample_id} } ) {

        ## Update the number of cores to be used in the analysis according to sequencing mode requirements
        $core_number = update_core_number_to_seq_mode(
            {
                core_number => $core_number,
                sequence_run_type =>
                  $sample_info_href->{sample}{$sample_id}{file}{$infile}
                  {sequence_run_type},
            }
        );
    }

    ## Limit number of cores requested to the maximum number of cores available per node
    $core_number = check_max_core_number(
        {
            core_number_requested => $core_number,
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_name) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $sample_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            program_directory               => lc($program_name),
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
        }
    );

    ## Assign directories
    my $insample_directory = $file_info_href->{$sample_id}{mip_infiles_dir};

    ## Assign suffix
    my $infile_suffix = $parameter_href->{$program_name}{infile_suffix};

    my $process_batches_count = 1;

# Used to print wait at the right times since infiles cannot be used (can be a mixture of .gz and .fast files)
    my $uncompressed_file_counter = 0;

    ## Gzip
    say {$FILEHANDLE} q{## } . $program_name;

  INFILE:
    foreach my $infile (@infiles) {

        ## For files ending with .fastq required since there can be a mixture (also .fastq.gz) within the sample dir
        if ( $infile =~ /$infile_suffix$/sxm ) {

            ## Using only $active_parameter{max_cores_per_node} cores
            if ( $uncompressed_file_counter == $process_batches_count *
                $active_parameter_href->{max_cores_per_node} )
            {

                say {$FILEHANDLE} q{wait}, $NEWLINE;
                $process_batches_count = $process_batches_count + 1;
            }

            ## Perl wrapper for writing gzip recipe to $FILEHANDLE
            gzip(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => catfile( $insample_directory, $infile ),
                }
            );
            say {$FILEHANDLE} q{&};
            $uncompressed_file_counter++;

            ## Add ".gz" to original fastq ending, since this will execute before fastQC and bwa.
            $infile .= $DOT . q{gz};
        }
    }
    print {$FILEHANDLE} $NEWLINE;
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    if ( $program_mode == 1 ) {

        slurm_submit_job_no_dependency_add_to_sample(
            {
                family_id        => $family_id,
                job_id_href      => $job_id_href,
                log              => $log,
                path             => $job_id_chain,
                sample_id        => $sample_id,
                sbatch_file_name => $file_name
            }
        );
    }
    return;
}

1;
