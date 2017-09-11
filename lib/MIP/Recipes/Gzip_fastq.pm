package MIP::Recipes::Gzip_fastq;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use autodie qw{:all};
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Basename qw(fileparse);
use File::Spec::Functions qw(catdir catfile);

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{Exporter};

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(analysis_gzip_fastq);

}

##Constants
Readonly my $DOT     => q{.};
Readonly my $NEWLINE => qq{\n};

sub analysis_gzip_fastq {

##analysis_gzip_fastq

##Function : Gzips fastq files
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $infile_href, $infile_lane_prefix_href, $job_id_href, $insample_directory, $sample_id, $program_name. $family_id
##         : $parameter_href          => Parameter hash {REF}
##         : $active_parameter_href   => Active parameters for this analysis hash {REF}
##         : $sample_info_href        => Info on samples and family hash {REF}
##         : $infile_href             => Infiles hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href             => Job id hash {REF}
##         : $insample_directory      => In sample directory
##         : $sample_id               => Sample id
##         : $program_name            => Program name
##         : $family_id               => Family id

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $infile_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $insample_directory;
    my $sample_id;
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
        infile_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_href
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
        insample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$insample_directory
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw(update_core_number_to_seq_mode);
    use MIP::Check::Cluster qw(check_max_core_number);
    use MIP::Script::Setup_script qw(setup_script);
    use MIP::Program::Compression::Gzip qw(gzip);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_no_dependency_add_to_sample);

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};

    ## Adjust according to number of infiles to process
    # One full lane on Hiseq takes approx. 2 h for gzip to process
    my $time =
      $active_parameter_href->{module_time}{$mip_program_name} *
      scalar @{ $infile_href->{$sample_id} };
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

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
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
            core_number_requested => $core_number,
        }
    );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_name) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $sample_id,
            program_name          => $program_name,
            program_directory     => lc($program_name),
            core_number           => $core_number,
            process_time          => $time,
        }
    );

    ## Assign suffix
    my $infile_suffix = $parameter_href->{$mip_program_name}{infile_suffix};

    my $process_batches_count = 1;

# Used to print wait at the right times since infiles cannot be used (can be a mixture of .gz and .fast files)
    my $uncompressed_file_counter = 0;

    ## Gzip
    say {$FILEHANDLE} q{## } . $program_name;

  INFILE:
    foreach my $infile ( @{ $infile_href->{$sample_id} } ) {

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
                    infile_path => catfile( $insample_directory, $infile ),
                    FILEHANDLE  => $FILEHANDLE,
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

    if ( $mip_program_mode == 1 ) {

        slurm_submit_job_no_dependency_add_to_sample(
            {
                job_id_href      => $job_id_href,
                family_id        => $family_id,
                sample_id        => $sample_id,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_name
            }
        );
    }
    return;
}

1;
