package MIP::Recipes::Analysis::Fastqc;

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
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(analysis_fastqc);

}

##Constants
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_fastqc {

##analysis_fastqc

##Function : Raw sequence quality analysis using FASTQC.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $infiles_ref, $infile_lane_prefix_href, $job_id_href, $insample_directory, $outsample_directory, $sample_id, $program_name, $temp_directory
##         : $parameter_href          => Parameter hash {REF}
##         : $active_parameter_href   => Active parameters for this analysis hash {REF}
##         : $sample_info_href        => Info on samples and family hash {REF}
##         : $infiles_ref             => Infiles {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href             => Job id hash {REF}
##         : $insample_directory      => In sample directory
##         : $outsample_directory     => Out sample directory
##         : $sample_id               => Sample id
##         : $program_name            => Program name
##         : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Default(s)
    my $temp_directory;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $infiles_ref;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $insample_directory;
    my $outsample_directory;
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
        infiles_ref => {
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infiles_ref
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
        outsample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outsample_directory
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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    use MIP::Check::Cluster qw{check_max_core_number};
    use MIP::Cluster qw{update_core_number_to_seq_mode};
    use MIP::Script::Setup_script qw{setup_script};
    use MIP::IO::Files qw{migrate_files};
    use MIP::Processmanagement::Processes qw{print_wait};
    use MIP::Gnu::Coreutils qw{gnu_cp};
    use MIP::Program::Qc::Fastqc qw{fastqc};
    use MIP::QC::Record qw{add_program_outfile_to_sample_info};
    use MIP::Processmanagement::Slurm_processes
      qw{slurm_submit_job_no_dependency_dead_end};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

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
            program_directory     => $program_name,
            core_number           => $core_number,
            process_time =>
              $active_parameter_href->{module_time}{$mip_program_name},
            temp_directory => $temp_directory,
        }
    );

    ## Assign suffix
    my $infile_suffix = $parameter_href->{$mip_program_name}{infile_suffix};

    ## Copies files from source to destination
    migrate_files(
        {
            infiles_ref  => \@{$infiles_ref},
            outfile_path => $temp_directory,
            FILEHANDLE   => $FILEHANDLE,
            indirectory  => $insample_directory,
            core_number  => $core_number,
        }
    );

    say {$FILEHANDLE} q{## } . $program_name;

    my $process_batches_count = 1;

    while ( my ( $index, $infile ) = each @{$infiles_ref} ) {

        $process_batches_count = print_wait(
            {
                process_counter       => $index,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                FILEHANDLE            => $FILEHANDLE,
            }
        );

        ## Removes ".file_ending" in filename.FILENDING(.gz)
        my $file_at_lane_level =
          fileparse( $infile, qr/$infile_suffix|$infile_suffix[.]gz/sxm );

        fastqc(
            {
                infile_path       => catfile( $temp_directory, $infile ),
                outdirectory_path => $temp_directory,
                extract           => 1,
                FILEHANDLE        => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} q{&}, $NEWLINE;

        ## Collect QC metadata info for active program for later use
        if ( $mip_program_mode == 1 ) {

            my $qc_fastqc_outdirectory =
              catdir( $outsample_directory,
                $file_at_lane_level . $UNDERSCORE . $program_name );
            add_program_outfile_to_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    sample_id        => $sample_id,
                    program_name     => $program_name,
                    infile           => $infile,
                    outdirectory     => $qc_fastqc_outdirectory,
                    outfile          => q{fastqc_data.txt},
                }
            );
        }
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Copies files from temporary folder to source.
    $process_batches_count = 1;
    while ( my ( $index, $infile ) = each @{$infiles_ref} ) {

        $process_batches_count = print_wait(
            {
                process_counter       => $index,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                FILEHANDLE            => $FILEHANDLE,
            }
        );

        ## Removes ".file_ending" in filename.FILENDING(.gz)
        my $file_at_lane_level =
          fileparse( $infile, qr/$infile_suffix|$infile_suffix[.]gz/sxm );

        my $infile_path = catfile( $temp_directory,
            $file_at_lane_level . $UNDERSCORE . $program_name );
        gnu_cp(
            {
                FILEHANDLE   => $FILEHANDLE,
                recursive    => 1,
                infile_path  => $infile_path,
                outfile_path => $outsample_directory,
            }
        );
        say {$FILEHANDLE} q{&}, $NEWLINE;
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {

        slurm_submit_job_no_dependency_dead_end(
            {
                job_id_href      => $job_id_href,
                sbatch_file_name => $file_name,
                log              => $log,
            }
        );
    }
    return;
}

1;
