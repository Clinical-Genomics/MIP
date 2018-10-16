package MIP::Recipes::Analysis::Fastqc;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw(fileparse);
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{:all};
use Readonly;

BEGIN {

    require Exporter;
    use base qw{Exporter};

    # Set the version for version checking
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_fastqc };

}

## Constants
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_fastqc {

## Function : Raw sequence quality analysis using FASTQC.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
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
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    use MIP::Check::Cluster qw{ check_max_core_number };
    use MIP::Cluster qw{ update_core_number_to_seq_mode };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_module_parameters get_program_attributes };
    use MIP::Gnu::Coreutils qw{ gnu_mkdir };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ print_wait };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_no_dependency_dead_end };
    use MIP::Program::Qc::Fastqc qw{ fastqc };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            program_name   => $program_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my @infile_paths         = @{ $io{in}{file_paths} };
    my @infile_name_prefixes = @{ $io{in}{file_name_prefixes} };

    my $job_id_chain = get_program_attributes(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
            attribute      => q{chain},
        }
    );
    my $program_mode = $active_parameter_href->{$program_name};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Outpaths
    my $outsample_directory =
      catdir( $active_parameter_href->{outdata_dir}, $sample_id,
        $program_name );
    my @outfile_paths =
      map { catdir( $outsample_directory, $_ . $UNDERSCORE . $program_name, q{fastqc_data.txt} ) }
      @infile_name_prefixes;

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id       => $job_id_chain,
                id             => $sample_id,
                file_info_href => $file_info_href,
                file_paths_ref => \@outfile_paths,
                parameter_href => $parameter_href,
                program_name   => $program_name,
                temp_directory => $temp_directory,
            }
        )
    );

    my $outdir_path_prefix    = $io{out}{dir_path_prefix};
    my @outfile_name_prefixes = @{ $io{out}{file_name_prefixes} };
    @outfile_paths = @{ $io{out}{file_paths} };

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
            program_directory               => $program_name,
            program_name                    => $program_name,
            sleep                           => 1,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    say {$FILEHANDLE} q{## Create output dir};
    gnu_mkdir(
        {
            FILEHANDLE       => $FILEHANDLE,
            indirectory_path => $outsample_directory,
            parents          => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## } . $program_name;

    my $process_batches_count = 1;

    while ( my ( $index, $infile_path ) = each @infile_paths ) {

        $process_batches_count = print_wait(
            {
                FILEHANDLE            => $FILEHANDLE,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                process_counter       => $index,
            }
        );

        fastqc(
            {
                extract           => 1,
                FILEHANDLE        => $FILEHANDLE,
                infile_path       => $infile_path,
                outdirectory_path => $outsample_directory,
            }
        );
        say {$FILEHANDLE} q{&}, $NEWLINE;

        ## Collect QC metadata info for active program for later use
        if ( $program_mode == 1 ) {

            add_program_outfile_to_sample_info(
                {
                    infile           => $outfile_name_prefixes[$index],
                    path             => $outfile_paths[$index],
                    directory        => $outdir_path_prefix,
                    program_name     => $program_name,
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );
        }
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        slurm_submit_job_no_dependency_dead_end(
            {
                job_id_href      => $job_id_href,
                log              => $log,
                sbatch_file_name => $file_name,
            }
        );
    }
    return;
}

1;
