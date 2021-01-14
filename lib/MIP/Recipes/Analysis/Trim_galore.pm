package MIP::Recipes::Analysis::Trim_galore;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use POSIX qw{ floor };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $AMPERSAND $DOT $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_trim_galore };

}

sub analysis_trim_galore {

## Function : Trim reads using Trim galore
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name             => Recipe name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
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
        profile_base_command => {
            default     => q{sbatch},
            store       => \$profile_base_command,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
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

    use MIP::Cluster qw{ update_memory_allocation };
    use MIP::File_info qw{ get_sample_file_attribute };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Program::Trim_galore qw{ trim_galore };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
        }
    );
    my @infile_paths         = @{ $io{in}{file_paths} };
    my @infile_names         = @{ $io{in}{file_names} };
    my @infile_name_prefixes = @{ $io{in}{file_name_prefixes} };

## Construct outfiles
    my $outsample_directory =
      catdir( $active_parameter_href->{outdata_dir}, $sample_id, $recipe_name );

    my @outfile_paths = _construct_trim_galore_outfile_paths(
        {
            file_info_href           => $file_info_href,
            infile_name_prefixes_ref => \@infile_name_prefixes,
            outsample_directory      => $outsample_directory,
            sample_id                => $sample_id,
        }
    );

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id       => $recipe{job_id_chain},
                id             => $sample_id,
                file_info_href => $file_info_href,
                file_paths_ref => \@outfile_paths,
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
            }
        )
    );
    my @outfile_name_prefixes = @{ $io{out}{file_name_prefixes} };
    my $outdata_dir_path      = $io{out}{dir_path};
    my @qc_outfile_paths =
      map { catfile( $outdata_dir_path, $_ . q{_trimming_report.txt} ) } @infile_names;

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    my %file_info_sample = get_sample_file_attribute(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    my $parallel_processes = scalar @{ $file_info_sample{no_direction_infile_prefixes} };
    my ( $process_core_number, $recipe_core_number ) = _get_cores_for_trimgalore(
        {
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
            parallel_processes => $parallel_processes,
        }
    );

    my $memory_allocation = update_memory_allocation(
        {
            node_ram_memory           => $active_parameter_href->{node_ram_memory},
            parallel_processes        => $recipe_core_number,
            process_memory_allocation => $recipe{memory},
        }
    );

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe_core_number,
            directory_id          => $sample_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $memory_allocation,
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
        }
    );

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    # Keep track of paired end files
    my $paired_end_tracker = 0;
    my %qc_files;

  INFILE_PREFIX:
    foreach my $infile_prefix ( @{ $file_info_sample{no_direction_infile_prefixes} } ) {

        my $sequence_run_type = get_sample_file_attribute(
            {
                attribute      => q{sequence_run_type},
                file_info_href => $file_info_href,
                file_name      => $infile_prefix,
                sample_id      => $sample_id,
            }
        );

        ## Get files
        my @fastq_files = $infile_paths[$paired_end_tracker];
        $qc_files{ $infile_names[$paired_end_tracker] } =
          $qc_outfile_paths[$paired_end_tracker];

        my $paired_reads;

        # If second read direction is present
        if ( $sequence_run_type eq q{paired-end} ) {

            $paired_reads = 1;

            # Increment to collect correct read 2 from %infile
            $paired_end_tracker++;
            push @fastq_files, $infile_paths[$paired_end_tracker];
            $qc_files{ $infile_names[$paired_end_tracker] } =
              $qc_outfile_paths[$paired_end_tracker];
        }

        my $stderrfile_path =
          catfile( $recipe_info_path . $DOT . $infile_prefix . $DOT . q{stderr.txt} );

        ## Trim galore
        trim_galore(
            {
                cores            => $process_core_number,
                filehandle       => $filehandle,
                infile_paths_ref => \@fastq_files,
                outdir_path      => $outsample_directory,
                paired_reads     => $paired_reads,
                stderrfile_path  => $stderrfile_path,
            }
        );
        say {$filehandle} $AMPERSAND . $NEWLINE;

        ## Increment paired end tracker
        $paired_end_tracker++;

    }
    say {$filehandle} q{wait};

    ## Close filehandleS
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

        ## Outfiles
        my $outfile_path        = $outfile_paths[0];
        my $outfile_name_prefix = $outfile_name_prefixes[0];

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => $outfile_name_prefix,
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

      FASTQ_INFILE:
        foreach my $fastq_infile ( keys %qc_files ) {

            ## Collect QC metadata info for later use
            set_recipe_outfile_in_sample_info(
                {
                    infile           => $fastq_infile,
                    path             => $qc_files{$fastq_infile},
                    recipe_name      => q{trim_galore_stats},
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );
        }

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_sample},
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
                log                               => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_id          => $sample_id,
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }

    return 1;
}

sub _construct_trim_galore_outfile_paths {

## Function : Construct outfile paths for Trim galore
## Returns  : @outfile_paths
## Arguments: $file_info_href           => File info hash {REF}
##          : $infile_name_prefixes_ref => Array with infile name prefixes {REF}
##          : $outsample_directory      => Outsample directory
##          : $sample_id                => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $infile_name_prefixes_ref;
    my $outsample_directory;
    my $sample_id;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        infile_name_prefixes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_name_prefixes_ref,
            strict_type => 1,
        },
        outsample_directory => {
            required    => 1,
            store       => \$outsample_directory,
            strict_type => 1,
        },
        sample_id => {
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info qw{ get_sample_file_attribute };

    my @outfile_paths;

    my $paired_end_tracker = 0;

    my %file_info_sample = get_sample_file_attribute(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );
  INFILE_PREFIX:
    foreach my $infile_prefix ( @{ $file_info_sample{no_direction_infile_prefixes} } ) {

        my $sequence_run_type = get_sample_file_attribute(
            {
                attribute      => q{sequence_run_type},
                file_info_href => $file_info_href,
                file_name      => $infile_prefix,
                sample_id      => $sample_id,
            }
        );

        my $outfile_path_prefix =
          catfile( $outsample_directory, ${$infile_name_prefixes_ref}[$paired_end_tracker] );

        ## The suffixes differs depending on whether the reads are paired or not
        if ( $sequence_run_type eq q{paired-end} ) {
            ## Add read 1
            push @outfile_paths, $outfile_path_prefix . q{_val_1.fq.gz};

            ## Add read 2
            $paired_end_tracker++;

            $outfile_path_prefix =
              catfile( $outsample_directory, ${$infile_name_prefixes_ref}[$paired_end_tracker] );
            push @outfile_paths, $outfile_path_prefix . q{_val_2.fq.gz};

        }
        else {
            push @outfile_paths, $outfile_path_prefix . q{_trimmed.fq.gz};
        }
        $paired_end_tracker++;
    }
    return @outfile_paths;
}

sub _get_cores_for_trimgalore {

## Function : Calculate nr of cores to be supplied to each Trim galore process and the recipe
## Returns  : $core_argument, $recipe_core_number
## Arguments: $max_cores_per_node => Cores available
##          : $parallel_processes => Run fastqc after trimming

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $max_cores_per_node;
    my $parallel_processes;

    my $tmpl = {
        max_cores_per_node => {
            allow       => [qr/\A \d+ \z/xms],
            required    => 1,
            store       => \$max_cores_per_node,
            strict_type => 1,
        },
        parallel_processes => {
            allow       => [qr/\A \d+ \z/xms],
            required    => 1,
            store       => \$parallel_processes,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    Readonly my $THREE => 3;

    ## Currently (Trim galore v0.6.5) the way to calculate the core argument to trim galore:
    ## Always three cores for overhead (1 for trim galore and 2 for cutadapt)
    ## the rest are splitted between the three processe (read, write and cutadapt).
    my $core_argument = floor( ( $max_cores_per_node / $parallel_processes - $THREE ) / $THREE );
    my $recipe_core_number = ( $core_argument * $THREE + $THREE ) * $parallel_processes;

    ## Only supply core argument if more than 1
    $core_argument      = $core_argument > 1 ? $core_argument      : 0;
    $recipe_core_number = $core_argument     ? $recipe_core_number : $parallel_processes;

    return ( $core_argument, $recipe_core_number );
}

1;
