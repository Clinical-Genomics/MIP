package MIP::Recipes::Analysis::Samtools_merge;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ fileparse };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use POSIX;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ %ANALYSIS $ASTERISK $DOT $EMPTY_STR $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_samtools_merge analysis_samtools_merge_panel };

}

sub analysis_samtools_merge {

## Function : Merges all bam files using Samtools merge within each sampleid and files generated previously.
##          : The merged files have to be sorted before attempting to merge.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $referencefile_path      => Human genome reference file path
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;
    my $sample_id;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
    my $referencefile_path;
    my $temp_directory;
    my $xargs_file_counter;

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
        referencefile_path => {
            default     => $arg_href->{active_parameter_href}{human_genome_reference},
            store       => \$referencefile_path,
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
        xargs_file_counter => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info qw{ get_io_files set_merged_infile_prefix };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Samtools qw{ samtools_merge samtools_view };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
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
    my @infile_paths = @{ $io{in}{file_paths} };

    my $xargs_file_path_prefix;
    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number       = $recipe{core_number};
    my $memory_allocation = $recipe{memory};

    ## Assign suffix
    my $outfile_suffix = $recipe{outfile_suffix};

    ## Extract lanes
    my $lanes_id = join $EMPTY_STR, @{ $file_info_href->{$sample_id}{lanes} };

    ## Outpaths
    my $outsample_directory =
      catdir( $active_parameter_href->{outdata_dir}, $sample_id, $recipe_name );
    my $outfile_tag =
      $file_info_href->{$sample_id}{$recipe_name}{file_tag};
    my @outfile_paths =
      map {
        catdir( $outsample_directory,
                $sample_id
              . $UNDERSCORE
              . q{lanes}
              . $UNDERSCORE
              . $lanes_id
              . $outfile_tag
              . $DOT
              . $_
              . $outfile_suffix )
      } @{ $file_info_href->{bam_contigs_size_ordered} };

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

    @outfile_paths = @{ $io{out}{file_paths} };
    my %outfile_path = %{ $io{out}{file_path_href} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle      = IO::Handle->new();
    my $xargsfilehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $sample_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $memory_allocation,
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            temp_directory        => $temp_directory,
        }
    );

    ## Set helper value for finding merged_infiles downstream
    my $merged_infile_prefix = $sample_id . $UNDERSCORE . q{lanes} . $UNDERSCORE . $lanes_id;
    set_merged_infile_prefix(
        {
            file_info_href       => $file_info_href,
            merged_infile_prefix => $merged_infile_prefix,
            sample_id            => $sample_id,
        }
    );

    ### SHELL:

    ## More than one file - we have something to merge
    if ( scalar @infile_paths > 1 ) {

        ## Samtools_merge
        say {$filehandle} q{## Merging alignment files};

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                core_number        => $core_number,
                filehandle         => $filehandle,
                file_path          => $recipe_file_path,
                recipe_info_path   => $recipe_info_path,
                xargsfilehandle    => $xargsfilehandle,
                xargs_file_counter => $xargs_file_counter,
            }
        );

      CONTIG:
        foreach my $contig ( @{ $file_info_href->{bam_contigs_size_ordered} } ) {

            ## Get parameters
            my $stderrfile_path = $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};

            samtools_merge(
                {
                    filehandle         => $xargsfilehandle,
                    force              => 1,
                    infile_paths_ref   => \@infile_paths,
                    outfile_path       => $outfile_path{$contig},
                    output_format      => q{bam},
                    referencefile_path => $referencefile_path,
                    region             => $contig,
                    stderrfile_path    => $stderrfile_path,
                    thread_number      => 2,
                    write_index        => 1,
                }
            );
            say {$xargsfilehandle} $NEWLINE;
        }
    }
    else {
        ## Only 1 infile - rename sample and index instead of merge to streamline handling of filenames downstream

        ## Rename samples
        say {$filehandle}
          q{## Split file into contigs instead of merge to streamline handling of files downstream};

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                core_number        => $core_number,
                filehandle         => $filehandle,
                file_path          => $recipe_file_path,
                recipe_info_path   => $recipe_info_path,
                xargsfilehandle    => $xargsfilehandle,
                xargs_file_counter => $xargs_file_counter,
            }
        );

      CONTIG:
        foreach my $contig ( @{ $file_info_href->{bam_contigs_size_ordered} } ) {

            samtools_view(
                {
                    filehandle    => $xargsfilehandle,
                    infile_path   => $infile_paths[0],
                    outfile_path  => $outfile_path{$contig},
                    output_format => q{bam},
                    regions_ref   => [$contig],
                    write_index   => 1,
                }
            );
            say {$xargsfilehandle} $NEWLINE;
        }
    }

    close $xargsfilehandle;
    close $filehandle;

    if ( $recipe{mode} == 1 ) {

        my $qc_outfile_path = $outfile_paths[0];
        set_recipe_outfile_in_sample_info(
            {
                infile           => $merged_infile_prefix,
                path             => $qc_outfile_path,
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

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

sub analysis_samtools_merge_panel {

## Function : Merges all bam files using Samtools merge within each sample id and files generated previously.
##          : The merged files have to be sorted before attempting to merge.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;
    my $sample_id;

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

    use MIP::File_info qw{ get_io_files set_merged_infile_prefix };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Samtools qw{ samtools_merge samtools_view };
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
    my @infile_paths = @{ $io{in}{file_paths} };

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number       = $recipe{core_number};
    my $memory_allocation = $recipe{memory};

    ## Assign suffix
    my $outfile_suffix = $recipe{outfile_suffix};

    ## Extract lanes
    my $lanes_id = join $EMPTY_STR, @{ $file_info_href->{$sample_id}{lanes} };

    ## Outpaths
    my $outsample_directory =
      catdir( $active_parameter_href->{outdata_dir}, $sample_id, $recipe_name );
    my $outfile_tag =
      $file_info_href->{$sample_id}{$recipe_name}{file_tag};
    my $outfile_path = catfile( $outsample_directory,
            $sample_id
          . $UNDERSCORE
          . q{lanes}
          . $UNDERSCORE
          . $lanes_id
          . $outfile_tag
          . $outfile_suffix );

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id       => $recipe{job_id_chain},
                id             => $sample_id,
                file_info_href => $file_info_href,
                file_paths_ref => [$outfile_path],
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
            }
        )
    );
    $outfile_path = $io{out}{file_path};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $sample_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $memory_allocation,
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
        }
    );

    ## Set helper value for finding merged_infiles downstream
    my $merged_infile_prefix = $sample_id . $UNDERSCORE . q{lanes} . $UNDERSCORE . $lanes_id;
    set_merged_infile_prefix(
        {
            file_info_href       => $file_info_href,
            merged_infile_prefix => $merged_infile_prefix,
            sample_id            => $sample_id,
        }
    );

    ### SHELL:

    ## More than one file - we have something to merge
    if ( scalar @infile_paths > 1 ) {

        ## Samtools_merge
        say {$filehandle} q{## Merging alignment files};
        samtools_merge(
            {
                filehandle         => $filehandle,
                force              => 1,
                infile_paths_ref   => \@infile_paths,
                outfile_path       => $outfile_path,
                output_format      => q{bam},
                referencefile_path => ${active_parameter_href}->{human_genome_reference},
                thread_number      => 2,
                write_index        => 1,
            }
        );
        say {$filehandle} $NEWLINE;
    }
    else {
        ## Only 1 infile - rename sample and index instead of merge to streamline handling of filenames downstream

        say {$filehandle} q{## Only one infile, rename to streamline handling of files downstream};
        samtools_view(
            {
                filehandle    => $filehandle,
                infile_path   => $infile_paths[0],
                outfile_path  => $outfile_path,
                output_format => q{bam},
                write_index   => 1,
            }
        );
        say {$filehandle} $NEWLINE;
    }

    close $filehandle;

    if ( $recipe{mode} == 1 ) {

        set_recipe_outfile_in_sample_info(
            {
                infile           => $merged_infile_prefix,
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

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

1;
