package MIP::Recipes::Analysis::Picardtools_mergesamfiles;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ fileparse };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use POSIX;
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ %ANALYSIS $ASTERISK $DOT $EMPTY_STR $LOG_NAME $NEWLINE $SEMICOLON $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.18;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_picardtools_mergesamfiles };

}

## Constants
Readonly my $JAVA_GUEST_OS_MEMORY => $ANALYSIS{JAVA_GUEST_OS_MEMORY};

sub analysis_picardtools_mergesamfiles {

## Function : Merges all bam files using Picardtools mergesamfiles within each sampleid and files generated previously (option if provided with '-picardtools_mergesamfiles_previous_bams'). The merged files have to be sorted before attempting to merge.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
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
    my $infile_lane_prefix_href;
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

    use MIP::Cluster qw{ get_parallel_processes update_memory_allocation };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Picardtools
      qw{ picardtools_gatherbamfiles picardtools_mergesamfiles };
    use MIP::Program::Alignment::Sambamba qw{ split_and_index_aligment_file };
    use MIP::Program::Samtools qw{ samtools_index };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_merged_infile_prefix };

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
    my $infile_suffix        = $io{in}{file_suffix};
    my @infile_name_prefixes = @{ $io{in}{file_name_prefixes} };
    my @infile_paths         = @{ $io{in}{file_paths} };

    my %rec_atr = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $job_id_chain            = $rec_atr{chain};
    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my $recipe_mode             = $active_parameter_href->{$recipe_name};
    my $xargs_file_path_prefix;
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number       = $recipe_resource{core_number};
    my $memory_allocation = $recipe_resource{memory};

    ## Assign suffix
    my $outfile_suffix = $rec_atr{outfile_suffix};

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
                chain_id       => $job_id_chain,
                id             => $sample_id,
                file_info_href => $file_info_href,
                file_paths_ref => \@outfile_paths,
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
            }
        )
    );

    my $outdir_path = $io{out}{dir_path};
    @outfile_paths = @{ $io{out}{file_paths} };
    my %outfile_path        = %{ $io{out}{file_path_href} };
    my $outfile_path_prefix = $io{out}{file_path_prefix};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle      = IO::Handle->new();
    my $xargsfilehandle = IO::Handle->new();

    ## Get recipe memory allocation
    Readonly my $JAVA_MEMORY_ALLOCATION => 4;
    my $process_memory_allocation = $JAVA_MEMORY_ALLOCATION + $JAVA_GUEST_OS_MEMORY;

    ## Modify memory allocation according to which action is taken in recipe
    if ( scalar @infile_paths > 1 ) {

        # Get recipe memory allocation
        $memory_allocation = update_memory_allocation(
            {
                node_ram_memory           => $active_parameter_href->{node_ram_memory},
                parallel_processes        => $core_number,
                process_memory_allocation => $process_memory_allocation,
            }
        );
    }

    # Constrain parallelization to match available memory
    my $parallel_processes = get_parallel_processes(
        {
            process_memory_allocation => $process_memory_allocation,
            recipe_memory_allocation  => $memory_allocation,
            core_number               => $core_number,
        }
    );

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $sample_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            log                             => $log,
            memory_allocation               => $memory_allocation,
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
            temp_directory                  => $temp_directory,
        }
    );

    ## Set helper value for finding merged_infiles downstream
    my $merged_infile_prefix =
      $sample_id . $UNDERSCORE . q{lanes} . $UNDERSCORE . $lanes_id;
    set_merged_infile_prefix(
        {
            file_info_href       => $file_info_href,
            merged_infile_prefix => $merged_infile_prefix,
            sample_id            => $sample_id,
        }
    );

    ### SHELL:

  INFILE:
    while ( my ( $infile_index, $infile_path ) = each @infile_paths ) {

        ## Split BAMs
        say {$filehandle} q{## Split alignment files per contig};
        ($xargs_file_counter) = split_and_index_aligment_file(
            {
                contigs_ref         => \@{ $file_info_href->{bam_contigs_size_ordered} },
                core_number         => $core_number,
                filehandle          => $filehandle,
                file_path           => $recipe_file_path,
                infile_path         => $infile_path,
                outfile_path_prefix => $outdir_path
                  . $infile_name_prefixes[$infile_index],
                output_format      => substr( $infile_suffix, 1 ),
                recipe_info_path   => $recipe_info_path,
                xargsfilehandle    => $xargsfilehandle,
                xargs_file_counter => $xargs_file_counter,
            }
        );
    }

    ## More than one file - we have something to merge
    if ( scalar @infile_paths > 1 ) {

        ## picardtools_mergesamfiles
        say {$filehandle} q{## Merging alignment files};

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                core_number        => $parallel_processes,
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
            # Assemble infile paths by adding directory and file suffixes
            my @merge_infile_paths =
              map { $outdir_path . $_ . $DOT . $contig . $infile_suffix }
              @infile_name_prefixes;
            my $stderrfile_path =
              $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};

            picardtools_mergesamfiles(
                {
                    create_index     => q{true},
                    filehandle       => $xargsfilehandle,
                    infile_paths_ref => \@merge_infile_paths,
                    java_jar         => catfile(
                        $active_parameter_href->{picardtools_path}, q{picard.jar}
                    ),
                    java_use_large_pages =>
                      $active_parameter_href->{java_use_large_pages},
                    memory_allocation  => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
                    outfile_path       => $outfile_path{$contig},
                    referencefile_path => $referencefile_path,
                    stderrfile_path    => $stderrfile_path,
                    temp_directory     => $temp_directory,
                    threading          => q{true},
                }
            );
            say {$xargsfilehandle} $NEWLINE;
        }
    }
    else {
        ## Only 1 infile - rename sample and index instead of merge to streamline handling of filenames downstream

        ## Rename samples
        say {$filehandle}
q{## Renaming sample instead of merge to streamline handling of filenames downstream};

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

          INFILES:
            foreach my $infile_name_prefix (@infile_name_prefixes) {

                ## Get parameters
                my $gnu_infile_path =
                  $outdir_path . $infile_name_prefix . $DOT . $contig . $infile_suffix;
                my $gnu_outfile_path = $outfile_path{$contig};

                ## Rename
                gnu_mv(
                    {
                        filehandle   => $xargsfilehandle,
                        infile_path  => $gnu_infile_path,
                        outfile_path => $gnu_outfile_path,
                    }
                );
                print {$xargsfilehandle} $SEMICOLON . $SPACE;

                ## Index
                samtools_index(
                    {
                        bai_format  => 1,
                        filehandle  => $xargsfilehandle,
                        infile_path => $gnu_outfile_path,
                    }
                );
            }
            say {$xargsfilehandle} $NEWLINE;
        }
    }

    ## Gather BAM files
    say {$filehandle} q{## Gather BAM files};

    ## Assemble infile paths in contig order and not per size
    my @gather_infile_paths =
      map { $outfile_path{$_} } @{ $file_info_href->{bam_contigs} };

    picardtools_gatherbamfiles(
        {
            create_index     => q{true},
            filehandle       => $filehandle,
            infile_paths_ref => \@gather_infile_paths,
            java_jar =>
              catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
            java_use_large_pages => $active_parameter_href->{java_use_large_pages},
            memory_allocation    => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
            outfile_path         => $outfile_path_prefix . $outfile_suffix,
            referencefile_path   => $referencefile_path,
            temp_directory       => $temp_directory,
        }
    );
    say {$filehandle} $NEWLINE;

    close $xargsfilehandle;
    close $filehandle;

    if ( $recipe_mode == 1 ) {

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
                base_command            => $profile_base_command,
                case_id                 => $case_id,
                dependency_method       => q{sample_to_sample},
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_chain            => $job_id_chain,
                job_id_href             => $job_id_href,
                log                     => $log,
                recipe_file_path        => $recipe_file_path,
                sample_id               => $sample_id,
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

1;
