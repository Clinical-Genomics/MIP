package MIP::Recipes::Analysis::Markduplicates;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
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
use MIP::Constants qw{ %ANALYSIS $ASTERISK $DOT $LOG_NAME $NEWLINE $SPACE $SEMICOLON $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.28;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ analysis_markduplicates analysis_markduplicates_panel analysis_markduplicates_rna };

}

## Constants
Readonly my $JAVA_MEMORY_ALLOCATION      => 6;
Readonly my $JAVA_MEMORY_RECIPE_ADDITION => 2;
Readonly my $JAVA_GUEST_OS_MEMORY => $ANALYSIS{JAVA_GUEST_OS_MEMORY} + $JAVA_MEMORY_RECIPE_ADDITION;

sub analysis_markduplicates {

## Function : Mark duplicated reads using Picardtools markduplicates in files generated from alignment (sorted, merged).
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
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
    my $profile_base_command;
    my $recipe_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
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

    use MIP::Get::File qw{ get_merged_infile_prefix get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Program::Gnu::Coreutils qw{ gnu_cat };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Picardtools qw{ picardtools_markduplicates };
    use MIP::Program::Samtools qw{ samtools_flagstat samtools_view };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
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
    my %infile_path = %{ $io{in}{file_path_href} };

    my %rec_atr = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $job_id_chain       = $rec_atr{chain};
    my $recipe_mode        = $active_parameter_href->{$recipe_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $xargs_file_path_prefix;
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number = $recipe_resource{core_number};

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Outpaths
    ## Assign suffix
    my $outfile_suffix = $rec_atr{outfile_suffix};
    my $outsample_directory =
      catdir( $active_parameter_href->{outdata_dir}, $sample_id, $recipe_name );
    my $outfile_tag =
      $file_info_href->{$sample_id}{$recipe_name}{file_tag};
    my @outfile_paths =
      map {
        catdir( $outsample_directory,
            $merged_infile_prefix . $outfile_tag . $DOT . $_ . $outfile_suffix )
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
    my $outfile_name_prefix = $io{out}{file_name_prefix};
    my %outfile_path        = %{ $io{out}{file_path_href} };
    my $outfile_path_prefix = $io{out}{file_path_prefix};

    ## Update memory and parallel processes for markduplicates
    my ( $recipe_memory, $contig_memory_href, $parallel_processes ) = _get_markdup_resources(
        {
            active_contigs_ref  => $file_info_href->{bam_contigs},
            node_memory         => $active_parameter_href->{node_ram_memory},
            primary_contigs_ref => $file_info_href->{primary_contigs},
            recipe_core_number  => $core_number,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle      = IO::Handle->new();
    my $xargsfilehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $parallel_processes,
            directory_id                    => $sample_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            log                             => $log,
            memory_allocation               => $recipe_memory,
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    ## Marking Duplicates
    say {$filehandle} q{## Marking Duplicates};

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

        my $stderrfile_path = $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        my $metrics_file    = $outfile_path_prefix . $DOT . $contig . $DOT . q{metric};

        picardtools_markduplicates(
            {
                create_index     => q{true},
                filehandle       => $xargsfilehandle,
                infile_paths_ref => [ $infile_path{$contig} ],
                java_jar => catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
                java_use_large_pages       => $active_parameter_href->{java_use_large_pages},
                memory_allocation          => q{Xmx} . $contig_memory_href->{$contig} . q{g},
                metrics_file               => $metrics_file,
                optical_duplicate_distance =>
                  $active_parameter_href->{markduplicates_picardtools_opt_dup_dist},
                outfile_path       => $outfile_path{$contig},
                referencefile_path => $referencefile_path,
                stderrfile_path    => $stderrfile_path,
                temp_directory     => $temp_directory,
            }
        );
        print {$xargsfilehandle} $SEMICOLON . $SPACE;

        ## Process BAM with samtools flagstat to produce metric file for downstream analysis
        samtools_flagstat(
            {
                filehandle      => $xargsfilehandle,
                infile_path     => $outfile_path{$contig},
                stdoutfile_path => $outfile_path_prefix . $DOT . $contig . $UNDERSCORE . q{metric},
                stderrfile_path_append => $stderrfile_path,
            }
        );
        say {$xargsfilehandle} $NEWLINE;
    }

    ## Concatenate all metric files
    gnu_cat(
        {
            filehandle       => $filehandle,
            infile_paths_ref =>
              [ $outfile_path_prefix . $DOT . $ASTERISK . $UNDERSCORE . q{metric} ],
            stdoutfile_path => $outfile_path_prefix . $UNDERSCORE . q{metric_all},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Collect duplicate reads and reads mapped across all metric contig files. Calculate fraction duplicates.
    ## Write it to stdout.
    _calculate_fraction_duplicates_for_all_metric_files(
        {
            filehandle          => $filehandle,
            outfile_path_prefix => $outfile_path_prefix,
        }
    );

    ## Close filehandles
    close $xargsfilehandle;
    close $filehandle;

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => $outfile_name_prefix,
                path             => catfile( $outfile_path_prefix . $UNDERSCORE . q{metric} ),
                recipe_name      => q{markduplicates},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_sample},
                job_id_chain                      => $job_id_chain,
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

sub analysis_markduplicates_panel {

## Function : Mark duplicated reads using Picardtools markduplicates in files generated from alignment (sorted, merged).
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
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
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Picardtools qw{ picardtools_markduplicates };
    use MIP::Program::Samtools qw{ samtools_flagstat samtools_view };
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
    my $infile_path        = $io{in}{file_path};
    my $infile_name_prefix = $io{in}{file_name_prefix};

    my %rec_atr = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $job_id_chain       = $rec_atr{chain};
    my $recipe_mode        = $active_parameter_href->{$recipe_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my %recipe_resource    = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number       = $recipe_resource{core_number};
    my $memory_allocation = $recipe_resource{memory};

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => [$infile_name_prefix],
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );

    #my $outfile_name_prefix = $io{out}{file_name_prefix};
    my $outfile_path        = $io{out}{file_path};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_name_prefix = $io{out}{file_name_prefix};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Update recipe memory allocation for picard
    ## Variables used downstream of if statment
    my $process_memory_allocation = $JAVA_MEMORY_ALLOCATION + $JAVA_GUEST_OS_MEMORY;

    # Get recipe memory allocation
    $memory_allocation = update_memory_allocation(
        {
            node_ram_memory           => $active_parameter_href->{node_ram_memory},
            parallel_processes        => $core_number,
            process_memory_allocation => $process_memory_allocation,
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
            temp_directory                  => ${active_parameter_href}->{temp_directory},
        }
    );

    ### SHELL:

    ## Marking Duplicates
    say {$filehandle} q{## Marking Duplicates};

    my $metrics_file = $outfile_path_prefix . $DOT . q{metric};
    picardtools_markduplicates(
        {
            create_index     => q{true},
            filehandle       => $filehandle,
            infile_paths_ref => [$infile_path],
            java_jar => catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
            java_use_large_pages       => $active_parameter_href->{java_use_large_pages},
            memory_allocation          => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
            metrics_file               => $metrics_file,
            optical_duplicate_distance =>
              $active_parameter_href->{markduplicates_picardtools_opt_dup_dist},
            outfile_path       => $outfile_path,
            referencefile_path => $referencefile_path,
            temp_directory     => ${active_parameter_href}->{temp_directory},
        }
    );
    print {$filehandle} $SEMICOLON . $SPACE;

    ## Process BAM with samtools flagstat to produce metric file for downstream analysis
    samtools_flagstat(
        {
            filehandle      => $filehandle,
            infile_path     => $outfile_path,
            stdoutfile_path => $outfile_path_prefix . $UNDERSCORE . q{metric_all},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Collect duplicate reads and reads mapped across all metric contig files. Calculate fraction duplicates.
    ## Write it to stdout.
    _calculate_fraction_duplicates_for_all_metric_files(
        {
            filehandle          => $filehandle,
            outfile_path_prefix => $outfile_path_prefix,
        }
    );

    ## Close filehandles
    close $filehandle;

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => $outfile_name_prefix,
                path             => catfile( $outfile_path_prefix . $UNDERSCORE . q{metric} ),
                recipe_name      => q{markduplicates},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_sample},
                job_id_chain                      => $job_id_chain,
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

sub analysis_markduplicates_rna {

## Function : Mark duplicated rna reads using Picardtools markduplicates on files generated from alignment (sorted, merged).
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
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
    my $profile_base_command;
    my $recipe_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
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

    use MIP::Get::File qw{ get_merged_infile_prefix get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Program::Gnu::Coreutils qw{ gnu_cat };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Picardtools qw{ picardtools_markduplicates picardtools_gatherbamfiles };
    use MIP::Program::Samtools qw{ samtools_flagstat samtools_index samtools_view };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Sample_info qw{
      set_file_path_to_store
      set_recipe_outfile_in_sample_info };
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
    my %infile_path = %{ $io{in}{file_path_href} };

    my %rec_atr = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $job_id_chain       = $rec_atr{chain};
    my $recipe_mode        = $active_parameter_href->{$recipe_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $xargs_file_path_prefix;
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number = $recipe_resource{core_number};

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Outpaths
    ## Assign suffix
    my $outfile_suffix = $rec_atr{outfile_suffix};
    my $outsample_directory =
      catdir( $active_parameter_href->{outdata_dir}, $sample_id, $recipe_name );
    my $outfile_tag =
      $file_info_href->{$sample_id}{$recipe_name}{file_tag};
    my @outfile_paths =
      map {
        catdir( $outsample_directory,
            $merged_infile_prefix . $outfile_tag . $DOT . $_ . $outfile_suffix )
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
    my $outfile_name_prefix = $io{out}{file_name_prefix};
    my %outfile_path        = %{ $io{out}{file_path_href} };
    my $outfile_path_prefix = $io{out}{file_path_prefix};

    ## Update memory and parallel processes for markduplicates
    my ( $recipe_memory, $contig_memory_href, $parallel_processes ) = _get_markdup_resources(
        {
            active_contigs_ref  => $file_info_href->{bam_contigs},
            node_memory         => $active_parameter_href->{node_ram_memory},
            primary_contigs_ref => $file_info_href->{primary_contigs},
            recipe_core_number  => $core_number,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle      = IO::Handle->new();
    my $xargsfilehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $parallel_processes,
            directory_id                    => $sample_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            log                             => $log,
            memory_allocation               => $recipe_memory,
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    ## Marking Duplicates
    say {$filehandle} q{## Marking Duplicates};

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

        my $stderrfile_path = $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        my $metrics_file    = $outfile_path_prefix . $DOT . $contig . $DOT . q{metric};
        picardtools_markduplicates(
            {
                create_index     => q{true},
                filehandle       => $xargsfilehandle,
                infile_paths_ref => [ $infile_path{$contig} ],
                java_jar => catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
                java_use_large_pages       => $active_parameter_href->{java_use_large_pages},
                memory_allocation          => q{Xmx} . $contig_memory_href->{$contig} . q{g},
                metrics_file               => $metrics_file,
                optical_duplicate_distance =>
                  $active_parameter_href->{markduplicates_picardtools_opt_dup_dist},
                outfile_path       => $outfile_path{$contig},
                referencefile_path => $referencefile_path,
                stderrfile_path    => $stderrfile_path,
                temp_directory     => $temp_directory,
            }
        );
        print {$xargsfilehandle} $SEMICOLON . $SPACE;

        ## Process BAM with samtools flagstat to produce metric file for downstream analysis
        samtools_flagstat(
            {
                filehandle      => $xargsfilehandle,
                infile_path     => $outfile_path{$contig},
                stdoutfile_path => $outfile_path_prefix . $DOT . $contig . $UNDERSCORE . q{metric},
                stderrfile_path_append => $stderrfile_path,
            }
        );
        say {$xargsfilehandle} $NEWLINE;
    }

    ## Concatenate all metric files
    gnu_cat(
        {
            filehandle       => $filehandle,
            infile_paths_ref =>
              [ $outfile_path_prefix . $DOT . $ASTERISK . $UNDERSCORE . q{metric} ],
            stdoutfile_path => $outfile_path_prefix . $UNDERSCORE . q{metric_all},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Collect duplicate reads and reads mapped across all metric contig files. Calculate fraction duplicates.
    ## Write it to stdout.
    _calculate_fraction_duplicates_for_all_metric_files(
        {
            filehandle          => $filehandle,
            outfile_path_prefix => $outfile_path_prefix,
        }
    );

    ## Gather bam files in contig order
    my @gather_infile_paths =
      map { $outfile_path{$_} } @{ $file_info_href->{bam_contigs} };
    my $store_outfile_path   = $outfile_path_prefix . $outfile_suffix;
    my $store_outfile_suffix = $outfile_suffix;

    picardtools_gatherbamfiles(
        {
            create_index     => q{true},
            filehandle       => $filehandle,
            infile_paths_ref => \@gather_infile_paths,
            java_jar => catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
            java_use_large_pages => $active_parameter_href->{java_use_large_pages},
            memory_allocation    => q{Xmx4g},
            outfile_path         => $outfile_path_prefix . $outfile_suffix,
            referencefile_path   => $referencefile_path,
            temp_directory       => $temp_directory,
        }
    );
    say {$filehandle} $NEWLINE;

    $store_outfile_path   = $outfile_path_prefix . $DOT . q{cram};
    $store_outfile_suffix = $DOT . q{cram};

    ## Create BAM to CRAM for long term storage
    say {$filehandle} q{## BAM to CRAM for long term storage};
    samtools_view(
        {
            filehandle         => $filehandle,
            infile_path        => $outfile_path_prefix . $outfile_suffix,
            outfile_path       => $store_outfile_path,
            output_format      => q{cram},
            referencefile_path => $referencefile_path,
            thread_number      => $core_number,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Index CRAM
    samtools_index(
        {
            filehandle  => $filehandle,
            infile_path => $store_outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Close filehandles
    close $xargsfilehandle;
    close $filehandle;

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => $outfile_name_prefix,
                path             => catfile( $outfile_path_prefix . $UNDERSCORE . q{metric} ),
                recipe_name      => q{markduplicates},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        set_recipe_outfile_in_sample_info(
            {
                infile           => $outfile_name_prefix . $store_outfile_suffix,
                path             => $store_outfile_path,
                recipe_name      => $recipe_name . $UNDERSCORE . q{alignment},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        set_file_path_to_store(
            {
                format           => q{cram},
                id               => $sample_id,
                path             => $store_outfile_path,
                path_index       => $store_outfile_path . $DOT . q{crai},
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_sample},
                job_id_chain                      => $job_id_chain,
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

sub _calculate_fraction_duplicates_for_all_metric_files {

## Function : Collect duplicate reads and reads mapped across all metric contig files. Calculate fraction duplicates. Write it to stdout.
## Returns  :
## Arguments: $filehandle          => Filehandle to write to
##          : $outfile_path_prefix => Outfile path prefix

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $outfile_path_prefix;

    my $tmpl = {
        filehandle          => { defined => 1, required => 1, store => \$filehandle, },
        outfile_path_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path_prefix,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Sums all mapped and duplicate reads and takes fraction of before finishing
    # Execute perl
    my $regexp = q?perl -nae'?;

    # Iniate has to store duplicates and read mapped
    $regexp .= q?my %feature; ?;

    # Read file line by line
    $regexp .= q?while (<>) { ?;

    # Find duplicate metric
    $regexp .= q?if($_=~/duplicates/ && $_=~/^(\d+)/) { ?;

    # Add to previous duplicate metrics
    $regexp .= q?$feature{dup} = $feature{dup} + $1 } ?;

    # Find reads mapped
    $regexp .= q?if($_=~/\d+\smapped/ && $_=~/^(\d+)/) { ?;

    # Add to previous reads mapped
    $regexp .= q?$feature{map} = $feature{map} + $1} ?;

    # End of while loop
    $regexp .= q?} ?;

    # Print metrics to stdout
    $regexp .= q?print "Read Mapped: ".$feature{map}."\nDuplicates: ".$feature{dup}."\n".?;

    # Print Fraction duplicates to stdout
    $regexp .= q?"Fraction Duplicates: ".$feature{dup}/$feature{map}, "\n"; ?;

    # Quit
    $regexp .= q?last;'?;

    ## Sum metric over concatenated file
    print {$filehandle} $regexp . $SPACE;
    print {$filehandle} $outfile_path_prefix . $UNDERSCORE . q{metric_all} . $SPACE;
    say   {$filehandle} q{>} . $SPACE . $outfile_path_prefix . $UNDERSCORE . q{metric} . $SPACE,
      $NEWLINE;
    return;
}

sub _get_markdup_resources {

## Function : Calculate and return resources for markduplicates recipe
## Returns  : $recipe_memory, $contig_memory_href, $parallel_processes
## Arguments: $active_contigs_ref  => Active contigs {REF}
##          : $node_memory         => Available memory
##          : $primary_contigs_ref => Primary contigs {REF}
##          : $recipe_core_number  => Allocated core number

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_contigs_ref;
    my $node_memory;
    my $primary_contigs_ref;
    my $recipe_core_number;

    my $tmpl = {
        active_contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$active_contigs_ref,
            strict_type => 1,
        },
        node_memory => {
            defined     => 1,
            required    => 1,
            store       => \$node_memory,
            strict_type => 1,
        },
        primary_contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$primary_contigs_ref,
            strict_type => 1,
        },
        recipe_core_number => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_core_number,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use List::Util qw{ sum };
    use MIP::Environment::Cluster qw{ check_recipe_memory_allocation };

    Readonly my $CONTIG_MEM   => 6;
    Readonly my $CONTIG_2_MEM => 12;

    ## Java memory per contig
    my %contig_mem_alloc = map { $_ => $CONTIG_MEM } @{$primary_contigs_ref};
    $contig_mem_alloc{ $primary_contigs_ref->[1] } = $CONTIG_2_MEM;

    ## Slice hash on active contigs
    %contig_mem_alloc = %contig_mem_alloc{ @{$active_contigs_ref} };

    ## Calculate max recipe memory in use
    my @process_mem_allocs        = map { $_ + $JAVA_GUEST_OS_MEMORY } values %contig_mem_alloc;
    my @sorted_process_mem_allocs = reverse sort { $a <=> $b } @process_mem_allocs;
    my $max_process_mem           = sum @sorted_process_mem_allocs[ 0 .. $recipe_core_number - 1 ];
    my $recipe_memory             = check_recipe_memory_allocation(
        {
            node_ram_memory          => $node_memory,
            recipe_memory_allocation => $max_process_mem,
        }
    );

    ## Get number of parallel processes given recipe memory
    my $parallel_processes;
    my $memory_requirement = 0;

  MEMORY:
    foreach my $memory_alloc (@sorted_process_mem_allocs) {

        last if ( $memory_requirement >= $recipe_memory );

        $parallel_processes++;
        $memory_requirement += $memory_alloc;
    }

    return $recipe_memory, \%contig_mem_alloc, $parallel_processes;
}
1;
