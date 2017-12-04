package MIP::Recipes::Analysis::Gatk_baserecalibration;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename };
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

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ analysis_gatk_baserecalibration analysis_gatk_baserecalibration_rio };

}

## Constants
Readonly my $ASTERIX    => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_gatk_baserecalibration {

## Function : GATK baserecalibrator/printreads to recalibrate bases before variant calling. Both BaseRecalibrator/PrintReads will be executed within the same sbatch script.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File info hash {REF}
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $insample_directory      => In sample directory
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outsample_directory     => Out sample directory
##          : $parameter_href          => Parameter hash {REF}
##          : $program_info_path       => The program info path
##          : $program_name            => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $file_path;
    my $infile_lane_prefix_href;
    my $insample_directory;
    my $job_id_href;
    my $outsample_directory;
    my $parameter_href;
    my $program_info_path;
    my $program_name;
    my $sample_info_href;
    my $sample_id;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
    my $temp_directory;
    my $xargs_file_counter;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        file_path          => { strict_type => 1, store => \$file_path },
        insample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$insample_directory,
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
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
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
            store       => \$parameter_href,
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        program_name      => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory,
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Cluster qw{ check_max_core_number };
    use MIP::Delete::File qw{ delete_contig_files };
    use MIP::File::Interval qw{ generate_contig_interval_file };
    use MIP::Get::File
      qw{ get_file_suffix get_merged_infile_prefix get_exom_target_bed_file};
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Language::Java qw{ java_core };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::Program::Alignment::Gatk
      qw{ gatk_baserecalibrator gatk_printreads };
    use MIP::Program::Alignment::Picardtools qw{ picardtools_gatherbamfiles };
    use MIP::QC::Record
      qw{ add_program_outfile_to_sample_info add_program_metafile_to_sample_info add_processing_metafile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain       = $parameter_href->{$mip_program_name}{chain};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $analysis_type = $active_parameter_href->{analysis_type}{$sample_id};
    my $xargs_file_path_prefix;
    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();
    my $FILEHANDLE      = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $sample_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            process_time                    => $time,
            program_directory               => catfile($outaligner_dir),
            program_name                    => $program_name,
            source_environment_commands_ref => [$source_environment_cmd],
            temp_directory                  => $temp_directory,
        }
    );

    $parameter_href->{$mip_program_name}{$sample_id}{indirectory} =
      $outsample_directory;    #Used downstream

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$sample_id}{pgatk_realigner}{file_tag};
    my $outfile_tag =
      $file_info_href->{$sample_id}{$mip_program_name}{file_tag};

    ## Files
    my $infile_prefix  = $merged_infile_prefix . $infile_tag;
    my $outfile_prefix = $merged_infile_prefix . $outfile_tag;

    ## Paths
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ## Assign suffix
    my $infile_suffix = my $outfile_suffix = get_file_suffix(
        {
            jobid_chain    => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
        }
    );

    ## Get exome_target_bed file for specfic sample_id and add file_ending from file_info hash if supplied
    my $exome_target_bed_file = get_exom_target_bed_file(
        {
            exome_target_bed_href => $active_parameter_href->{exome_target_bed},
            file_ending           => $file_info_href->{exome_target_bed}[0],
            log                   => $log,
            sample_id             => $sample_id,
        }
    );

    ## Exome analysis
    if ( $analysis_type eq q{wes} ) {

        ## Generate contig specific interval_list
        generate_contig_interval_file(
            {
                contigs_ref => \@{ $file_info_href->{contigs_size_ordered} },
                exome_target_bed_file => $exome_target_bed_file,
                FILEHANDLE            => $FILEHANDLE,
                file_ending           => $DOT . q{intervals},
                max_cores_per_node    => $core_number,
                outdirectory          => $temp_directory,
                reference_dir => $active_parameter_href->{reference_dir},
            }
        );

        ## Add required GATK ending and reroute to only filename
        $exome_target_bed_file =
          basename($exome_target_bed_file) . $DOT . q{intervals};
    }

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    ($xargs_file_counter) = xargs_migrate_contig_files(
        {
            contigs_ref        => \@{ $file_info_href->{contigs_size_ordered} },
            core_number        => $core_number,
            indirectory        => $insample_directory,
            infile             => $infile_prefix,
            FILEHANDLE         => $FILEHANDLE,
            file_ending        => substr( $infile_suffix, 0, 2 ) . $ASTERIX,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
            temp_directory     => $temp_directory,
        }
    );

    ## Division by X according to the java heap
    Readonly my $JAVA_MEMORY_ALLOCATION => 6;
    $core_number = floor(
        $active_parameter_href->{node_ram_memory} / $JAVA_MEMORY_ALLOCATION );

    ## Limit number of cores requested to the maximum number of cores available per node
    $core_number = check_max_core_number(
        {
            core_number_requested => $core_number,
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    ## GATK BaseRecalibrator
    say {$FILEHANDLE} q{## GATK BaseRecalibrator};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number   => $core_number,
            FILEHANDLE    => $FILEHANDLE,
            file_path     => $file_path,
            first_command => q{java},
            java_jar      => catfile(
                $active_parameter_href->{gatk_path},
                q{GenomeAnalysisTK.jar},
            ),
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            memory_allocation  => q{Xmx6g},
            program_info_path  => $program_info_path,
            temp_directory     => $temp_directory,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## Get parameters
        # Exome analysis
        my @intervals;
        if ( $analysis_type eq q{wes} ) {

            ## Limit to targets kit target file
            @intervals = (
                catfile(
                    $temp_directory,
                    $contig . $UNDERSCORE . $exome_target_bed_file
                )
            );
        }
        else {
            ## wgs

            ## Per contig
            @intervals = ($contig);
        }

        my $infile_path =
          $file_path_prefix . $UNDERSCORE . $contig . $infile_suffix;
        my $outfile_path =
          $file_path_prefix . $UNDERSCORE . $contig . $DOT . q{grp};
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        gatk_baserecalibrator(
            {
                covariates_ref => \@{
                    $active_parameter_href->{gatk_baserecalibration_covariates}
                },
                downsample_to_coverage =>
                  $active_parameter_href->{gatk_downsample_to_coverage},
                FILEHANDLE                            => $XARGSFILEHANDLE,
                gatk_disable_auto_index_and_file_lock => $active_parameter_href
                  ->{gatk_disable_auto_index_and_file_lock},
                infile_path       => $infile_path,
                intervals_ref     => \@intervals,
                known_alleles_ref => \@{
                    $active_parameter_href->{gatk_baserecalibration_known_sites}
                },
                logging_level => $active_parameter_href->{gatk_logging_level},
                num_cpu_threads_per_data_thread =>
                  $active_parameter_href->{max_cores_per_node},
                outfile_path       => $outfile_path,
                referencefile_path => $referencefile_path,
                stderrfile_path    => $stderrfile_path,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## GATK PrintReads
    say {$FILEHANDLE} q{## GATK PrintReads};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number   => $core_number,
            FILEHANDLE    => $FILEHANDLE,
            file_path     => $file_path,
            first_command => q{java},
            java_jar      => catfile(
                $active_parameter_href->{gatk_path},
                q{GenomeAnalysisTK.jar},
            ),
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            memory_allocation  => q{Xmx6g},
            program_info_path  => $program_info_path,
            temp_directory     => $temp_directory,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## Get parameters
        # Exome  analysis
        my @intervals;
        if ( $analysis_type eq q{wes} ) {

            ## Limit to targets kit target file
            @intervals = (
                catfile(
                    $temp_directory,
                    $contig . $UNDERSCORE . $exome_target_bed_file
                )
            );
        }
        else {
            ## wgs

            ## Per contig
            @intervals = ($contig);
        }

        my $infile_path =
          $file_path_prefix . $UNDERSCORE . $contig . $infile_suffix;
        my $outfile_path =
          $outfile_path_prefix . $UNDERSCORE . $contig . $outfile_suffix;
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        my $base_quality_score_recalibration_file =
          $file_path_prefix . $UNDERSCORE . $contig . $DOT . q{grp};
        gatk_printreads(
            {
                base_quality_score_recalibration_file =>
                  $base_quality_score_recalibration_file,
                disable_indel_qual => $active_parameter_href
                  ->{gatk_baserecalibration_disable_indel_qual},
                downsample_to_coverage =>
                  $active_parameter_href->{gatk_downsample_to_coverage},
                FILEHANDLE                            => $XARGSFILEHANDLE,
                gatk_disable_auto_index_and_file_lock => $active_parameter_href
                  ->{gatk_disable_auto_index_and_file_lock},
                infile_path   => $infile_path,
                intervals_ref => \@intervals,
                logging_level => $active_parameter_href->{gatk_logging_level},
                num_cpu_threads_per_data_thread =>
                  $active_parameter_href->{max_cores_per_node},
                read_filters_ref => \@{
                    $active_parameter_href
                      ->{gatk_baserecalibration_read_filters}
                },
                static_quantized_quals_ref => \@{
                    $active_parameter_href
                      ->{gatk_baserecalibration_static_quantized_quals}
                },
                outfile_path       => $outfile_path,
                referencefile_path => $referencefile_path,
                stderrfile_path    => $stderrfile_path,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## Copies file from temporary directory. Per contig for variant callers.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    ($xargs_file_counter) = xargs_migrate_contig_files(
        {
            contigs_ref        => \@{ $file_info_href->{contigs_size_ordered} },
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_ending        => substr( $infile_suffix, 0, 2 ) . $ASTERIX,
            file_path          => $file_path,
            outfile            => $outfile_prefix,
            outdirectory       => $outsample_directory,
            program_info_path  => $program_info_path,
            temp_directory     => $temp_directory,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Gather BAM files
    say {$FILEHANDLE} q{## Gather BAM files};

    ## Assemble infile paths by adding directory and file ending
    my @infile_paths = map {
        catfile( $temp_directory,
            $outfile_prefix . $UNDERSCORE . $_ . $outfile_suffix )
    } @{ $file_info_href->{contigs} };

    picardtools_gatherbamfiles(
        {
            create_index     => q{true},
            FILEHANDLE       => $FILEHANDLE,
            infile_paths_ref => \@infile_paths,
            java_jar         => catfile(
                $active_parameter_href->{picardtools_path},
                q{picard.jar}
            ),
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            memory_allocation => q{Xmx4g},
            outfile_path =>
              catfile( $temp_directory, $outfile_prefix . $outfile_suffix ),
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $outfile_path_prefix
              . substr( $infile_suffix, 0, 2 )
              . $ASTERIX,
            outfile_path => $outsample_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $XARGSFILEHANDLE;
    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {

        my $program_outfile_path =
          catfile( $outsample_directory, $outfile_prefix . $outfile_suffix );

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                infile           => $merged_infile_prefix,
                path             => $program_outfile_path,
                program_name     => q{gatk_baserecalibration},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );
        my $most_complete_format_key =
          q{most_complete} . $UNDERSCORE . substr $outfile_suffix, 1;
        add_processing_metafile_to_sample_info(
            {
                metafile_tag     => $most_complete_format_key,
                path             => $program_outfile_path,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        slurm_submit_job_sample_id_dependency_add_to_sample(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_id               => $sample_id,
                sbatch_file_name        => $file_path
            }
        );
    }
    return;
}

sub analysis_gatk_baserecalibration_rio {

## Function : GATK baserecalibrator/printreads to recalibrate bases before variant calling. Both BaseRecalibrator/PrintReads will be executed within the same sbatch script.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $FILEHANDLE              => Filehandle to write to
##          : $file_info_href          => File info hash {REF}
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $outsample_directory     => Out sample directory
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $program_info_path       => The program info path
##          : $program_name            => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $FILEHANDLE;
    my $file_info_href;
    my $file_path;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $outsample_directory;
    my $parameter_href;
    my $program_name;
    my $program_info_path;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
    my $temp_directory;
    my $xargs_file_counter;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        FILEHANDLE     => { store => \$FILEHANDLE, },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        file_path               => { strict_type => 1, store => \$file_path },
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
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
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
            store       => \$parameter_href,
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        program_name      => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory,
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Cluster qw{ check_max_core_number };
    use MIP::Delete::File qw{ delete_contig_files };
    use MIP::File::Interval qw{ generate_contig_interval_file };
    use MIP::Get::File
      qw{ get_file_suffix get_merged_infile_prefix get_exom_target_bed_file};
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Language::Java qw{ java_core };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::Program::Alignment::Gatk
      qw{ gatk_baserecalibrator gatk_printreads };
    use MIP::Program::Alignment::Picardtools qw{ picardtools_gatherbamfiles };
    use MIP::QC::Record
      qw{ add_program_outfile_to_sample_info add_program_metafile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain       = $parameter_href->{$mip_program_name}{chain};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $analysis_type = $active_parameter_href->{analysis_type}{$sample_id};
    my $xargs_file_path_prefix;
    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();

    $parameter_href->{$mip_program_name}{$sample_id}{indirectory} =
      $outsample_directory;    #Used downstream

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$sample_id}{pgatk_realigner}{file_tag};
    my $outfile_tag =
      $file_info_href->{$sample_id}{$mip_program_name}{file_tag};

    ## Files
    my $infile_prefix  = $merged_infile_prefix . $infile_tag;
    my $outfile_prefix = $merged_infile_prefix . $outfile_tag;

    ## Paths
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ## Assign suffix
    my $infile_suffix = my $outfile_suffix = get_file_suffix(
        {
            jobid_chain    => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
        }
    );

    ## Get exome_target_bed file for specfic sample_id and add file_ending from file_info hash if supplied
    my $exome_target_bed_file = get_exom_target_bed_file(
        {
            exome_target_bed_href => $active_parameter_href->{exome_target_bed},
            file_ending           => $file_info_href->{exome_target_bed}[0],
            log                   => $log,
            sample_id             => $sample_id,
        }
    );

    ## Exome analysis
    if ( $analysis_type eq q{wes} ) {

        ## Generate contig specific interval_list
        generate_contig_interval_file(
            {
                contigs_ref => \@{ $file_info_href->{contigs_size_ordered} },
                exome_target_bed_file => $exome_target_bed_file,
                FILEHANDLE            => $FILEHANDLE,
                file_ending           => $DOT . q{intervals},
                max_cores_per_node    => $core_number,
                outdirectory          => $temp_directory,
                reference_dir => $active_parameter_href->{reference_dir},
            }
        );

        ## Add required GATK ending and reroute to only filename
        $exome_target_bed_file =
          basename($exome_target_bed_file) . $DOT . q{intervals};
    }

    ## Division by X according to the java heap
    Readonly my $JAVA_MEMORY_ALLOCATION => 6;
    $core_number = floor(
        $active_parameter_href->{node_ram_memory} / $JAVA_MEMORY_ALLOCATION );

    ## Limit number of cores requested to the maximum number of cores available per node
    $core_number = check_max_core_number(
        {
            core_number_requested => $core_number,
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    ## GATK BaseRecalibrator
    say {$FILEHANDLE} q{## GATK BaseRecalibrator};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number   => $core_number,
            FILEHANDLE    => $FILEHANDLE,
            file_path     => $file_path,
            first_command => q{java},
            java_jar      => catfile(
                $active_parameter_href->{gatk_path},
                q{GenomeAnalysisTK.jar},
            ),
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            memory_allocation  => q{Xmx6g},
            program_info_path  => $program_info_path,
            temp_directory     => $temp_directory,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## Get parameters
        # Exome analysis
        my @intervals;
        if ( $analysis_type eq q{wes} ) {

            ## Limit to targets kit target file
            @intervals = (
                catfile(
                    $temp_directory,
                    $contig . $UNDERSCORE . $exome_target_bed_file
                )
            );
        }
        else {
            ## wgs

            ## Per contig
            @intervals = ($contig);
        }

        my $infile_path =
          $file_path_prefix . $UNDERSCORE . $contig . $infile_suffix;
        my $outfile_path =
          $file_path_prefix . $UNDERSCORE . $contig . $DOT . q{grp};
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        gatk_baserecalibrator(
            {
                covariates_ref => \@{
                    $active_parameter_href->{gatk_baserecalibration_covariates}
                },
                downsample_to_coverage =>
                  $active_parameter_href->{gatk_downsample_to_coverage},
                FILEHANDLE                            => $XARGSFILEHANDLE,
                gatk_disable_auto_index_and_file_lock => $active_parameter_href
                  ->{gatk_disable_auto_index_and_file_lock},
                infile_path       => $infile_path,
                intervals_ref     => \@intervals,
                known_alleles_ref => \@{
                    $active_parameter_href->{gatk_baserecalibration_known_sites}
                },
                logging_level => $active_parameter_href->{gatk_logging_level},
                num_cpu_threads_per_data_thread =>
                  $active_parameter_href->{max_cores_per_node},
                outfile_path       => $outfile_path,
                referencefile_path => $referencefile_path,
                stderrfile_path    => $stderrfile_path,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## GATK PrintReads
    say {$FILEHANDLE} q{## GATK PrintReads};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number   => $core_number,
            FILEHANDLE    => $FILEHANDLE,
            file_path     => $file_path,
            first_command => q{java},
            java_jar      => catfile(
                $active_parameter_href->{gatk_path},
                q{GenomeAnalysisTK.jar},
            ),
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            memory_allocation  => q{Xmx6g},
            program_info_path  => $program_info_path,
            temp_directory     => $temp_directory,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## Get parameters
        # Exome  analysis
        my @intervals;
        if ( $analysis_type eq q{wes} ) {

            ## Limit to targets kit target file
            @intervals = (
                catfile(
                    $temp_directory,
                    $contig . $UNDERSCORE . $exome_target_bed_file
                )
            );
        }
        else {
            ## wgs

            ## Per contig
            @intervals = ($contig);
        }

        my $infile_path =
          $file_path_prefix . $UNDERSCORE . $contig . $infile_suffix;
        my $outfile_path =
          $outfile_path_prefix . $UNDERSCORE . $contig . $outfile_suffix;
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        my $base_quality_score_recalibration_file =
          $file_path_prefix . $UNDERSCORE . $contig . $DOT . q{grp};
        gatk_printreads(
            {
                base_quality_score_recalibration_file =>
                  $base_quality_score_recalibration_file,
                disable_indel_qual => $active_parameter_href
                  ->{gatk_baserecalibration_disable_indel_qual},
                downsample_to_coverage =>
                  $active_parameter_href->{gatk_downsample_to_coverage},
                FILEHANDLE                            => $XARGSFILEHANDLE,
                gatk_disable_auto_index_and_file_lock => $active_parameter_href
                  ->{gatk_disable_auto_index_and_file_lock},
                infile_path   => $infile_path,
                intervals_ref => \@intervals,
                logging_level => $active_parameter_href->{gatk_logging_level},
                num_cpu_threads_per_data_thread =>
                  $active_parameter_href->{max_cores_per_node},
                outfile_path     => $outfile_path,
                read_filters_ref => \@{
                    $active_parameter_href
                      ->{gatk_baserecalibration_read_filters}
                },
                referencefile_path         => $referencefile_path,
                static_quantized_quals_ref => \@{
                    $active_parameter_href
                      ->{gatk_baserecalibration_static_quantized_quals}
                },
                stderrfile_path => $stderrfile_path,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## Copies file from temporary directory. Per contig for variant callers.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    ($xargs_file_counter) = xargs_migrate_contig_files(
        {
            contigs_ref        => \@{ $file_info_href->{contigs_size_ordered} },
            FILEHANDLE         => $FILEHANDLE,
            file_ending        => substr( $infile_suffix, 0, 2 ) . $ASTERIX,
            file_path          => $file_path,
            core_number        => $core_number,
            outdirectory       => $outsample_directory,
            outfile            => $outfile_prefix,
            program_info_path  => $program_info_path,
            temp_directory     => $temp_directory,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Remove file at temporary directory
    delete_contig_files(
        {
            core_number       => $core_number,
            FILEHANDLE        => $FILEHANDLE,
            file_elements_ref => \@{ $file_info_href->{contigs_size_ordered} },
            file_ending       => substr( $infile_suffix, 0, 2 ) . $ASTERIX,
            file_name         => $infile_prefix,
            indirectory       => $temp_directory,
        }
    );

    ## Gather BAM files
    say {$FILEHANDLE} q{## Gather BAM files};

    ## Assemble infile paths by adding directory and file ending
    my @infile_paths = map {
        catfile( $temp_directory,
            $outfile_prefix . $UNDERSCORE . $_ . $outfile_suffix )
    } @{ $file_info_href->{contigs} };

    picardtools_gatherbamfiles(
        {
            create_index     => q{true},
            FILEHANDLE       => $FILEHANDLE,
            infile_paths_ref => \@infile_paths,
            java_jar         => catfile(
                $active_parameter_href->{picardtools_path},
                q{picard.jar}
            ),
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            memory_allocation => q{Xmx4g},
            outfile_path =>
              catfile( $temp_directory, $outfile_prefix . $outfile_suffix ),
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $outfile_path_prefix
              . substr( $infile_suffix, 0, 2 )
              . $ASTERIX,
            outfile_path => $outsample_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $XARGSFILEHANDLE;
    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {

        my $program_outfile_path =
          catfile( $outsample_directory, $outfile_prefix . $outfile_suffix );

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                infile           => $merged_infile_prefix,
                path             => $program_outfile_path,
                program_name     => q{gatk_baserecalibration},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );
        my $most_complete_format_key =
          q{most_complete} . $UNDERSCORE . substr $outfile_suffix, 1;
        add_processing_metafile_to_sample_info(
            {
                metafile_tag     => $most_complete_format_key,
                path             => $program_outfile_path,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        slurm_submit_job_sample_id_dependency_add_to_sample(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_id               => $sample_id,
                sbatch_file_name        => $file_path
            }
        );
    }
    return;
}

1;
