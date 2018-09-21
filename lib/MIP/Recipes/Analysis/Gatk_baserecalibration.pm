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
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ analysis_gatk_baserecalibration analysis_gatk_baserecalibration_rio };

}

## Constants
Readonly my $ASTERIX    => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};
Readonly my $MINUS_ONE  => -1;

sub analysis_gatk_baserecalibration {

## Function : GATK baserecalibrator/ApplyBQSR to recalibrate bases before variant calling. Both BaseRecalibrator/ApplyBQSR will be executed within the same sbatch script.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File info hash {REF}
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
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
    my $job_id_href;
    my $parameter_href;
    my $program_info_path;
    my $program_name;
    my $sample_info_href;
    my $sample_id;

    ## Default(s)
    my $family_id;
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
    use MIP::File::Interval qw{ generate_contig_interval_file };
    use MIP::Get::File
      qw{ get_exom_target_bed_file get_file_suffix get_merged_infile_prefix get_io_files };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::Program::Alignment::Gatk
      qw{ gatk_baserecalibrator gatk_applybqsr };
    use MIP::Program::Alignment::Picardtools qw{ picardtools_gatherbamfiles };
    use MIP::QC::Record
      qw{ add_program_outfile_to_sample_info add_program_metafile_to_sample_info add_processing_metafile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
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
    my $indir_path_prefix  = $io{in}{dir_path_prefix};
    my $infile_suffix      = $io{in}{file_suffix};
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my @temp_infile_paths  = @{ $io{temp}{file_paths} };

    my $job_id_chain       = $parameter_href->{$program_name}{chain};
    my $program_mode       = $active_parameter_href->{$program_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $analysis_type = $active_parameter_href->{analysis_type}{$sample_id};
    my $xargs_file_path_prefix;
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Outpaths
    ## Assign suffix
    my $outfile_suffix = get_file_suffix(
        {
            jobid_chain    => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
        }
    );
    my $outsample_directory =
      catdir( $active_parameter_href->{outdata_dir}, $sample_id,
        $program_name );
    my $outfile_tag =
      $file_info_href->{$sample_id}{$program_name}{file_tag};
    my @outfile_paths =
      map {
        catdir( $outsample_directory,
            $merged_infile_prefix . $outfile_tag . $DOT . $_ . $outfile_suffix )
      } @{ $file_info_href->{contigs_size_ordered} };

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

    my $outdir_path_prefix       = $io{out}{dir_path_prefix};
    my $outfile_name_prefix      = $io{out}{file_name_prefix};
    my $temp_outfile_path_prefix = $io{temp}{file_path_prefix};
    my @temp_outfile_paths       = @{ $io{temp}{file_paths} };

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
            log                             => $log,
            process_time                    => $time,
            program_directory               => $program_name,
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
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

    ### SHELL:

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
            indirectory        => $indir_path_prefix,
            infile             => $infile_name_prefix,
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
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    while ( my ( $infile_index, $contig ) =
        each @{ $file_info_href->{contigs_size_ordered} } )
    {

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

        my $base_quality_score_recalibration_file =
          $temp_outfile_paths[$infile_index] . $DOT . q{grp};
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        gatk_baserecalibrator(
            {
                FILEHANDLE    => $XARGSFILEHANDLE,
                infile_path       => $temp_infile_paths[$infile_index],
                intervals_ref => \@intervals,
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                memory_allocation => q{Xmx6g},
                known_sites_ref   => \@{
                    $active_parameter_href->{gatk_baserecalibration_known_sites}
                },
                outfile_path       => $base_quality_score_recalibration_file,
                referencefile_path => $referencefile_path,
                stderrfile_path    => $stderrfile_path,
                temp_directory     => $temp_directory,
                verbosity    => $active_parameter_href->{gatk_logging_level},
                xargs_mode         => 1,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## GATK ApplyBQSR
    say {$FILEHANDLE} q{## GATK ApplyBQSR};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    while ( my ( $infile_index, $contig ) =
        each @{ $file_info_href->{contigs_size_ordered} } )
    {

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

        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        my $base_quality_score_recalibration_file =
          $temp_file_path_prefix . $UNDERSCORE . $contig . $DOT . q{grp};
        gatk_applybqsr(
            {
                base_quality_score_recalibration_file =>
                  $base_quality_score_recalibration_file,
                FILEHANDLE    => $XARGSFILEHANDLE,
                infile_path   => $temp_infile_paths[$infile_index],
                intervals_ref => \@intervals,
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                memory_allocation => q{Xmx6g},
                verbosity => $active_parameter_href->{gatk_logging_level},
                read_filters_ref => \@{
                    $active_parameter_href
                      ->{gatk_baserecalibration_read_filters}
                },
                referencefile_path         => $referencefile_path,
                static_quantized_quals_ref => \@{
                    $active_parameter_href
                      ->{gatk_baserecalibration_static_quantized_quals}
                },
                outfile_path     => $temp_outfile_paths[$infile_index],
                referencefile_path => $referencefile_path,
                stderrfile_path    => $stderrfile_path,
                temp_directory     => $temp_directory,
                xargs_mode         => 1,
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
            file_ending        => substr( $outfile_suffix, 0, 2 ) . $ASTERIX,
            file_path          => $file_path,
            outdirectory       => $outdir_path_prefix,
            outfile            => $outfile_name_prefix,
            program_info_path  => $program_info_path,
            temp_directory     => $temp_directory,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Gather BAM files
    say {$FILEHANDLE} q{## Gather BAM files};

    ## Assemble infile paths in contig order and not per size
    my @gather_infile_paths =
      map { catdir( $temp_outfile_path_prefix . $DOT . $_ . $outfile_suffix ) }
      @{ $file_info_href->{contigs} };

    picardtools_gatherbamfiles(
        {
            create_index     => q{true},
            FILEHANDLE       => $FILEHANDLE,
            infile_paths_ref => \@gather_infile_paths,
            java_jar         => catfile(
                $active_parameter_href->{picardtools_path},
                q{picard.jar}
            ),
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            memory_allocation  => q{Xmx4g},
            outfile_path       => $temp_outfile_path_prefix . $outfile_suffix,
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
            infile_path => $temp_outfile_path_prefix
              . substr( $outfile_suffix, 0, 2 )
              . $ASTERIX,
            outfile_path => $outdir_path_prefix,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $XARGSFILEHANDLE;
    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        my $gathered_outfile_path =
          catfile( $outdir_path_prefix,
            $outfile_name_prefix . $outfile_suffix );

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                infile           => $outfile_name_prefix,
                path             => $gathered_outfile_path,
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
                path             => $gathered_outfile_path,
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

## Function : GATK baserecalibrator/ApplyBQSR to recalibrate bases before variant calling. Both BaseRecalibrator/ApplyBQSR will be executed within the same sbatch script.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $FILEHANDLE              => Filehandle to write to
##          : $file_info_href          => File info hash {REF}
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
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
    my $FILEHANDLE;
    my $file_info_href;
    my $file_path;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $program_info_path;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
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
      qw{ get_exom_target_bed_file get_file_suffix get_merged_infile_prefix get_io_files};
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::Program::Alignment::Gatk
      qw{ gatk_baserecalibrator gatk_applybqsr };
    use MIP::Program::Alignment::Picardtools qw{ picardtools_gatherbamfiles };
    use MIP::QC::Record
      qw{ add_program_outfile_to_sample_info add_program_metafile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };

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
    my $indir_path_prefix       = $io{in}{dir_path_prefix};
    my $infile_suffix           = $io{in}{file_suffix};
    my $infile_name_prefix      = $io{in}{file_name_prefix};
    my $temp_infile_name_prefix = $io{temp}{file_name_prefix};
    my @temp_infile_paths       = @{ $io{temp}{file_paths} };

    my $job_id_chain       = $parameter_href->{$program_name}{chain};
    my $program_mode       = $active_parameter_href->{$program_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $analysis_type = $active_parameter_href->{analysis_type}{$sample_id};
    my $xargs_file_path_prefix;
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Outpaths
    ## Assign suffix
    my $outfile_suffix = get_file_suffix(
        {
            jobid_chain    => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
        }
    );
    my $outsample_directory =
      catdir( $active_parameter_href->{outdata_dir}, $sample_id,
        $program_name );
    my $outfile_tag =
      $file_info_href->{$sample_id}{$program_name}{file_tag};
    my @outfile_paths =
      map {
        catdir( $outsample_directory,
            $merged_infile_prefix . $outfile_tag . $DOT . $_ . $outfile_suffix )
      } @{ $file_info_href->{contigs_size_ordered} };

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

    my $outdir_path_prefix       = $io{out}{dir_path_prefix};
    my $outfile_name_prefix      = $io{out}{file_name_prefix};
    my $temp_outfile_path_prefix = $io{temp}{file_path_prefix};
    my @temp_outfile_paths       = @{ $io{temp}{file_paths} };

    ## Filehandles
    # Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get exome_target_bed file for specfic sample_id and add file_ending from file_info hash if supplied
    my $exome_target_bed_file = get_exom_target_bed_file(
        {
            exome_target_bed_href => $active_parameter_href->{exome_target_bed},
            file_ending           => $file_info_href->{exome_target_bed}[0],
            log                   => $log,
            sample_id             => $sample_id,
        }
    );

    ### SHELL:

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
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    while ( my ( $infile_index, $contig ) =
        each @{ $file_info_href->{contigs_size_ordered} } )
    {

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

        my $base_quality_score_recalibration_file =
          $temp_outfile_paths[$infile_index] . $DOT . q{grp};
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        gatk_baserecalibrator(
            {
                FILEHANDLE    => $XARGSFILEHANDLE,
                infile_path       => $temp_infile_paths[$infile_index],
                intervals_ref => \@intervals,
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                memory_allocation => q{Xmx6g},
                known_sites_ref   => \@{
                    $active_parameter_href->{gatk_baserecalibration_known_sites}
                },
                verbosity    => $active_parameter_href->{gatk_logging_level},
                outfile_path       => $base_quality_score_recalibration_file,
                referencefile_path => $referencefile_path,
                stderrfile_path    => $stderrfile_path,
                temp_directory     => $temp_directory,
                xargs_mode         => 1,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## GATK PrintReads
    say {$FILEHANDLE} q{## GATK ApplyBQSR};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    while ( my ( $infile_index, $contig ) =
        each @{ $file_info_href->{contigs_size_ordered} } )
    {

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

        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        my $base_quality_score_recalibration_file =
          $temp_outfile_paths[$infile_index] . $DOT . q{grp};
        gatk_applybqsr(
            {
                base_quality_score_recalibration_file =>
                  $base_quality_score_recalibration_file,
                FILEHANDLE    => $XARGSFILEHANDLE,
                infile_path   => $temp_infile_paths[$infile_index],
                intervals_ref => \@intervals,
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                memory_allocation => q{Xmx6g},
                verbosity => $active_parameter_href->{gatk_logging_level},
                read_filters_ref => \@{
                    $active_parameter_href
                      ->{gatk_baserecalibration_read_filters}
                },
                static_quantized_quals_ref => \@{
                    $active_parameter_href
                      ->{gatk_baserecalibration_static_quantized_quals}
                },
                outfile_path     => $temp_outfile_paths[$infile_index],
                referencefile_path => $referencefile_path,
                stderrfile_path    => $stderrfile_path,
                temp_directory     => $temp_directory,
                xargs_mode         => 1,
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
            file_ending        => substr( $outfile_suffix, 0, 2 ) . $ASTERIX,
            file_path          => $file_path,
            core_number        => $core_number,
            outdirectory       => $outdir_path_prefix,
            outfile            => $outfile_name_prefix,
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
            file_name         => $temp_infile_name_prefix,
            indirectory       => $temp_directory,
        }
    );

    ## Gather BAM files
    say {$FILEHANDLE} q{## Gather BAM files};

    ## Assemble infile paths in contig order and not per size
    my @gather_infile_paths =
      map { catdir( $temp_outfile_path_prefix . $DOT . $_ . $outfile_suffix ) }
      @{ $file_info_href->{contigs} };

    picardtools_gatherbamfiles(
        {
            create_index     => q{true},
            FILEHANDLE       => $FILEHANDLE,
            infile_paths_ref => \@gather_infile_paths,
            java_jar         => catfile(
                $active_parameter_href->{picardtools_path},
                q{picard.jar}
            ),
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            memory_allocation  => q{Xmx4g},
            outfile_path       => $temp_outfile_path_prefix . $outfile_suffix,
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
            infile_path => $temp_outfile_path_prefix
              . substr( $outfile_suffix, 0, 2 )
              . $ASTERIX,
            outfile_path => $outdir_path_prefix,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $XARGSFILEHANDLE;
    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        my $gathered_outfile_path =
          catfile( $outdir_path_prefix,
            $outfile_name_prefix . $outfile_suffix );

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                infile           => $outfile_name_prefix,
                path             => $gathered_outfile_path,
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
                path             => $gathered_outfile_path,
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
