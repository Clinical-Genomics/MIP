package MIP::Recipes::Analysis::Gatk_realigner;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use autodie qw{ :all };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Basename qw{ basename };
use File::Spec::Functions qw{ catdir catfile devnull };
use POSIX;

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_gatk_realigner analysis_gatk_realigner_rio };

}

##Constants
Readonly my $ASTERIX    => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_gatk_realigner {

## Function : GATK ReAlignerTargetCreator/IndelRealigner to rearrange reads around INDELs. Both ReAlignerTargetCreator and IndelRealigner will be executed within the same sbatch script.
## Returns  :
## Arguments: $parameter_href          => Parameter hash {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $sample_id               => Sample id
##          : $insample_directory      => In sample directory
##          : $outsample_directory     => Out sample directory
##          : $program_name            => Program name
##          : $program_info_path       => The program info path
##          : $file_path               => File path
##          : $family_id               => Family id
##          : $temp_directory          => Temporary directory
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $sample_id;
    my $insample_directory;
    my $outsample_directory;
    my $program_name;
    my $program_info_path;
    my $file_path;

    ## Default(s)
    my $family_id;
    my $temp_directory;
    my $outaligner_dir;
    my $xargs_file_counter;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
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
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
        },
        insample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$insample_directory,
        },
        outsample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outsample_directory,
        },
        program_name      => { strict_type => 1, store => \$program_name },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        file_path         => { strict_type => 1, store => \$file_path },
        family_id         => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
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
      qw{ get_file_suffix get_merged_infile_prefix get_exom_target_bed_file };
    use MIP::IO::Files qw{ xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::Program::Alignment::Gatk
      qw{ gatk_realignertargetcreator gatk_indelrealigner };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $analysis_type = $active_parameter_href->{analysis_type}{$sample_id};
    my $xargs_file_path_prefix;

    ## Filehandles
    # Create anonymous filehandel
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $sample_id,
            program_name          => $program_name,
            program_directory     => catfile($outaligner_dir),
            core_number           => $core_number,
            process_time          => $time,
            temp_directory        => $temp_directory,
        }
    );

    ## Used downstream
    $parameter_href->{$mip_program_name}{$sample_id}{indirectory} =
      $outsample_directory;

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$sample_id}{pmarkduplicates}{file_tag};
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
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
            jobid_chain    => $job_id_chain,
        }
    );

    ## Get exome_target_bed file for specfic sample_id and add file_ending from file_info hash if supplied
    my $exome_target_bed_file = get_exom_target_bed_file(
        {
            exome_target_bed_href => $active_parameter_href->{exome_target_bed},
            sample_id             => $sample_id,
            log                   => $log,
            file_ending           => $file_info_href->{exome_target_bed}[2],
        }
    );

    ## Exome analysis
    if ( $analysis_type eq q{wes} ) {

        ## Generate contig specific interval_list
        generate_contig_interval_file(
            {
                contigs_ref => \@{ $file_info_href->{contigs_size_ordered} },
                exome_target_bed_file => $exome_target_bed_file,
                reference_dir      => $active_parameter_href->{reference_dir},
                outdirectory       => $temp_directory,
                max_cores_per_node => $core_number,
                FILEHANDLE         => $FILEHANDLE,
            }
        );

        ## Reroute to only filename
        $exome_target_bed_file = basename($exome_target_bed_file);
    }

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    ($xargs_file_counter) = xargs_migrate_contig_files(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            contigs_ref        => \@{ $file_info_href->{contigs_size_ordered} },
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
            infile             => $infile_prefix,
            indirectory        => $insample_directory,
            file_ending        => substr( $infile_suffix, 0, 2 ) . $ASTERIX,
            temp_directory     => $temp_directory,
        }
    );

    ## GATK ReAlignerTargetCreator
    say {$FILEHANDLE} q{## GATK ReAlignerTargetCreator};

    Readonly my $JAVA_MEMORY_ALLOCATION => 4;

    ## Division by 4 since the java heap is 4GB
    $core_number = floor(
        $active_parameter_href->{node_ram_memory} / $JAVA_MEMORY_ALLOCATION );

    ## Limit number of cores requested to the maximum number of cores available per node
    $core_number = check_max_core_number(
        {
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
            core_number_requested => $core_number
        }
    );

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
            first_command      => q{java},
            memory_allocation  => q{Xmx4g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $temp_directory,
            java_jar       => catfile(
                $active_parameter_href->{gatk_path},
                q{GenomeAnalysisTK.jar}
            ),
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

            # Per contig
            @intervals = ($contig);
        }

        my $infile_path =
          $file_path_prefix . $UNDERSCORE . $contig . $infile_suffix;
        my $realign_outfile_path =
          $outfile_path_prefix . $UNDERSCORE . $contig . $DOT . q{intervals};
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        gatk_realignertargetcreator(
            {
                known_alleles_ref => \@{
                    $active_parameter_href->{gatk_realigner_indel_known_sites}
                },
                intervals_ref      => \@intervals,
                infile_path        => $infile_path,
                outfile_path       => $realign_outfile_path,
                stderrfile_path    => $stderrfile_path,
                referencefile_path => $referencefile_path,
                logging_level => $active_parameter_href->{gatk_logging_level},
                downsample_to_coverage =>
                  $active_parameter_href->{gatk_downsample_to_coverage},
                gatk_disable_auto_index_and_file_lock => $active_parameter_href
                  ->{gatk_disable_auto_index_and_file_lock},
                FILEHANDLE => $XARGSFILEHANDLE,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## GATK IndelRealigner
    say {$FILEHANDLE} q{## GATK IndelRealigner};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
            first_command      => q{java},
            memory_allocation  => q{Xmx4g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $temp_directory,
            java_jar       => catfile(
                $active_parameter_href->{gatk_path},
                q{GenomeAnalysisTK.jar}
            ),
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

            # Per contig
            @intervals = ($contig);
        }

        my $infile_path =
          $file_path_prefix . $UNDERSCORE . $contig . $infile_suffix;
        my $outfile_path =
          $outfile_path_prefix . $UNDERSCORE . $contig . $outfile_suffix;
        my $target_intervals_file =
          $outfile_path_prefix . $UNDERSCORE . $contig . $DOT . q{intervals};
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        gatk_indelrealigner(
            {
                known_alleles_ref => \@{
                    $active_parameter_href->{gatk_realigner_indel_known_sites}
                },
                intervals_ref         => \@intervals,
                infile_path           => $infile_path,
                outfile_path          => $outfile_path,
                target_intervals_file => $target_intervals_file,
                stderrfile_path       => $stderrfile_path,
                referencefile_path    => $referencefile_path,
                logging_level => $active_parameter_href->{gatk_logging_level},
                downsample_to_coverage =>
                  $active_parameter_href->{gatk_downsample_to_coverage},
                gatk_disable_auto_index_and_file_lock => $active_parameter_href
                  ->{gatk_disable_auto_index_and_file_lock},
                FILEHANDLE => $XARGSFILEHANDLE,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## Copies file from temporary directory. Per contig
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    ($xargs_file_counter) = xargs_migrate_contig_files(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            contigs_ref        => \@{ $file_info_href->{contigs_size_ordered} },
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
            outfile            => $outfile_prefix,
            outdirectory       => $outsample_directory,
            temp_directory     => $temp_directory,
            file_ending        => substr( $infile_suffix, 0, 2 ) . $ASTERIX,
        }
    );

    close $XARGSFILEHANDLE;
    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {

        slurm_submit_job_sample_id_dependency_add_to_sample(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                family_id               => $family_id,
                sample_id               => $sample_id,
                path                    => $job_id_chain,
                log                     => $log,
                sbatch_file_name        => $file_path
            }
        );
    }
    return;
}

sub analysis_gatk_realigner_rio {

## Function : GATK ReAlignerTargetCreator/IndelRealigner to rearrange reads around INDELs. Both ReAlignerTargetCreator and IndelRealigner will be executed within the same sbatch script.
## Returns  : $xargs_file_counter
## Arguments: $parameter_href          => Parameter hash {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $sample_id               => Sample id
##          : $program_name            => Program name
##          : $program_info_path       => The program info path
##          : $file_path               => File path
##          : $FILEHANDLE              => Filehandle to write to
##          : $family_id               => Family id
##          : $temp_directory          => Temporary directory
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $file_info_href;
    my $sample_id;
    my $program_name;
    my $program_info_path;
    my $file_path;
    my $FILEHANDLE;

    ## Default(s)
    my $family_id;
    my $temp_directory;
    my $outaligner_dir;
    my $xargs_file_counter;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
        },
        program_name      => { strict_type => 1, store => \$program_name },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        file_path         => { strict_type => 1, store => \$file_path },
        FILEHANDLE => { store => \$FILEHANDLE, },
        family_id  => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
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
      qw{ get_file_suffix get_merged_infile_prefix get_exom_target_bed_file };
    use MIP::IO::Files qw{ xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::Program::Alignment::Gatk
      qw{ gatk_realignertargetcreator gatk_indelrealigner };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;

    ## Alias
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $analysis_type = $active_parameter_href->{analysis_type}{$sample_id};
    my $xargs_file_path_prefix;

    ## Filehandles
    # Create anonymous filehandel
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Assign directories
    my $outsample_directory = catdir( $active_parameter_href->{outdata_dir},
        $sample_id, $outaligner_dir );

    ## Used downstream
    $parameter_href->{$mip_program_name}{$sample_id}{indirectory} =
      $outsample_directory;

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$sample_id}{pmarkduplicates}{file_tag};
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
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
            jobid_chain    => $job_id_chain,
        }
    );

    ## Get exome_target_bed file for specfic sample_id and add file_ending from file_info hash if supplied
    my $exome_target_bed_file = get_exom_target_bed_file(
        {
            exome_target_bed_href => $active_parameter_href->{exome_target_bed},
            sample_id             => $sample_id,
            log                   => $log,
            file_ending           => $file_info_href->{exome_target_bed}[2],
        }
    );

    ## Exome analysis
    if ( $analysis_type eq q{wes} ) {

        ## Generate contig specific interval_list
        generate_contig_interval_file(
            {
                contigs_ref => \@{ $file_info_href->{contigs_size_ordered} },
                exome_target_bed_file => $exome_target_bed_file,
                reference_dir      => $active_parameter_href->{reference_dir},
                outdirectory       => $temp_directory,
                max_cores_per_node => $core_number,
                FILEHANDLE         => $FILEHANDLE,
            }
        );

        ## Reroute to only filename
        $exome_target_bed_file = basename($exome_target_bed_file);
    }

    ## GATK ReAlignerTargetCreator
    say {$FILEHANDLE} q{## GATK ReAlignerTargetCreator};

    Readonly my $JAVA_MEMORY_ALLOCATION => 4;

    ## Division by 4 since the java heap is 4GB
    $core_number = floor(
        $active_parameter_href->{node_ram_memory} / $JAVA_MEMORY_ALLOCATION );

    ## Limit number of cores requested to the maximum number of cores available per node
    $core_number = check_max_core_number(
        {
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
            core_number_requested => $core_number
        }
    );

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
            first_command      => q{java},
            memory_allocation  => q{Xmx4g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $temp_directory,
            java_jar       => catfile(
                $active_parameter_href->{gatk_path},
                q{GenomeAnalysisTK.jar}
            ),
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

            # Per contig
            @intervals = ($contig);
        }

        my $infile_path =
          $file_path_prefix . $UNDERSCORE . $contig . $infile_suffix;
        my $realign_outfile_path =
          $outfile_path_prefix . $UNDERSCORE . $contig . $DOT . q{intervals};
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        gatk_realignertargetcreator(
            {
                known_alleles_ref => \@{
                    $active_parameter_href->{gatk_realigner_indel_known_sites}
                },
                intervals_ref      => \@intervals,
                infile_path        => $infile_path,
                outfile_path       => $realign_outfile_path,
                stderrfile_path    => $stderrfile_path,
                referencefile_path => $referencefile_path,
                logging_level => $active_parameter_href->{gatk_logging_level},
                downsample_to_coverage =>
                  $active_parameter_href->{gatk_downsample_to_coverage},
                gatk_disable_auto_index_and_file_lock => $active_parameter_href
                  ->{gatk_disable_auto_index_and_file_lock},
                FILEHANDLE => $XARGSFILEHANDLE,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## GATK IndelRealigner
    say {$FILEHANDLE} q{## GATK IndelRealigner};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
            first_command      => q{java},
            memory_allocation  => q{Xmx4g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $temp_directory,
            java_jar       => catfile(
                $active_parameter_href->{gatk_path},
                q{GenomeAnalysisTK.jar}
            ),
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

            # Per contig
            @intervals = ($contig);
        }

        my $infile_path =
          $file_path_prefix . $UNDERSCORE . $contig . $infile_suffix;
        my $outfile_path =
          $outfile_path_prefix . $UNDERSCORE . $contig . $outfile_suffix;
        my $target_intervals_file =
          $outfile_path_prefix . $UNDERSCORE . $contig . $DOT . q{intervals};
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        gatk_indelrealigner(
            {
                known_alleles_ref => \@{
                    $active_parameter_href->{gatk_realigner_indel_known_sites}
                },
                intervals_ref         => \@intervals,
                infile_path           => $infile_path,
                outfile_path          => $outfile_path,
                target_intervals_file => $target_intervals_file,
                stderrfile_path       => $stderrfile_path,
                referencefile_path    => $referencefile_path,
                logging_level => $active_parameter_href->{gatk_logging_level},
                downsample_to_coverage =>
                  $active_parameter_href->{gatk_downsample_to_coverage},
                gatk_disable_auto_index_and_file_lock => $active_parameter_href
                  ->{gatk_disable_auto_index_and_file_lock},
                FILEHANDLE => $XARGSFILEHANDLE,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## Remove file at temporary directory
    delete_contig_files(
        {
            file_elements_ref => \@{ $file_info_href->{contigs_size_ordered} },
            FILEHANDLE        => $FILEHANDLE,
            core_number       => $core_number,
            file_name         => $infile_prefix,
            file_ending       => substr( $infile_suffix, 0, 2 ) . $ASTERIX,
            indirectory       => $temp_directory,
        }
    );

    close $XARGSFILEHANDLE;

    ## Track the number of created xargs scripts per module for Block algorithm
    return $xargs_file_counter;
}

1;
