package MIP::Recipes::Analysis::Picardtools_mergesamfiles;

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

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ analysis_picardtools_mergesamfiles analysis_picardtools_mergesamfiles_rio };

}

## Constants
Readonly my $ASTERIX      => q{*};
Readonly my $DOT          => q{.};
Readonly my $EMPTY_STRING => q{};
Readonly my $NEWLINE      => qq{\n};
Readonly my $SEMICOLON    => q{;};
Readonly my $SPACE        => q{ };
Readonly my $UNDERSCORE   => q{_};

sub analysis_picardtools_mergesamfiles {

## Function : Merges all bam files using Picardtools mergesamfiles within each sampleid and files generated previously (option if provided with '-picardtools_mergesamfiles_previous_bams'). The merged files have to be sorted before attempting to merge.
## Returns  : |$xargs_file_counter
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $insample_directory      => In sample directory
##          : $job_id_href             => Job id hash {REF}
##          : $lane_href               => The lane info hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $program_info_path       => The program info path
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outsample_directory     => Out sample directory
##          : $referencefile_path      => Human genome reference file path
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $file_path;
    my $insample_directory;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $lane_href;
    my $outsample_directory;
    my $parameter_href;
    my $program_info_path;
    my $program_name;
    my $sample_info_href;
    my $sample_id;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
    my $referencefile_path;
    my $temp_directory;
    my $xargs_file_counter;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href
        },
        file_path          => { strict_type => 1, store => \$file_path },
        insample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$insample_directory
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
        lane_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$lane_href
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir
        },
        outsample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outsample_directory
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        program_name      => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name
        },
        referencefile_path => {
            default =>
              $arg_href->{active_parameter_href}{human_genome_reference},
            strict_type => 1,
            store       => \$referencefile_path
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Cluster qw{ check_max_core_number };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::IO::Files qw{ migrate_files xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::Program::Alignment::Picardtools qw{ picardtools_mergesamfiles };
    use MIP::Program::Alignment::Sambamba qw{ split_and_index_aligment_file };
    use MIP::Program::Alignment::Samtools qw{ samtools_index };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_merged_infile_prefix };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $xargs_file_path_prefix;
    my ( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    # Extract lanes
    my $lanes_id = join $EMPTY_STRING, @{ $lane_href->{$sample_id} };

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $sample_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            process_time                    => $time,
            program_directory               => $outaligner_dir,
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    # Used downstream
    $parameter_href->{$mip_program_name}{$sample_id}{indirectory} =
      $outsample_directory;

    ## Assign file_tags
    my $infile_tag = $file_info_href->{$sample_id}
      { $parameter_href->{active_aligner} }{file_tag};
    my $outfile_tag =
      $file_info_href->{$sample_id}{ppicardtools_mergesamfiles}{file_tag};

    ## Assign suffix
    my $infile_suffix = my $outfile_suffix = get_file_suffix(
        {
            jobid_chain    => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
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

    ## Copies files from source to destination
    migrate_files(
        {
            core_number  => $core_number,
            FILEHANDLE   => $FILEHANDLE,
            file_ending  => $infile_tag . $infile_suffix . $ASTERIX,
            indirectory  => $insample_directory,
            infiles_ref  => \@{ $infile_lane_prefix_href->{$sample_id} },
            outfile_path => $temp_directory,
        }
    );

  INFILE:
    foreach my $infile ( @{ $infile_lane_prefix_href->{$sample_id} } ) {

        ## Split BAMs using Samtools
        say {$FILEHANDLE} q{## Split alignment files per contig};
        ($xargs_file_counter) = split_and_index_aligment_file(
            {
                active_parameter_href => $active_parameter_href,
                contigs_ref   => \@{ $file_info_href->{contigs_size_ordered} },
                core_number   => $core_number,
                FILEHANDLE    => $FILEHANDLE,
                file_path     => $file_path,
                infile        => $infile . $infile_tag,
                output_format => substr( $infile_suffix, 1 ),
                program_info_path  => $program_info_path,
                temp_directory     => $temp_directory,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                xargs_file_counter => $xargs_file_counter,
            }
        );
    }

    ## More than one file - we have something to merge
    if ( scalar @{ $infile_lane_prefix_href->{$sample_id} } > 1 ) {

        ## picardtools_mergesamfiles
        say {$FILEHANDLE} q{## Merging alignment files};

        Readonly my $JAVA_MEMORY_ALLOCATION => 4;

        # Division by X according to java Heap size
        $core_number = floor( $active_parameter_href->{node_ram_memory} /
              $JAVA_MEMORY_ALLOCATION );

        ## Limit number of cores requested to the maximum number of cores available per node
        $core_number = check_max_core_number(
            {
                core_number_requested => $core_number,
                max_cores_per_node =>
                  $active_parameter_href->{max_cores_per_node},
            }
        );

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                core_number   => $core_number,
                FILEHANDLE    => $FILEHANDLE,
                file_path     => $file_path,
                first_command => q{java},
                java_jar      => catfile(
                    $active_parameter_href->{picardtools_path},
                    q{picard.jar}
                ),
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                memory_allocation  => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
                program_info_path  => $program_info_path,
                temp_directory     => $temp_directory,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                xargs_file_counter => $xargs_file_counter,
            }
        );

      CONTIG:
        foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

            ## Get parameters
            # Assemble infile paths by adding directory and file ending
            my @infile_paths = map {
                    catfile( $temp_directory, $_ )
                  . $infile_tag
                  . $UNDERSCORE
                  . $contig
                  . $infile_suffix
            } @{ $infile_lane_prefix_href->{$sample_id} };

            my $outfile_path = catfile( $temp_directory,
                    $sample_id
                  . $UNDERSCORE
                  . q{lanes}
                  . $UNDERSCORE
                  . $lanes_id
                  . $outfile_tag
                  . $UNDERSCORE
                  . $contig
                  . $outfile_suffix );
            my $stderrfile_path =
              $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};

            picardtools_mergesamfiles(
                {
                    FILEHANDLE         => $XARGSFILEHANDLE,
                    create_index       => q{true},
                    infile_paths_ref   => \@infile_paths,
                    outfile_path       => $outfile_path,
                    referencefile_path => $referencefile_path,
                    stderrfile_path    => $stderrfile_path,
                    threading          => q{true},
                }
            );
            say {$XARGSFILEHANDLE} $NEWLINE;
        }
    }
    else {
        ## Only 1 infile - rename sample and index instead of merge to streamline handling of filenames downstream

        ## Rename samples
        say {$FILEHANDLE}
q{## Renaming sample instead of merge to streamline handling of filenames downstream};

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
        foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

          INFILES:
            foreach my $infile ( @{ $infile_lane_prefix_href->{$sample_id} } ) {

                ## Get parameters
                my $gnu_mv_infile_path = catfile( $temp_directory,
                        $infile
                      . $infile_tag
                      . $UNDERSCORE
                      . $contig
                      . $infile_suffix );
                my $outfile_name =
                    $sample_id
                  . $UNDERSCORE
                  . q{lanes}
                  . $UNDERSCORE
                  . $lanes_id
                  . $outfile_tag
                  . $UNDERSCORE
                  . $contig
                  . $outfile_suffix;

                my $gnu_mv_outfile_path =
                  catfile( $temp_directory, $outfile_name );
                ## Rename
                gnu_mv(
                    {
                        FILEHANDLE   => $XARGSFILEHANDLE,
                        infile_path  => $gnu_mv_infile_path,
                        outfile_path => $gnu_mv_outfile_path,
                    }
                );
                print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

                ## Index
                samtools_index(
                    {
                        bai_format  => 1,
                        FILEHANDLE  => $XARGSFILEHANDLE,
                        infile_path => $gnu_mv_outfile_path,
                    }
                );
            }
            say {$XARGSFILEHANDLE} $NEWLINE;
        }
    }

    ## Copies file from temporary directory. Per contig
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    ($xargs_file_counter) = xargs_migrate_contig_files(
        {
            contigs_ref  => \@{ $file_info_href->{contigs_size_ordered} },
            core_number  => $core_number,
            FILEHANDLE   => $FILEHANDLE,
            file_ending  => $outfile_suffix . $ASTERIX,
            file_path    => $file_path,
            outdirectory => $outsample_directory,
            outfile      => $sample_id
              . $UNDERSCORE
              . q{lanes}
              . $UNDERSCORE
              . $lanes_id
              . $outfile_tag,
            program_info_path  => $program_info_path,
            temp_directory     => $temp_directory,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    close $XARGSFILEHANDLE;
    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {

        my $qc_outfile_path = catfile( $outsample_directory,
                $sample_id
              . $UNDERSCORE
              . q{lanes}
              . $UNDERSCORE
              . $lanes_id
              . $outfile_tag
              . $UNDERSCORE
              . $file_info_href->{contigs_size_ordered}[0]
              . $outfile_suffix );
        add_program_outfile_to_sample_info(
            {
                infile           => $merged_infile_prefix,
                path             => $qc_outfile_path,
                program_name     => $program_name,
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

sub analysis_picardtools_mergesamfiles_rio {

## Function : Merges all bam files using Picardtools mergesamfiles within each sampleid and files generated previously (option if provided with '-picardtools_mergesamfiles_previous_bams'). The merged files have to be sorted before attempting to merge.
## Returns  : |$xargs_file_counter
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $FILEHANDLE              => Filehandle to write to
##          : $file_info_href          => The file_info hash {REF}
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $insample_directory      => In sample directory
##          : $job_id_href             => Job id hash {REF}
##          : $lane_href               => The lane info hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $referencefile_path      => Human genome reference file path
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $program_info_path       => The program info path
##          : $program_name            => Program name
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $FILEHANDLE;
    my $file_info_href;
    my $file_path;
    my $infile_lane_prefix_href;
    my $insample_directory;
    my $job_id_href;
    my $lane_href;
    my $parameter_href;
    my $program_info_path;
    my $program_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
    my $referencefile_path;
    my $temp_directory;
    my $xargs_file_counter;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id
        },
        FILEHANDLE     => { store => \$FILEHANDLE },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href
        },
        file_path               => { strict_type => 1, store => \$file_path },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href
        },
        insample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$insample_directory
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        lane_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$lane_href
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        program_name      => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name
        },
        referencefile_path => {
            default =>
              $arg_href->{active_parameter_href}{human_genome_reference},
            strict_type => 1,
            store       => \$referencefile_path
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Cluster qw{ check_max_core_number };
    use MIP::Delete::File qw{ delete_files };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::IO::Files qw{ migrate_files xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::Program::Alignment::Picardtools qw{ picardtools_mergesamfiles };
    use MIP::Program::Alignment::Sambamba qw{ split_and_index_aligment_file };
    use MIP::Program::Alignment::Samtools qw{ samtools_index };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Set::File qw{ set_merged_infile_prefix };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $job_id_chain  = $parameter_href->{$mip_program_name}{chain};
    my $reduce_io_ref = \$active_parameter_href->{reduce_io};
    my $xargs_file_path_prefix;
    my ( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    # Extract lanes
    my $lanes_id = join $EMPTY_STRING, @{ $lane_href->{$sample_id} };

    ## Filehandles
    # Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Assign file_tags
    my $infile_tag = $file_info_href->{$sample_id}
      { $parameter_href->{active_aligner} }{file_tag};
    my $outfile_tag =
      $file_info_href->{$sample_id}{ppicardtools_mergesamfiles}{file_tag};

    ## Assign suffix
    my $infile_suffix = my $outfile_suffix = get_file_suffix(
        {
            jobid_chain    => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
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

    ## Copies files from source to destination
    migrate_files(
        {
            core_number  => $core_number,
            FILEHANDLE   => $FILEHANDLE,
            file_ending  => $infile_tag . $infile_suffix . $ASTERIX,
            indirectory  => $insample_directory,
            infiles_ref  => \@{ $infile_lane_prefix_href->{$sample_id} },
            outfile_path => $temp_directory,
        }
    );

  INFILE:
    foreach my $infile ( @{ $infile_lane_prefix_href->{$sample_id} } ) {

        ## Split BAMs using Samtools
        say {$FILEHANDLE} q{## Split alignment files per contig};
        ($xargs_file_counter) = split_and_index_aligment_file(
            {
                active_parameter_href => $active_parameter_href,
                contigs_ref   => \@{ $file_info_href->{contigs_size_ordered} },
                core_number   => $core_number,
                FILEHANDLE    => $FILEHANDLE,
                file_path     => $file_path,
                infile        => $infile . $infile_tag,
                output_format => substr( $infile_suffix, 1 ),
                program_info_path  => $program_info_path,
                temp_directory     => $temp_directory,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                xargs_file_counter => $xargs_file_counter,
            }
        );
    }

    ## More than one file - we have something to merge
    if ( scalar @{ $infile_lane_prefix_href->{$sample_id} } > 1 ) {

        ## picardtools_mergesamfiles
        say {$FILEHANDLE} q{## Merging alignment files};

        Readonly my $JAVA_MEMORY_ALLOCATION => 4;

        # Division by X according to java Heap size
        $core_number = floor( $active_parameter_href->{node_ram_memory} /
              $JAVA_MEMORY_ALLOCATION );

        ## Limit number of cores requested to the maximum number of cores available per node
        $core_number = check_max_core_number(
            {
                core_number_requested => $core_number,
                max_cores_per_node =>
                  $active_parameter_href->{max_cores_per_node},
            }
        );

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                core_number   => $core_number,
                FILEHANDLE    => $FILEHANDLE,
                file_path     => $file_path,
                first_command => q{java},
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                java_jar => catfile(
                    $active_parameter_href->{picardtools_path},
                    q{picard.jar}
                ),
                memory_allocation  => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
                program_info_path  => $program_info_path,
                temp_directory     => $temp_directory,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                xargs_file_counter => $xargs_file_counter,
            }
        );

      CONTIG:
        foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

            ## Get parameters
            # Assemble infile paths by adding directory and file ending
            my @infile_paths = map {
                    catfile( $temp_directory, $_ )
                  . $infile_tag
                  . $UNDERSCORE
                  . $contig
                  . $infile_suffix
            } @{ $infile_lane_prefix_href->{$sample_id} };

            my $outfile_path = catfile( $temp_directory,
                    $sample_id
                  . $UNDERSCORE
                  . q{lanes}
                  . $UNDERSCORE
                  . $lanes_id
                  . $outfile_tag
                  . $UNDERSCORE
                  . $contig
                  . $outfile_suffix );
            my $stderrfile_path =
              $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};

            picardtools_mergesamfiles(
                {
                    create_index       => q{true},
                    FILEHANDLE         => $XARGSFILEHANDLE,
                    infile_paths_ref   => \@infile_paths,
                    outfile_path       => $outfile_path,
                    referencefile_path => $referencefile_path,
                    stderrfile_path    => $stderrfile_path,
                    threading          => q{true},
                }
            );
            say {$XARGSFILEHANDLE} $NEWLINE;
        }
    }
    else {
        ## Only 1 infile - rename sample and index instead of merge to streamline handling of filenames downstream

        ## Rename samples
        say {$FILEHANDLE}
q{## Renaming sample instead of merge to streamline handling of filenames downstream};

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
        foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

          INFILES:
            foreach my $infile ( @{ $infile_lane_prefix_href->{$sample_id} } ) {

                ## Get parameters
                my $gnu_infile_path = catfile( $temp_directory,
                        $infile
                      . $infile_tag
                      . $UNDERSCORE
                      . $contig
                      . $infile_suffix );
                my $outfile_name =
                    $sample_id
                  . $UNDERSCORE
                  . q{lanes}
                  . $UNDERSCORE
                  . $lanes_id
                  . $outfile_tag
                  . $UNDERSCORE
                  . $contig
                  . $outfile_suffix;

                my $gnu_outfile_path =
                  catfile( $temp_directory, $outfile_name );
                ## Rename
                gnu_mv(
                    {
                        FILEHANDLE   => $XARGSFILEHANDLE,
                        infile_path  => $gnu_infile_path,
                        outfile_path => $gnu_outfile_path,
                    }
                );
                print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

                ## Index
                samtools_index(
                    {
                        bai_format  => 1,
                        FILEHANDLE  => $XARGSFILEHANDLE,
                        infile_path => $gnu_outfile_path,
                    }
                );
            }
            say {$XARGSFILEHANDLE} $NEWLINE;
        }
    }

    close $XARGSFILEHANDLE;

    delete_files(
        {
            core_number => $core_number,
            FILEHANDLE  => $FILEHANDLE,
            file_ending => $infile_tag . $ASTERIX,
            indirectory => $temp_directory,
            infiles_ref => \@{ $infile_lane_prefix_href->{$sample_id} },
        }
    );

    # Track the number of created xargs scripts per module
    return $xargs_file_counter;
}

1;
