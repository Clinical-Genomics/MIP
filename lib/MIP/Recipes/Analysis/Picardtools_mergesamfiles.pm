package MIP::Recipes::Analysis::Picardtools_mergesamfiles;

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
use File::Basename qw{ fileparse };
use File::Spec::Functions qw{ catdir catfile };
use POSIX;

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ analysis_picardtools_mergesamfiles analysis_picardtools_mergesamfiles_rio };

}

##Constants
Readonly my $ASTERIX      => q{*};
Readonly my $DOT          => q{.};
Readonly my $EMPTY_STRING => q{};
Readonly my $NEWLINE      => qq{\n};
Readonly my $SEMICOLON    => q{;};
Readonly my $SPACE        => q{ };
Readonly my $UNDERSCORE   => q{_};

sub analysis_picardtools_mergesamfiles {

##analysis_picardtools_mergesamfiles

##Function : Merges all bam files using Picardtools mergesamfiles within each sampleid and files generated previously (option if provided with '-picardtools_mergesamfiles_previous_bams'). The merged files have to be sorted before attempting to merge.
##Returns  : |$xargs_file_counter
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $lane_href, $job_id_href, $insample_directory, $outsample_directory, $sample_id, $program_name, $program_info_path, $file_path, $family_id, $outaligner_dir, $referencefile_path, $xargs_file_counter
##         : $parameter_href          => Parameter hash {REF}
##         : $active_parameter_href   => Active parameters for this analysis hash {REF}
##         : $sample_info_href        => Info on samples and family hash {REF}
##         : $file_info_href          => File_info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $lane_href               => The lane info hash {REF}
##         : $job_id_href             => Job id hash {REF}
##         : $insample_directory      => In sample directory
##         : $outsample_directory     => Out sample directory
##         : $sample_id               => Sample id
##         : $program_name            => Program name
##         : $program_info_path       => The program info path
##         : $file_path               => File path
##         : $family_id               => Family id
##         : $temp_directory          => Temporary directory
##         : $outaligner_dir          => Outaligner_dir used in the analysis
##         : $referencefile_path      => Human genome reference file path
##         : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id;
    my $temp_directory;
    my $outaligner_dir;
    my $referencefile_path;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $lane_href;
    my $job_id_href;
    my $insample_directory;
    my $outsample_directory;
    my $sample_id;
    my $program_name;
    my $program_info_path;
    my $file_path;

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
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href
        },
        lane_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$lane_href
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
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
            store       => \$insample_directory
        },
        outsample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outsample_directory
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        file_path         => { strict_type => 1, store => \$file_path },
        family_id         => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir
        },
        referencefile_path => {
            default =>
              $arg_href->{active_parameter_href}{human_genome_reference},
            strict_type => 1,
            store       => \$referencefile_path
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
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::IO::Files qw{ migrate_files xargs_migrate_contig_files };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Program::Alignment::Sambamba qw{ split_and_index_aligment_file };
    use MIP::Program::Alignment::Picardtools qw{ picardtools_mergesamfiles };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Program::Alignment::Samtools qw{ samtools_index };
    use MIP::Set::File qw{ set_merged_infile_prefix };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };

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
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $xargs_file_path_prefix;

    # Extract lanes
    my $lanes_id = join $EMPTY_STRING, @{ $lane_href->{$sample_id} };

    ## Filehandles
    # Create anonymous filehandle
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
            program_directory     => $outaligner_dir,
            core_number           => $core_number,
            process_time          => $time,
            temp_directory        => $temp_directory,
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
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
            jobid_chain    => $job_id_chain,
        }
    );

    ## Set helper value for finding merged_infiles downstream
    my $merged_infile_prefix =
      $sample_id . $UNDERSCORE . q{lanes} . $UNDERSCORE . $lanes_id;
    set_merged_infile_prefix(
        {
            file_info_href       => $file_info_href,
            sample_id            => $sample_id,
            merged_infile_prefix => $merged_infile_prefix,
        }
    );

    ## Copies files from source to destination
    migrate_files(
        {
            infiles_ref  => \@{ $infile_lane_prefix_href->{$sample_id} },
            outfile_path => $temp_directory,
            FILEHANDLE   => $FILEHANDLE,
            indirectory  => $insample_directory,
            core_number  => $core_number,
            file_ending  => $infile_tag . $infile_suffix . $ASTERIX,
        }
    );

  INFILE:
    foreach my $infile ( @{ $infile_lane_prefix_href->{$sample_id} } ) {

        ## Split BAMs using Samtools
        say {$FILEHANDLE} q{## Split alignment files per contig};
        ($xargs_file_counter) = split_and_index_aligment_file(
            {
                active_parameter_href => $active_parameter_href,
                contigs_ref => \@{ $file_info_href->{contigs_size_ordered} },
                FILEHANDLE  => $FILEHANDLE,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                file_path          => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
                xargs_file_counter => $xargs_file_counter,
                temp_directory     => $temp_directory,
                infile             => $infile . $infile_tag,
                output_format      => substr( $infile_suffix, 1 ),
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
                max_cores_per_node =>
                  $active_parameter_href->{max_cores_per_node},
                core_number_requested => $core_number,
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
                first_command      => q{java},
                xargs_file_counter => $xargs_file_counter,
                memory_allocation  => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                temp_directory => $temp_directory,
                java_jar       => catfile(
                    $active_parameter_href->{picardtools_path},
                    q{picard.jar}
                ),
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
                    infile_paths_ref   => \@infile_paths,
                    outfile_path       => $outfile_path,
                    stderrfile_path    => $stderrfile_path,
                    threading          => q{true},
                    create_index       => q{true},
                    referencefile_path => $referencefile_path,
                    FILEHANDLE         => $XARGSFILEHANDLE,
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
                FILEHANDLE         => $FILEHANDLE,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                file_path          => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
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
                        infile_path  => $gnu_mv_infile_path,
                        outfile_path => $gnu_mv_outfile_path,
                        FILEHANDLE   => $XARGSFILEHANDLE,
                    }
                );
                print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

                ## Index
                samtools_index(
                    {
                        infile_path => $gnu_mv_outfile_path,
                        bai_format  => 1,
                        FILEHANDLE  => $XARGSFILEHANDLE,
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
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            contigs_ref        => \@{ $file_info_href->{contigs_size_ordered} },
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
            outfile            => $sample_id
              . $UNDERSCORE
              . q{lanes}
              . $UNDERSCORE
              . $lanes_id
              . $outfile_tag,
            outdirectory   => $outsample_directory,
            temp_directory => $temp_directory,
            file_ending    => $outfile_suffix . $ASTERIX,
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
              . $UNDERSCORE
              . $file_info_href->{contigs_size_ordered}[0]
              . $outfile_tag
              . $outfile_suffix );
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                sample_id        => $sample_id,
                program_name     => $program_name,
                infile           => $merged_infile_prefix,
                path             => $qc_outfile_path,
            }
        );

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

sub analysis_picardtools_mergesamfiles_rio {

##analysis_picardtools_mergesamfiles_rio

##Function : Merges all bam files using Picardtools mergesamfiles within each sampleid and files generated previously (option if provided with '-picardtools_mergesamfiles_previous_bams'). The merged files have to be sorted before attempting to merge.
##Returns  : |$xargs_file_counter
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $lane_href, $job_id_href, $insample_directory, $sample_id, $program_name, $program_info_path, $file_path,, $FILEHANDLE, $family_id, $outaligner_dir, $referencefile_path, $xargs_file_counter
##         : $parameter_href          => Parameter hash {REF}
##         : $active_parameter_href   => Active parameters for this analysis hash {REF}
##         : $sample_info_href        => Info on samples and family hash {REF}
##         : $file_info_href          => The file_info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $lane_href               => The lane info hash {REF}
##         : $job_id_href             => Job id hash {REF}
##         : $insample_directory      => In sample directory
##         : $sample_id               => Sample id
##         : $program_name            => Program name
##         : $program_info_path       => The program info path
##         : $file_path               => File path
##         : $FILEHANDLE              => Filehandle to write to
##         : $family_id               => Family id
##         : $temp_directory          => Temporary directory
##         : $outaligner_dir          => Outaligner_dir used in the analysis
##         : $referencefile_path      => Human genome reference file path
##         : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id;
    my $temp_directory;
    my $outaligner_dir;
    my $referencefile_path;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $lane_href;
    my $job_id_href;
    my $insample_directory;
    my $sample_id;
    my $program_name;
    my $program_info_path;
    my $file_path;
    my $FILEHANDLE;

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
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href
        },
        lane_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$lane_href
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
        program_info_path => { strict_type => 1, store => \$program_info_path },
        file_path         => { strict_type => 1, store => \$file_path },
        FILEHANDLE => { store => \$FILEHANDLE },
        family_id  => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir
        },
        referencefile_path => {
            default =>
              $arg_href->{active_parameter_href}{human_genome_reference},
            strict_type => 1,
            store       => \$referencefile_path
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
    use MIP::IO::Files qw{ migrate_files xargs_migrate_contig_files };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Program::Alignment::Sambamba qw{ split_and_index_aligment_file };
    use MIP::Program::Alignment::Picardtools qw{ picardtools_mergesamfiles };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Program::Alignment::Samtools qw{ samtools_index };
    use MIP::Set::File qw{ set_merged_infile_prefix };
    use MIP::Delete::File qw{ delete_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };

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
    my $reduce_io_ref = \$active_parameter_href->{reduce_io};
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $xargs_file_path_prefix;

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
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
            jobid_chain    => $job_id_chain,
        }
    );

    ## Set helper value for finding merged_infiles downstream
    my $merged_infile_prefix =
      $sample_id . $UNDERSCORE . q{lanes} . $UNDERSCORE . $lanes_id;
    set_merged_infile_prefix(
        {
            file_info_href       => $file_info_href,
            sample_id            => $sample_id,
            merged_infile_prefix => $merged_infile_prefix,
        }
    );

    ## Copies files from source to destination
    migrate_files(
        {
            infiles_ref  => \@{ $infile_lane_prefix_href->{$sample_id} },
            outfile_path => $temp_directory,
            FILEHANDLE   => $FILEHANDLE,
            indirectory  => $insample_directory,
            core_number  => $core_number,
            file_ending  => $infile_tag . $infile_suffix . $ASTERIX,
        }
    );

  INFILE:
    foreach my $infile ( @{ $infile_lane_prefix_href->{$sample_id} } ) {

        ## Split BAMs using Samtools
        say {$FILEHANDLE} q{## Split alignment files per contig};
        ($xargs_file_counter) = split_and_index_aligment_file(
            {
                active_parameter_href => $active_parameter_href,
                contigs_ref => \@{ $file_info_href->{contigs_size_ordered} },
                FILEHANDLE  => $FILEHANDLE,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                file_path          => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
                xargs_file_counter => $xargs_file_counter,
                temp_directory     => $temp_directory,
                infile             => $infile . $infile_tag,
                output_format      => substr( $infile_suffix, 1 ),
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
                max_cores_per_node =>
                  $active_parameter_href->{max_cores_per_node},
                core_number_requested => $core_number,
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
                first_command      => q{java},
                xargs_file_counter => $xargs_file_counter,
                memory_allocation  => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                temp_directory => $temp_directory,
                java_jar       => catfile(
                    $active_parameter_href->{picardtools_path},
                    q{picard.jar}
                ),
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
                    infile_paths_ref   => \@infile_paths,
                    outfile_path       => $outfile_path,
                    stderrfile_path    => $stderrfile_path,
                    threading          => q{true},
                    create_index       => q{true},
                    referencefile_path => $referencefile_path,
                    FILEHANDLE         => $XARGSFILEHANDLE,
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
                FILEHANDLE         => $FILEHANDLE,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                file_path          => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
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
                        infile_path  => $gnu_infile_path,
                        outfile_path => $gnu_outfile_path,
                        FILEHANDLE   => $XARGSFILEHANDLE,
                    }
                );
                print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

                ## Index
                samtools_index(
                    {
                        infile_path => $gnu_outfile_path,
                        bai_format  => 1,
                        FILEHANDLE  => $XARGSFILEHANDLE,
                    }
                );
            }
            say {$XARGSFILEHANDLE} $NEWLINE;
        }
    }

    close $XARGSFILEHANDLE;

    delete_files(
        {
            infiles_ref => \@{ $infile_lane_prefix_href->{$sample_id} },
            FILEHANDLE  => $FILEHANDLE,
            indirectory => $temp_directory,
            core_number => $core_number,
            file_ending => $infile_tag . $ASTERIX,
        }
    );

    # Track the number of created xargs scripts per module
    return $xargs_file_counter;
}

1;
