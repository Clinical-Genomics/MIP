package MIP::Recipes::Analysis::Markduplicates;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile devnull };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
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
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_markduplicates analysis_markduplicates_rio };

}

## Constants
Readonly my $ASTERIX    => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $SEMICOLON  => q{;};
Readonly my $UNDERSCORE => q{_};

sub analysis_markduplicates {

## Function : Mark duplicated reads using Picardtools markduplicates or Sambamba markduplicates in files generated from alignment (sorted, merged).
## Returns  : |$xargs_file_counter
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
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_path               => { store => \$file_path, strict_type => 1, },
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
        program_info_path =>
          { store => \$program_info_path, strict_type => 1, },
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
        xargs_file_counter => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Delete::File qw{ delete_contig_files };
    use MIP::Get::File
      qw{ get_file_suffix get_merged_infile_prefix get_io_files };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Gnu::Coreutils qw{ gnu_cat };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::Program::Alignment::Sambamba
      qw{ sambamba_flagstat sambamba_markdup };
    use MIP::Program::Alignment::Picardtools qw{ picardtools_markduplicates };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_io_files };
    use MIP::QC::Record
      qw{ add_program_metafile_to_sample_info add_program_outfile_to_sample_info };

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

    my $outdir_path_prefix    = $io{out}{dir_path_prefix};
    my $outfile_name_prefix   = $io{out}{file_name_prefix};
    my $temp_file_path_prefix = $io{temp}{file_path_prefix};
    my @temp_outfile_paths    = @{ $io{temp}{file_paths} };

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    # Store which program performed the markduplication
    my $markduplicates_program;

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

    ### SHELL:

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    ($xargs_file_counter) = xargs_migrate_contig_files(
        {
            contigs_ref        => \@{ $file_info_href->{contigs_size_ordered} },
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_ending        => substr( $infile_suffix, 0, 2 ) . $ASTERIX,
            file_path          => $file_path,
            indirectory        => $indir_path_prefix,
            infile             => $infile_name_prefix,
            program_info_path  => $program_info_path,
            temp_directory     => $temp_directory,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Marking Duplicates
    say {$FILEHANDLE} q{## Marking Duplicates};

    ## Picardtools
    if ( $active_parameter_href->{markduplicates_picardtools_markduplicates} ) {

        $markduplicates_program = q{picardtools_markduplicates};

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
                memory_allocation  => q{Xmx4g},
                program_info_path  => $program_info_path,
                temp_directory     => $temp_directory,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                xargs_file_counter => $xargs_file_counter,
            }
        );

      CONTIG:
        while ( my ( $infile_index, $contig ) =
            each @{ $file_info_href->{contigs_size_ordered} } )
        {

            my $stderrfile_path =
              $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
            my $metrics_file =
              $temp_file_path_prefix . $DOT . $contig . $DOT . q{metric};
            picardtools_markduplicates(
                {
                    create_index       => q{true},
                    FILEHANDLE         => $XARGSFILEHANDLE,
                    infile_paths_ref   => [ $temp_infile_paths[$infile_index] ],
                    metrics_file       => $metrics_file,
                    outfile_path       => $temp_outfile_paths[$infile_index],
                    referencefile_path => $referencefile_path,
                    stderrfile_path    => $stderrfile_path,
                }
            );
            print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

            ## Process BAM with sambamba flagstat to produce metric file for downstream analysis
            sambamba_flagstat(
                {
                    FILEHANDLE   => $XARGSFILEHANDLE,
                    infile_path  => $temp_outfile_paths[$infile_index],
                    outfile_path => $temp_file_path_prefix
                      . $DOT
                      . $contig
                      . $UNDERSCORE
                      . q{metric},
                    stderrfile_path => $stderrfile_path,
                }
            );
            say {$XARGSFILEHANDLE} $NEWLINE;
        }
    }

    ## Sambamba
    if ( $active_parameter_href->{markduplicates_sambamba_markdup} ) {

        $markduplicates_program = q{sambamba_markdup};

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

            my $stderrfile_path =
              $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
            sambamba_markdup(
                {
                    FILEHANDLE      => $XARGSFILEHANDLE,
                    hash_table_size => $active_parameter_href
                      ->{markduplicates_sambamba_markdup_hash_table_size},
                    infile_path    => [ $temp_infile_paths[$infile_index] ],
                    io_buffer_size => $active_parameter_href
                      ->{markduplicates_sambamba_markdup_io_buffer_size},
                    outfile_path       => $temp_outfile_paths[$infile_index],
                    overflow_list_size => $active_parameter_href
                      ->{markduplicates_sambamba_markdup_overflow_list_size},
                    temp_directory  => $temp_directory,
                    show_progress   => 1,
                    stderrfile_path => $stderrfile_path,
                }
            );
            print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

            ## Process BAM with sambamba flagstat to produce metric file for downstream analysis
            sambamba_flagstat(
                {
                    FILEHANDLE   => $XARGSFILEHANDLE,
                    infile_path  => $temp_outfile_paths[$infile_index],
                    outfile_path => $temp_file_path_prefix
                      . $DOT
                      . $contig
                      . $UNDERSCORE
                      . q{metric},
                    stderrfile_path => $stderrfile_path,
                }
            );
            say {$XARGSFILEHANDLE} $NEWLINE;
        }
    }

    ## Concatenate all metric files
    gnu_cat(
        {
            FILEHANDLE       => $FILEHANDLE,
            infile_paths_ref => [
                    $temp_file_path_prefix
                  . $DOT
                  . $ASTERIX
                  . $UNDERSCORE
                  . q{metric}
            ],
            outfile_path => $temp_file_path_prefix
              . $UNDERSCORE
              . q{metric_all},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Collect duplicate reads and reads mapped across all metric contig files. Calculate fraction duplicates.
    ## Write it to stdout.
    _calculate_fraction_duplicates_for_all_metric_files(
        {
            FILEHANDLE          => $FILEHANDLE,
            outfile_path_prefix => $temp_file_path_prefix,
        }
    );

    migrate_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $temp_file_path_prefix . $UNDERSCORE . q{metric},
            outfile_path => $outdir_path_prefix,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Copies file from temporary directory. Per contig
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    ($xargs_file_counter) = xargs_migrate_contig_files(
        {
            contigs_ref        => \@{ $file_info_href->{contigs_size_ordered} },
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $file_path,
            file_ending        => substr( $outfile_suffix, 0, 2 ) . $ASTERIX,
            outdirectory       => $outdir_path_prefix,
            outfile            => $outfile_name_prefix,
            program_info_path  => $program_info_path,
            temp_directory     => $temp_directory,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Close FILEHANDLES
    close $XARGSFILEHANDLE;
    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                infile => $outfile_name_prefix,
                path   => catfile(
                    $outdir_path_prefix,
                    $outfile_name_prefix . $UNDERSCORE . q{metric}
                ),
                program_name     => q{markduplicates},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

# Markduplicates can be processed by either picardtools markduplicates or sambamba markdup
        add_program_metafile_to_sample_info(
            {
                infile           => $outfile_name_prefix,
                metafile_tag     => q{marking_duplicates},
                processed_by     => $markduplicates_program,
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

sub analysis_markduplicates_rio {

## Function : Mark duplicated reads using Picardtools markduplicates or Sambamba markduplicates in files generated from alignment (sorted, merged).
## Returns  : |$xargs_file_counter
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $FILEHANDLE              => Filehandle to write to
##          : $file_info_href          => File info hash {REF}
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $parameter_href          => Parameter hash {REF}
##          : $program_info_path       => The program info path
##          : $program_name            => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $active_parameter_href;
    my $file_info_href;
    my $file_path;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_info_path;
    my $program_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
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
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
            strict_type => 1,
        },
        FILEHANDLE     => { required => 1, store => \$FILEHANDLE, },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_path               => { store => \$file_path, strict_type => 1, },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        program_info_path =>
          { store => \$program_info_path, strict_type => 1, },
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
        xargs_file_counter => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Delete::File qw{ delete_contig_files };
    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Gnu::Coreutils qw{ gnu_cat };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::Program::Alignment::Sambamba
      qw{ sambamba_flagstat sambamba_markdup };
    use MIP::Program::Alignment::Picardtools qw{ picardtools_markduplicates };
    use MIP::QC::Record
      qw{ add_program_metafile_to_sample_info add_program_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $program_mode = $active_parameter_href->{$program_name};

    ## Unpack parameters
    my $job_id_chain       = $parameter_href->{$program_name}{chain};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $xargs_file_path_prefix;
    my ($core_number) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Assign directories
    my $outsample_directory = catdir( $active_parameter_href->{outdata_dir},
        $sample_id, $outaligner_dir );

    # Used downstream
    $parameter_href->{$program_name}{$sample_id}{indirectory} =
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
      $file_info_href->{$sample_id}{picardtools_mergesamfiles}{file_tag};
    my $outfile_tag =
      $file_info_href->{$sample_id}{$program_name}{file_tag};

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

    # Store which program performed the markduplication
    my $markduplicates_program;

    ## Marking Duplicates
    say {$FILEHANDLE} q{## Marking Duplicates};

    ##Picardtools
    if ( $active_parameter_href->{markduplicates_picardtools_markduplicates} ) {

        $markduplicates_program = q{picardtools_markduplicates};

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
                memory_allocation  => q{Xmx4g},
                program_info_path  => $program_info_path,
                temp_directory     => $temp_directory,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                xargs_file_counter => $xargs_file_counter,
            }
        );

      CONTIG:
        foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

            my $outfile_path =
              $outfile_path_prefix . $UNDERSCORE . $contig . $outfile_suffix;
            my $stderrfile_path =
              $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
            my $metrics_file =
              $outfile_path_prefix . $UNDERSCORE . $contig . $DOT . q{metric};
            picardtools_markduplicates(
                {
                    create_index     => q{true},
                    FILEHANDLE       => $XARGSFILEHANDLE,
                    infile_paths_ref => [
                            $file_path_prefix
                          . $UNDERSCORE
                          . $contig
                          . $infile_suffix
                    ],
                    metrics_file       => $metrics_file,
                    outfile_path       => $outfile_path,
                    referencefile_path => $referencefile_path,
                    stderrfile_path    => $stderrfile_path,
                }
            );
            print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

            ## Process BAM with sambamba flagstat to produce metric file for downstream analysis
            sambamba_flagstat(
                {
                    FILEHANDLE   => $XARGSFILEHANDLE,
                    infile_path  => $outfile_path,
                    outfile_path => $outfile_path_prefix
                      . $UNDERSCORE
                      . $contig
                      . $UNDERSCORE
                      . q{metric},
                    stderrfile_path => $stderrfile_path,
                }
            );
            say {$XARGSFILEHANDLE} $NEWLINE;
        }
    }

    ## Sambamba
    if ( $active_parameter_href->{markduplicates_sambamba_markdup} ) {

        $markduplicates_program = q{sambamba_markdup};

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

            my $infile_path =
              $file_path_prefix . $UNDERSCORE . $contig . $infile_suffix;
            my $outfile_path =
              $outfile_path_prefix . $UNDERSCORE . $contig . $outfile_suffix;
            my $stderrfile_path =
              $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
            sambamba_markdup(
                {
                    FILEHANDLE      => $XARGSFILEHANDLE,
                    hash_table_size => $active_parameter_href
                      ->{markduplicates_sambamba_markdup_hash_table_size},
                    infile_path    => $infile_path,
                    io_buffer_size => $active_parameter_href
                      ->{markduplicates_sambamba_markdup_io_buffer_size},
                    outfile_path       => $outfile_path,
                    overflow_list_size => $active_parameter_href
                      ->{markduplicates_sambamba_markdup_overflow_list_size},
                    show_progress   => 1,
                    stderrfile_path => $stderrfile_path,
                    temp_directory  => $temp_directory,
                }
            );
            print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

            ## Process BAM with sambamba flagstat to produce metric file for downstream analysis
            sambamba_flagstat(
                {
                    FILEHANDLE   => $XARGSFILEHANDLE,
                    infile_path  => $outfile_path,
                    outfile_path => $outfile_path_prefix
                      . $UNDERSCORE
                      . $contig
                      . $UNDERSCORE
                      . q{metric},
                    stderrfile_path => $stderrfile_path,
                }
            );
            say {$XARGSFILEHANDLE} $NEWLINE;
        }
    }

    ## Concatenate all metric files
    gnu_cat(
        {
            FILEHANDLE       => $FILEHANDLE,
            infile_paths_ref => [
                    $outfile_path_prefix
                  . $UNDERSCORE
                  . $ASTERIX
                  . $UNDERSCORE
                  . q{metric}
            ],
            outfile_path => $outfile_path_prefix . $UNDERSCORE . q{metric_all},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Collect duplicate reads and reads mapped across all metric contig files. Calculate fraction duplicates.
    ## Write it to stdout.
    _calculate_fraction_duplicates_for_all_metric_files(
        {
            FILEHANDLE          => $FILEHANDLE,
            outfile_path_prefix => $outfile_path_prefix,
        }
    );

    migrate_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $outfile_path_prefix . $UNDERSCORE . q{metric},
            outfile_path => $outsample_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Remove file at temporary Directory
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

    close $XARGSFILEHANDLE;

    if ( $program_mode == 1 ) {

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                infile => $merged_infile_prefix,
                path   => catfile(
                    $outsample_directory,
                    $outfile_prefix . $UNDERSCORE . q{metric}
                ),
                program_name     => q{markduplicates},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

# Markduplicates can be processed by either picardtools markduplicates or sambamba markdup
        add_program_metafile_to_sample_info(
            {
                infile           => $merged_infile_prefix,
                metafile_tag     => q{marking_duplicates},
                processed_by     => $markduplicates_program,
                program_name     => $program_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );
    }

    # Track the number of created xargs scripts per module for Block algorithm
    return $xargs_file_counter;
}

sub _calculate_fraction_duplicates_for_all_metric_files {

## Function : Collect duplicate reads and reads mapped across all metric contig files. Calculate fraction duplicates. Write it to stdout.
## Returns  :
## Arguments: $FILEHANDLE          => Filehandle to write to
##          : $outfile_path_prefix => Outfile path prefix

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $outfile_path_prefix;

    my $tmpl = {
        FILEHANDLE => { defined => 1, required => 1, store => \$FILEHANDLE, },
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
    $regexp .=
q?print "Read Mapped: ".$feature{map}."\nDuplicates: ".$feature{dup}."\n".?;

    # Print Fraction duplicates to stdout
    $regexp .= q?"Fraction Duplicates: ".$feature{dup}/$feature{map}, "\n"; ?;

    # Quit
    $regexp .= q?last;'?;

    ## Sum metric over concatenated file
    print {$FILEHANDLE} $regexp . $SPACE;
    print {$FILEHANDLE} $outfile_path_prefix
      . $UNDERSCORE
      . q{metric_all}
      . $SPACE;
    say {$FILEHANDLE} q{>}
      . $SPACE
      . $outfile_path_prefix
      . $UNDERSCORE
      . q{metric}
      . $SPACE,
      $NEWLINE;
    return;
}

1;
