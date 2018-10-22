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
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ analysis_picardtools_mergesamfiles analysis_picardtools_mergesamfiles_rio };

}

## Constants
Readonly my $ASTERISK     => q{*};
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
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $program_info_path       => The program info path
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
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_info_path;
    my $program_name;
    my $sample_info_href;
    my $sample_id;

    ## Default(s)
    my $family_id;
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
        program_info_path => { store => \$program_info_path, strict_type => 1, },
        program_name      => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
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

    use MIP::Check::Cluster qw{ check_max_core_number };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_module_parameters  get_program_attributes };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::IO::Files qw{ migrate_files xargs_migrate_contig_files };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Alignment::Picardtools qw{ picardtools_mergesamfiles };
    use MIP::Program::Alignment::Sambamba qw{ split_and_index_aligment_file };
    use MIP::Program::Alignment::Samtools qw{ samtools_index };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_merged_infile_prefix };

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
    my $indir_path_prefix         = $io{in}{dir_path_prefix};
    my $infile_suffix             = $io{in}{file_suffix};
    my @infile_name_prefixes      = @{ $io{in}{file_name_prefixes} };
    my @temp_infile_path_prefixes = @{ $io{temp}{file_path_prefixes} };

    my %prg_atr = get_program_attributes(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
        }
    );
    my $job_id_chain = $prg_atr{chain};
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $program_mode = $active_parameter_href->{$program_name};
    my $xargs_file_path_prefix;
    my ( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
    );

    ## Assign suffix
    my $outfile_suffix = $prg_atr{outfile_suffix};

    ## Extract lanes
    my $lanes_id = join $EMPTY_STRING, @{ $file_info_href->{$sample_id}{lanes} };

    ## Outpaths
    my $outsample_directory =
      catdir( $active_parameter_href->{outdata_dir}, $sample_id, $program_name );
    my $outfile_tag =
      $file_info_href->{$sample_id}{$program_name}{file_tag};
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

    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_name_prefix = $io{out}{file_name_prefix};
    @outfile_paths = @{ $io{out}{file_paths} };
    my @outfile_suffixes  = @{ $io{out}{file_suffixes} };
    my %temp_outfile_path = %{ $io{temp}{file_path_href} };

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
            log                             => $log,
            process_time                    => $time,
            program_directory               => $program_name,
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
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

    ## Copies files from source to destination
    migrate_files(
        {
            core_number  => $core_number,
            FILEHANDLE   => $FILEHANDLE,
            file_ending  => substr( $infile_suffix, 0, 2 ) . $ASTERISK,
            indirectory  => $indir_path_prefix,
            infiles_ref  => \@infile_name_prefixes,
            outfile_path => $temp_directory,
        }
    );

  INFILE:
    foreach my $infile (@infile_name_prefixes) {

        ## Split BAMs
        say {$FILEHANDLE} q{## Split alignment files per contig};
        ($xargs_file_counter) = split_and_index_aligment_file(
            {
                active_parameter_href => $active_parameter_href,
                contigs_ref           => \@{ $file_info_href->{contigs_size_ordered} },
                core_number           => $core_number,
                FILEHANDLE            => $FILEHANDLE,
                file_path             => $file_path,
                infile                => $infile,
                output_format         => substr( $infile_suffix, 1 ),
                program_info_path     => $program_info_path,
                temp_directory        => $temp_directory,
                XARGSFILEHANDLE       => $XARGSFILEHANDLE,
                xargs_file_counter    => $xargs_file_counter,
            }
        );
    }

    ## More than one file - we have something to merge
    if ( scalar @infile_name_prefixes > 1 ) {

        ## picardtools_mergesamfiles
        say {$FILEHANDLE} q{## Merging alignment files};

        Readonly my $JAVA_MEMORY_ALLOCATION => 4;

        # Division by X according to java Heap size
        $core_number =
          floor( $active_parameter_href->{node_ram_memory} / $JAVA_MEMORY_ALLOCATION );

        ## Limit number of cores requested to the maximum number of cores available per node
        $core_number = check_max_core_number(
            {
                core_number_requested => $core_number,
                max_cores_per_node    => $active_parameter_href->{max_cores_per_node},
            }
        );

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                core_number   => $core_number,
                FILEHANDLE    => $FILEHANDLE,
                file_path     => $file_path,
                first_command => q{java},
                java_jar =>
                  catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                memory_allocation    => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
                program_info_path    => $program_info_path,
                temp_directory       => $temp_directory,
                XARGSFILEHANDLE      => $XARGSFILEHANDLE,
                xargs_file_counter   => $xargs_file_counter,
            }
        );

      CONTIG:
        foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

            ## Get parameters
            # Assemble infile paths by adding directory and file suffixes
            my @merge_temp_infile_paths =
              map { $_ . $DOT . $contig . $infile_suffix } @temp_infile_path_prefixes;
            my $stderrfile_path =
              $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};

            picardtools_mergesamfiles(
                {
                    create_index       => q{true},
                    FILEHANDLE         => $XARGSFILEHANDLE,
                    infile_paths_ref   => \@merge_temp_infile_paths,
                    outfile_path       => $temp_outfile_path{$contig},
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
            foreach my $temp_infile_path_prefix (@temp_infile_path_prefixes) {

                ## Get parameters
                my $gnu_temp_infile_path =
                  $temp_infile_path_prefix . $DOT . $contig . $infile_suffix;
                my $gnu_temp_outfile_path = $temp_outfile_path{$contig};

                ## Rename
                gnu_mv(
                    {
                        FILEHANDLE   => $XARGSFILEHANDLE,
                        infile_path  => $gnu_temp_infile_path,
                        outfile_path => $gnu_temp_outfile_path,
                    }
                );
                print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

                ## Index
                samtools_index(
                    {
                        bai_format  => 1,
                        FILEHANDLE  => $XARGSFILEHANDLE,
                        infile_path => $gnu_temp_outfile_path,
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
            contigs_ref        => \@{ $file_info_href->{contigs_size_ordered} },
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_ending        => substr( $outfile_suffix, 0, 2 ) . $ASTERISK,
            file_path          => $file_path,
            outdirectory       => $outdir_path_prefix,
            outfile            => $outfile_name_prefix,
            program_info_path  => $program_info_path,
            temp_directory     => $temp_directory,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    close $XARGSFILEHANDLE;
    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        my $qc_outfile_path = $outfile_paths[0];
        add_program_outfile_to_sample_info(
            {
                infile           => $merged_infile_prefix,
                path             => $qc_outfile_path,
                program_name     => $program_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                active_parameter_href   => $active_parameter_href,
                dependency_method       => q{sample_to_sample},
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_chain            => $job_id_chain,
                job_id_href             => $job_id_href,
                log                     => $log,
                recipe_file_name        => $file_path,
                sample_id               => $sample_id,
            }
        );
    }
    return;
}

sub analysis_picardtools_mergesamfiles_rio {

## Function :  Merges all bam files using Picardtools mergesamfiles within each sampleid. The merged files have to be sorted before attempting to merge.
## Returns  : |$xargs_file_counter
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $FILEHANDLE              => Filehandle to write to
##          : $file_info_href          => The file_info hash {REF}
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
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
    my $job_id_href;
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
        FILEHANDLE     => { store => \$FILEHANDLE },
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
        program_info_path => { store => \$program_info_path, strict_type => 1, },
        program_name      => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
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

    use MIP::Check::Cluster qw{ check_max_core_number };
    use MIP::Delete::File qw{ delete_files };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::IO::Files qw{ migrate_files xargs_migrate_contig_files };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::Program::Alignment::Picardtools qw{ picardtools_mergesamfiles };
    use MIP::Program::Alignment::Sambamba qw{ split_and_index_aligment_file };
    use MIP::Program::Alignment::Samtools qw{ samtools_index };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Set::File qw{ set_merged_infile_prefix };

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
    my $indir_path_prefix         = $io{in}{dir_path_prefix};
    my $infile_suffix             = $io{in}{file_suffix};
    my @infile_name_prefixes      = @{ $io{in}{file_name_prefixes} };
    my @temp_infile_name_prefixes = @{ $io{temp}{file_name_prefixes} };
    my @temp_infile_path_prefixes = @{ $io{temp}{file_path_prefixes} };

    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my %prg_atr = get_program_attributes(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
        }
    );
    my $job_id_chain  = $prg_atr{chain};
    my $program_mode  = $active_parameter_href->{$program_name};
    my $reduce_io_ref = \$active_parameter_href->{reduce_io};
    my $xargs_file_path_prefix;
    my ( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
    );

    ## Assign suffix
    my $outfile_suffix = $prg_atr{outfile_suffix};

    # Extract lanes
    my $lanes_id = join $EMPTY_STRING, @{ $file_info_href->{$sample_id}{lanes} };

    ## Outpaths
    my $outsample_directory =
      catdir( $active_parameter_href->{outdata_dir}, $sample_id, $program_name );
    my $outfile_tag =
      $file_info_href->{$sample_id}{$program_name}{file_tag};
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

    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_name_prefix = $io{out}{file_name_prefix};
    @outfile_paths = @{ $io{out}{file_paths} };
    my @outfile_suffixes  = @{ $io{out}{file_suffixes} };
    my %temp_outfile_path = %{ $io{temp}{file_path_href} };

    ## Filehandles
    # Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();

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

    ## Copies files from source to destination
    migrate_files(
        {
            core_number  => $core_number,
            FILEHANDLE   => $FILEHANDLE,
            file_ending  => substr( $infile_suffix, 0, 2 ) . $ASTERISK,
            indirectory  => $indir_path_prefix,
            infiles_ref  => \@infile_name_prefixes,
            outfile_path => $temp_directory,
        }
    );

  INFILE:
    foreach my $infile (@infile_name_prefixes) {

        ## Split BAMs using Samtools
        say {$FILEHANDLE} q{## Split alignment files per contig};
        ($xargs_file_counter) = split_and_index_aligment_file(
            {
                active_parameter_href => $active_parameter_href,
                contigs_ref           => \@{ $file_info_href->{contigs_size_ordered} },
                core_number           => $core_number,
                FILEHANDLE            => $FILEHANDLE,
                file_path             => $file_path,
                infile                => $infile,
                output_format         => substr( $infile_suffix, 1 ),
                program_info_path     => $program_info_path,
                temp_directory        => $temp_directory,
                XARGSFILEHANDLE       => $XARGSFILEHANDLE,
                xargs_file_counter    => $xargs_file_counter,
            }
        );
    }

    ## More than one file - we have something to merge
    if ( scalar @infile_name_prefixes > 1 ) {

        ## picardtools_mergesamfiles
        say {$FILEHANDLE} q{## Merging alignment files};

        Readonly my $JAVA_MEMORY_ALLOCATION => 4;

        # Division by X according to java Heap size
        $core_number =
          floor( $active_parameter_href->{node_ram_memory} / $JAVA_MEMORY_ALLOCATION );

        ## Limit number of cores requested to the maximum number of cores available per node
        $core_number = check_max_core_number(
            {
                core_number_requested => $core_number,
                max_cores_per_node    => $active_parameter_href->{max_cores_per_node},
            }
        );

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                core_number          => $core_number,
                FILEHANDLE           => $FILEHANDLE,
                file_path            => $file_path,
                first_command        => q{java},
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                java_jar =>
                  catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
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
            my @merge_temp_infile_paths =
              map { $_ . $DOT . $contig . $infile_suffix } @temp_infile_path_prefixes;
            my $stderrfile_path =
              $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};

            picardtools_mergesamfiles(
                {
                    create_index       => q{true},
                    FILEHANDLE         => $XARGSFILEHANDLE,
                    infile_paths_ref   => \@merge_temp_infile_paths,
                    outfile_path       => $temp_outfile_path{$contig},
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
            foreach my $temp_infile_path_prefix (@temp_infile_path_prefixes) {

                ## Get parameters
                my $gnu_temp_infile_path =
                  $temp_infile_path_prefix . $DOT . $contig . $infile_suffix;
                my $gnu_temp_outfile_path = $temp_outfile_path{$contig};

                ## Rename
                gnu_mv(
                    {
                        FILEHANDLE   => $XARGSFILEHANDLE,
                        infile_path  => $gnu_temp_infile_path,
                        outfile_path => $gnu_temp_outfile_path,
                    }
                );
                print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

                ## Index
                samtools_index(
                    {
                        bai_format  => 1,
                        FILEHANDLE  => $XARGSFILEHANDLE,
                        infile_path => $gnu_temp_outfile_path,
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
            file_ending => substr( $infile_suffix, 0, 2 ) . $ASTERISK,
            indirectory => $temp_directory,
            infiles_ref => \@temp_infile_name_prefixes,
        }
    );

    # Track the number of created xargs scripts per module
    return $xargs_file_counter;
}
1;
