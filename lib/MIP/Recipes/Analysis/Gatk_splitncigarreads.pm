package MIP::Recipes::Analysis::Gatk_splitncigarreads;

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
    our @EXPORT_OK = qw{ analysis_gatk_splitncigarreads };

}

## Constants
Readonly my $ASTERISK => q{*};
Readonly my $DOT      => q{.};
Readonly my $NEWLINE  => qq{\n};

sub analysis_gatk_splitncigarreads {

## Function : GATK SplitNCigarReads to splits reads into exon segments and hard-clip any sequences overhanging into the intronic regions and reassign mapping qualities from STAR.
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
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
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
    use MIP::Get::File qw{ get_io_files get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Gnu::Coreutils qw{ gnu_cp };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::Program::Alignment::Gatk qw{ gatk_splitncigarreads };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
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
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my %temp_infile_path   = %{ $io{temp}{file_path_href} };
    my $infile_path_prefix = $io{in}{file_path_prefix};
    my $program_mode       = $active_parameter_href->{$program_name};
    my $job_id_chain       = $parameter_href->{$program_name}{chain};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Build outfile paths
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
    my $program_tag   = $parameter_href->{$program_name}{file_tag};
    my @outfile_paths = map {
        catfile( $outsample_directory,
            $infile_name_prefix . $program_tag . $DOT . $_ . $outfile_suffix )
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
    my $outdir_path              = $io{out}{dir_path};
    my %temp_outfile_path        = %{ $io{temp}{file_path_href} };
    my $temp_outfile_path_prefix = $io{temp}{file_path_prefix};
    my $xargs_file_path_prefix;

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

    ### SHELL

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    gnu_cp(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $infile_path_prefix . $ASTERISK,
            outfile_path => $temp_directory,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Division by X according to the java heap
    Readonly my $JAVA_MEMORY_ALLOCATION => 12;
    Readonly my $MEMORY_SPLIT           => 1.25;
    my $node_ram_split =
      $active_parameter_href->{node_ram_memory} / $MEMORY_SPLIT;
    $core_number = floor( $node_ram_split / $JAVA_MEMORY_ALLOCATION );

    ## Limit number of cores requested to the maximum number of cores available per node
    $core_number = check_max_core_number(
        {
            core_number_requested => $core_number,
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    ## GATK SplitNCigarReads
    say {$FILEHANDLE} q{## GATK SplitNCigarReads};

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

        my $infile_path  = $temp_infile_path{$contig};
        my $outfile_path = $temp_outfile_path{$contig};
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};

        gatk_splitncigarreads(
            {
                FILEHANDLE  => $XARGSFILEHANDLE,
                infile_path => $infile_path,
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                memory_allocation => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
                outfile_path      => $outfile_path,
                referencefile_path =>
                  $active_parameter_href->{human_genome_reference},
                stderrfile_path => $stderrfile_path,
                temp_directory  => $temp_directory,
                verbosity       => $active_parameter_href->{gatk_logging_level},
                xargs_mode      => 1,
            }
        );
        print {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## Copies out files from temporary directory.
    say {$FILEHANDLE} q{## Copy file(s) from temporary directory};
    gnu_cp(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $temp_outfile_path_prefix . $ASTERISK,
            outfile_path => $outdir_path,
        }
    );

    close $FILEHANDLE;
    close $XARGSFILEHANDLE;

    if ( $program_mode == 1 ) {

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
