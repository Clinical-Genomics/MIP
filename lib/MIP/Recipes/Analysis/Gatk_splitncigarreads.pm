package MIP::Recipes::Analysis::Gatk_splitncigarreads;

use 5.026;
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

## MIPs lib/
use MIP::Constants qw{ %ANALYSIS $ASTERISK $DOT $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.16;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_gatk_splitncigarreads };

}

## Constants
Readonly my $JAVA_GUEST_OS_MEMORY => $ANALYSIS{JAVA_GUEST_OS_MEMORY};

sub analysis_gatk_splitncigarreads {

## Function : GATK SplitNCigarReads to splits reads into exon segments and hard-clip any sequences overhanging into the intronic regions and reassign mapping qualities from STAR.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
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
    my $file_path;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;
    my $sample_id;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
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

    use MIP::Cluster qw{ get_parallel_processes };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Gnu::Coreutils qw{ gnu_cp };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Gatk qw{ gatk_splitncigarreads };
    use MIP::Sample_info
      qw{ set_recipe_outfile_in_sample_info set_processing_metafile_in_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my %io = get_io_files(
        {
            file_info_href => $file_info_href,
            id             => $sample_id,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my %infile_path        = %{ $io{in}{file_path_href} };
    my $recipe_mode        = $active_parameter_href->{$recipe_name};
    my $job_id_chain       = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number = $recipe_resource{core_number};

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $sample_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => $file_info_href->{bam_contigs_size_ordered},
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
                temp_directory   => $temp_directory,
            }
        )
    );
    my $outfile_name_prefix = $io{out}{file_name_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my @outfile_paths       = @{ $io{out}{file_paths} };
    my %outfile_path        = %{ $io{out}{file_path_href} };
    my $xargs_file_path_prefix;

    ## Filehandles
    # Create anonymous filehandle
    my $xargsfilehandle = IO::Handle->new();
    my $filehandle      = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $sample_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            log                             => $log,
            memory_allocation               => $recipe_resource{memory},
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL

    Readonly my $JAVA_MEMORY_ALLOCATION => 12;
    my $process_memory_allocation = $JAVA_MEMORY_ALLOCATION + $JAVA_GUEST_OS_MEMORY;

    # Constrain parallelization to match available memory
    my $parallel_processes = get_parallel_processes(
        {
            process_memory_allocation => $process_memory_allocation,
            recipe_memory_allocation  => $recipe_resource{memory},
            core_number               => $core_number,
        }
    );

    ## GATK SplitNCigarReads
    say {$filehandle} q{## GATK SplitNCigarReads};

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
    while ( my ( $infile_index, $contig ) =
        each @{ $file_info_href->{bam_contigs_size_ordered} } )
    {

        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};

        gatk_splitncigarreads(
            {
                filehandle           => $xargsfilehandle,
                infile_path          => $infile_path{$contig},
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                memory_allocation    => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
                outfile_path         => $outfile_path{$contig},
                referencefile_path   => $active_parameter_href->{human_genome_reference},
                stderrfile_path      => $stderrfile_path,
                temp_directory       => $temp_directory,
                verbosity            => $active_parameter_href->{gatk_logging_level},
                xargs_mode           => 1,
            }
        );
        print {$xargsfilehandle} $NEWLINE;
    }

    close $filehandle;
    close $xargsfilehandle;

    if ( $recipe_mode == 1 ) {

        my $first_outfile_path = $outfile_paths[0];

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => $outfile_name_prefix,
                path             => $first_outfile_path,
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        my $most_complete_format_key =
          q{most_complete} . $UNDERSCORE . substr $outfile_suffix, 1;
        set_processing_metafile_in_sample_info(
            {
                metafile_tag     => $most_complete_format_key,
                path             => $first_outfile_path,
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
                job_reservation_name    => $active_parameter_href->{job_reservation_name},
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
