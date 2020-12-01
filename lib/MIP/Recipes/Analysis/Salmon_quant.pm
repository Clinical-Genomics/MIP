package MIP::Recipes::Analysis::Salmon_quant;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw{ uniq };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $NEWLINE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.15;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_salmon_quant };

}

sub analysis_salmon_quant {

## Function : Transcript quantification using salmon quant
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

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
    my $temp_directory;

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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info qw{ get_sample_file_attribute };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources};
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Program::Salmon qw{ salmon_quant };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Sample_info qw{ set_file_path_to_store set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    ## Get the io infiles per chain and id
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
    my @infile_paths = @{ $io{in}{file_paths} };
    my $recipe_mode  = $active_parameter_href->{$recipe_name};
    my $referencefile_dir_path =
        $active_parameter_href->{salmon_quant_reference_genome}
      . $file_info_href->{salmon_quant_reference_genome}[0];
    my %rec_atr = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $job_id_chain    = $rec_atr{chain};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set outfile
    my $recipe_dir =
      catdir( $active_parameter_href->{outdata_dir}, $sample_id, $recipe_name );
    my $file_path = catfile( $recipe_dir, $rec_atr{file_tag} . $rec_atr{outfile_suffix} );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id       => $job_id_chain,
                id             => $sample_id,
                file_info_href => $file_info_href,
                file_paths_ref => [$file_path],
                outdata_dir    => $active_parameter_href->{outdata_dir},
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
                temp_directory => $temp_directory,
            }
        )
    );
    my $outdir_path  = $io{out}{dir_path};
    my $outfile_name = $io{out}{file_names}->[0];
    my $outfile_path = $io{out}{file_paths}->[0];

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe_resource{core_number},
            directory_id          => $sample_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe_resource{memory},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            process_time          => $recipe_resource{time},
            temp_directory        => $temp_directory,
        }
    );

    ### SHELL

    ## Salmon quant
    say {$filehandle} q{## Quantifying transcripts using } . $recipe_name;

    my @infile_prefixes = get_sample_file_attribute(
        {
            file_info_href => $file_info_href,
            file_name      => q{no_direction_infile_prefixes},
            sample_id      => $sample_id,
        }
    );

    my $sequence_run_type = get_sample_file_attribute(
        {
            attribute      => q{sequence_run_type},
            file_info_href => $file_info_href,
            file_name      => $infile_prefixes[0],
            sample_id      => $sample_id,
        }
    );

    ## For paired end
    if ( $sequence_run_type eq q{paired-end} ) {

        ## Grep every other read file and place in new arrays
        # Even array indexes get a 0 remainder and are evalauted as false
        my @read_1_fastq_paths =
          @infile_paths[ grep { !( $_ % 2 ) } 0 .. $#infile_paths ];

        # Odd array indexes get a 1 remainder and are evalauted as true
        my @read_2_fastq_paths = @infile_paths[ grep { $_ % 2 } 0 .. $#infile_paths ];

        salmon_quant(
            {
                filehandle             => $filehandle,
                gc_bias                => 1,
                index_path             => $referencefile_dir_path,
                outdir_path            => $outdir_path,
                read_1_fastq_paths_ref => \@read_1_fastq_paths,
                read_2_fastq_paths_ref => \@read_2_fastq_paths,
            }
        );
        say {$filehandle} $NEWLINE;
    }
    ## For single end
    else {

        salmon_quant(
            {
                filehandle             => $filehandle,
                gc_bias                => 1,
                index_path             => $referencefile_dir_path,
                outdir_path            => $outdir_path,
                read_1_fastq_paths_ref => \@infile_paths,
            }
        );
        say {$filehandle} $NEWLINE;
    }

    ## Close filehandleS
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => $outfile_name,
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        set_file_path_to_store(
            {
                format           => q{meta},
                id               => $sample_id,
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command         => $profile_base_command,
                case_id              => $case_id,
                dependency_method    => q{sample_to_sample},
                job_id_chain         => $job_id_chain,
                job_id_href          => $job_id_href,
                job_reservation_name => $active_parameter_href->{job_reservation_name},
                log                  => $log,
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

1;
