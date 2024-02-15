package MIP::Recipes::Analysis::Fastp;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname fileparse };
use File::Spec::Functions qw{ catdir catfile devnull splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_fastp };

}

sub analysis_fastp {

    ## Function : Trim incoming fastq files using fastp.
    ## Returns  :
    ## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
    ##          : $case_id                 => Family id
    ##          : $file_info_href          => File info hash {REF}
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info qw{ get_io_files get_sample_file_attribute parse_io_outfiles set_io_files };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Fastp qw{ fastp };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{
      set_recipe_metafile_in_sample_info
      set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
        }
    );
    my @infile_paths = @{ $io{in}{file_paths} };

    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my %recipe             = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    my %file_info_sample = get_sample_file_attribute(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $recipe{job_id_chain},
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => $file_info_sample{no_direction_infile_prefixes},
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );

    my $outfile_suffix        = $io{out}{file_constant_suffix};
    my @outfile_name_prefixes = @{ $io{out}{file_name_prefixes} };
    my @outfile_path_prefixes = @{ $io{out}{file_path_prefixes} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    # Too avoid adjusting infile_index in submitting to jobs
    my $paired_end_tracker = 0;

    my @outfile_paths;

    ## Perform per single-end or read pair
  INFILE_PREFIX:
    while ( my ( $infile_index, $infile_prefix ) =
        each @{ $file_info_sample{no_direction_infile_prefixes} } )
    {

        ## Assign file features
        my $outfile_name_prefix = $outfile_name_prefixes[$infile_index];

        # Collect interleaved status for fastq file
        my $sequence_run_type = get_sample_file_attribute(
            {
                attribute      => q{sequence_run_type},
                file_info_href => $file_info_href,
                file_name      => $infile_prefix,
                sample_id      => $sample_id,
            }
        );
        my $is_interleaved_fastq = $sequence_run_type eq q{interleaved} ? 1 : 0;

        ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
        my ( $recipe_file_path, $recipe_info_path ) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                core_number           => $recipe{core_number},
                directory_id          => $sample_id,
                filehandle            => $filehandle,
                job_id_href           => $job_id_href,
                memory_allocation     => $recipe{memory},
                recipe_directory      => $recipe_name,
                recipe_name           => $recipe_name,
                process_time          => $recipe{time},
                temp_directory        => $temp_directory,
            }
        );

        ### SHELL:

        ### Fastp
        say {$filehandle} q{## Trimming reads with } . $recipe_name;

        ### Get parameters

        ## Infile(s)
        my $fastq_file_path = $infile_paths[$paired_end_tracker];
        my $second_fastq_file_path;

        ## Outfile(s)
        my $first_outfile_path = $outfile_path_prefixes[$infile_index] . $outfile_suffix;
        my $second_outfile_path;
        my $report_html = $outfile_path_prefixes[$infile_index] . q{.html};
        my $report_json = $outfile_path_prefixes[$infile_index] . q{.json};

        # If second read direction is present
        if ( $sequence_run_type eq q{paired-end} ) {

            ## Set paired end output filenames
            $first_outfile_path  = $outfile_path_prefixes[$infile_index] . q{_1} . $outfile_suffix;
            $second_outfile_path = $outfile_path_prefixes[$infile_index] . q{_2} . $outfile_suffix;

            # Increment to collect correct read 2
            $paired_end_tracker     = $paired_end_tracker + 1;
            $second_fastq_file_path = $infile_paths[$paired_end_tracker];

            push @outfile_paths, ( $first_outfile_path, $second_outfile_path );
        }
        else {

            push @outfile_paths, $first_outfile_path;
        }

        fastp(
            {
                detect_pe_adapter     => $active_parameter_href->{fastp_detect_pe_adapter},
                filehandle            => $filehandle,
                first_infile_path     => $fastq_file_path,
                first_outfile_path    => $first_outfile_path,
                interleaved_in        => $is_interleaved_fastq,
                length_required       => $active_parameter_href->{fastp_length_required},
                low_complexity_filter => $active_parameter_href->{fastp_low_complexity_filter},
                overrepresentation_analysis =>
                  $active_parameter_href->{fastp_overrepresentation_analysis},
                report_html         => $report_html,
                report_json         => $report_json,
                second_infile_path  => $second_fastq_file_path,
                second_outfile_path => $second_outfile_path,
                threads             => $recipe{core_number},
                trim_poly_g         => $active_parameter_href->{fastp_trim_poly_g},
            }
        );
        say {$filehandle} $NEWLINE;

        ## Increment paired end tracker
        $paired_end_tracker++;

        close $filehandle;

        if ( $recipe{mode} == 1 ) {

            set_recipe_outfile_in_sample_info(
                {
                    infile           => $outfile_name_prefix,
                    path             => $first_outfile_path,
                    recipe_name      => $recipe_name,
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );

            submit_recipe(
                {
                    base_command         => $profile_base_command,
                    case_id              => $case_id,
                    dependency_method    => q{sample_to_sample_parallel},
                    job_id_chain         => $recipe{job_id_chain},
                    job_id_href          => $job_id_href,
                    job_reservation_name => $active_parameter_href->{job_reservation_name},
                    log                  => $log,
                    max_parallel_processes_count_href =>
                      $file_info_href->{max_parallel_processes_count},
                    recipe_file_path     => $recipe_file_path,
                    recipe_files_tracker => $infile_index,
                    sample_id            => $sample_id,
                    submission_profile   => $active_parameter_href->{submission_profile},
                }
            );
        }
    }

    ## Set input files for next module
    set_io_files(
        {
            chain_id       => $recipe{job_id_chain},
            id             => $sample_id,
            file_info_href => $file_info_href,
            file_paths_ref => \@outfile_paths,
            recipe_name    => $recipe_name,
            stream         => q{out},
        }
    );

    return 1;
}

1;
