package MIP::Recipes::Analysis::Split_fastq_file;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{fileparse};
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{:all};
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $ASTERISK $DOT $EMPTY_STR $LOG_NAME $NEWLINE $PIPE $SPACE $TAB $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_split_fastq_file };

}

sub analysis_split_fastq_file {

## Function : Split input fastq files into batches of reads, versions and compress. Moves original file to subdirectory.
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

    use MIP::File_info qw{ get_io_files };
    use MIP::Program::Gnu::Coreutils qw{ gnu_cp gnu_mkdir gnu_mv gnu_rm gnu_split };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Pigz qw{ pigz };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constants
    Readonly my $FASTQC_SEQUENCE_LINE_BLOCK => 4;
    Readonly my $SUFFIX_LENGTH              => 4;

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
            temp_directory => $temp_directory,
        }
    );
    my $indir_path_prefix         = $io{in}{dir_path_prefix};
    my @infile_names              = @{ $io{in}{file_names} };
    my @infile_paths              = @{ $io{in}{file_paths} };
    my $infile_suffix             = $io{in}{file_constant_suffix};
    my @temp_infile_path_prefixes = @{ $io{temp}{file_path_prefixes} };

    my $sequence_read_batch = $active_parameter_href->{split_fastq_file_read_batch};
    my %recipe              = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

  INFILE:
    while ( my ( $infile_index, $infile_path ) = each @infile_paths ) {

        ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
        my ($recipe_file_path) = setup_script(
            {
                active_parameter_href           => $active_parameter_href,
                core_number                     => $recipe{core_number},
                directory_id                    => $sample_id,
                filehandle                      => $filehandle,
                job_id_href                     => $job_id_href,
                memory_allocation               => $recipe{memory},
                process_time                    => $recipe{time},
                recipe_directory                => $recipe_name,
                recipe_name                     => $recipe_name,
                source_environment_commands_ref => $recipe{load_env_ref},
                temp_directory                  => $temp_directory,
            }
        );

        say {$filehandle} q{## } . $recipe_name;

        my %fastq_file_info;

        ## Detect fastq file info for later rebuild of filename
        if (
            $infile_path =~ qr{
	      (\d+)_ # Lane
	      (\d+)_ # Date
	      ([^_]+)_ # Flowcell
	      ([^_]+)_ # Sample id
	      ([^_]+)_ # Index
	      (\d) # direction
	      $infile_suffix
	    }sxm
          )
        {

            %fastq_file_info = (
                lane      => $1,
                date      => $2,
                flowcell  => $3,
                sample_id => $4,
                index     => $5,
                direction => $6,
            );
        }

        ### SHELL:

        ## Decompress file and split
        pigz(
            {
                decompress  => 1,
                filehandle  => $filehandle,
                infile_path => $infile_path,
                processes   => $recipe{core_number},
                stdout      => 1,
            }
        );
        print {$filehandle} $PIPE . $SPACE;    #Pipe

        gnu_split(
            {
                filehandle       => $filehandle,
                infile_path      => q{-},
                lines            => ( $sequence_read_batch * $FASTQC_SEQUENCE_LINE_BLOCK ),
                numeric_suffixes => 1,
                prefix           => $temp_infile_path_prefixes[$infile_index]
                  . $UNDERSCORE
                  . q{splitted}
                  . $UNDERSCORE,
                suffix_length => $SUFFIX_LENGTH,
            }
        );
        say {$filehandle} $NEWLINE;

        my $splitted_suffix               = q{fastq};
        my $splitted_flowcell_name_prefix = catfile( $indir_path_prefix,
                $fastq_file_info{lane}
              . $UNDERSCORE
              . $fastq_file_info{date}
              . $UNDERSCORE
              . $fastq_file_info{flowcell} );

        _list_all_splitted_files(
            {
                fastq_file_info_href => \%fastq_file_info,
                filehandle           => $filehandle,
                infile_suffix        => $DOT . $splitted_suffix,
                temp_directory       => $temp_directory,
            }
        );

        say {$filehandle} $NEWLINE . q{## Compress file(s) again};
        ## Compress file again
        pigz(
            {
                filehandle  => $filehandle,
                infile_path => $splitted_flowcell_name_prefix . q{*-SP*} . $DOT . $splitted_suffix,
            }
        );
        say {$filehandle} $NEWLINE;

        ## Copies files from temporary folder to source
        gnu_cp(
            {
                filehandle   => $filehandle,
                infile_path  => $splitted_flowcell_name_prefix . q{*-SP*} . $infile_suffix,
                outfile_path => $indir_path_prefix,
            }
        );
        say {$filehandle} $NEWLINE;

        gnu_mkdir(
            {
                filehandle       => $filehandle,
                indirectory_path => catfile( $indir_path_prefix, q{original_fastq_files}, ),
                parents          => 1,
            }
        );
        say {$filehandle} $NEWLINE;

        ## Move original file to not be included in subsequent analysis
        say {$filehandle} q{## Move original file to not be included in subsequent analysis};
        gnu_mv(
            {
                filehandle   => $filehandle,
                infile_path  => $infile_path,
                outfile_path => catfile(
                    $indir_path_prefix, q{original_fastq_files}, $infile_names[$infile_index]
                ),
            }
        );
        say {$filehandle} $NEWLINE;

        if ( $recipe{mode} == 1 ) {

            submit_recipe(
                {
                    base_command         => $profile_base_command,
                    case_id              => $case_id,
                    dependency_method    => q{sample_to_island},
                    job_id_chain         => $recipe{job_id_chain},
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
    }
    close $filehandle;
    return 1;
}

sub _list_all_splitted_files {

## Function : List all splitted files
## Returns  :
## Arguments: $fastq_file_info_href => Fastq file info {REF}
##          : $filehandle           => Filehandle to write to
##          : $infile_suffix        => Infile suffix
##          : $temp_directory       => Temprary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fastq_file_info_href;
    my $filehandle;
    my $infile_suffix;
    my $temp_directory;

    my $tmpl = {
        fastq_file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$fastq_file_info_href
        },
        filehandle    => { store => \$filehandle },
        infile_suffix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_suffix
        },
        temp_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$temp_directory
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Shell qw{quote_bash_variable};

    ## Double quote incoming variables in string
    my $temp_directory_quoted =
      quote_bash_variable( { string_with_variable_to_quote => $temp_directory, } );

    ## Find all splitted files
    say {$filehandle} q{splitted_files=(}
      . catfile( $temp_directory_quoted, q{*_splitted_*} )
      . q{)}, $NEWLINE;

    ## Iterate through array using a counter
    say {$filehandle}
      q?for ((file_counter=0; file_counter<${#splitted_files[@]}; file_counter++)); do ?;

    ## Rename each element of array to include splitted suffix in flowcell id
    print {$filehandle} $TAB . q?mv "${splitted_files[$file_counter]}" ?;
    print {$filehandle} catfile( $temp_directory_quoted, $EMPTY_STR );
    print {$filehandle} $fastq_file_info_href->{lane} . $UNDERSCORE;
    print {$filehandle} $fastq_file_info_href->{date} . $UNDERSCORE;
    print {$filehandle} $fastq_file_info_href->{flowcell} . q?-SP"$file_counter"?;
    print {$filehandle} $UNDERSCORE . $fastq_file_info_href->{sample_id} . $UNDERSCORE;
    print {$filehandle} $fastq_file_info_href->{index} . $UNDERSCORE;
    print {$filehandle} $fastq_file_info_href->{direction} . $infile_suffix;
    say   {$filehandle} $NEWLINE;

    say {$filehandle} $TAB . q?echo "${splitted_files[$file_counter]}" ?;
    say {$filehandle} q{done};
    return;
}
1;
