package MIP::Recipes::Analysis::Split_fastq_file;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{fileparse};
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{:all};
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_split_fastq_file };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $DOT        => q{.};
Readonly my $EMPTY_STR  => q{};
Readonly my $NEWLINE    => qq{\n};
Readonly my $PIPE       => q{|};
Readonly my $SPACE      => q{ };
Readonly my $TAB        => qq{\t};
Readonly my $UNDERSCORE => q{_};

sub analysis_split_fastq_file {

## Function : Split input fastq files into batches of reads, versions and compress. Moves original file to subdirectory.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $temp_directory;

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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_module_parameters get_program_attributes };
    use MIP::Gnu::Coreutils qw{ gnu_cp gnu_mkdir gnu_mv gnu_rm gnu_split };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Compression::Pigz qw{ pigz };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constants
    Readonly my $FASTQC_SEQUENCE_LINE_BLOCK => 4;
    Readonly my $SUFFIX_LENGTH              => 4;

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
    my @infile_names              = @{ $io{in}{file_names} };
    my @infile_paths              = @{ $io{in}{file_paths} };
    my $infile_suffix             = $io{in}{file_constant_suffix};
    my @temp_infile_path_prefixes = @{ $io{temp}{file_path_prefixes} };

    my $job_id_chain = get_program_attributes(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
            attribute      => q{chain},
        }
    );
    my $program_mode        = $active_parameter_href->{$program_name};
    my $sequence_read_batch = $active_parameter_href->{split_fastq_file_read_batch};
    my ( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

  INFILE:
    while ( my ( $infile_index, $infile_path ) = each @infile_paths ) {

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        my ($file_name) = setup_script(
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

        say {$FILEHANDLE} q{## } . $program_name;

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
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $infile_path,
                processes   => $core_number,
                stdout      => 1,
            }
        );
        print {$FILEHANDLE} $PIPE . $SPACE;    #Pipe

        gnu_split(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => q{-},
                lines       => ( $sequence_read_batch * $FASTQC_SEQUENCE_LINE_BLOCK ),
                numeric_suffixes => 1,
                prefix           => $temp_infile_path_prefixes[$infile_index]
                  . $UNDERSCORE
                  . q{splitted}
                  . $UNDERSCORE,
                suffix_length => $SUFFIX_LENGTH,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

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
                FILEHANDLE           => $FILEHANDLE,
                infile_suffix        => $DOT . $splitted_suffix,
                temp_directory       => $temp_directory,
            }
        );

        say {$FILEHANDLE} $NEWLINE . q{## Compress file(s) again};
        ## Compress file again
        pigz(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $splitted_flowcell_name_prefix
                  . q{*-SP*}
                  . $DOT
                  . $splitted_suffix,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Copies files from temporary folder to source
        gnu_cp(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $splitted_flowcell_name_prefix . q{*-SP*} . $infile_suffix,
                outfile_path => $indir_path_prefix,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        gnu_mkdir(
            {
                FILEHANDLE => $FILEHANDLE,
                indirectory_path =>
                  catfile( $indir_path_prefix, q{original_fastq_files}, ),
                parents => 1,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Move original file to not be included in subsequent analysis
        say {$FILEHANDLE}
          q{## Move original file to not be included in subsequent analysis};
        gnu_mv(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $infile_path,
                outfile_path => catfile(
                    $indir_path_prefix, q{original_fastq_files},
                    $infile_names[$infile_index]
                ),
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        if ( $program_mode == 1 ) {

            submit_recipe(
                {
                    active_parameter_href   => $active_parameter_href,
                    dependency_method       => q{sample_to_island},
                    family_id               => $family_id,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    log                     => $log,
                    job_id_chain            => $job_id_chain,
                    sample_id               => $sample_id,
                    recipe_file_name        => $file_name,
                }
            );
        }
    }
    close $FILEHANDLE;
    return;
}

sub _list_all_splitted_files {

## Function : List all splitted files
## Returns  :
## Arguments: $fastq_file_info_href => Fastq file info {REF}
##          : $FILEHANDLE           => Filehandle to write to
##          : $infile_suffix        => Infile suffix
##          : $temp_directory       => Temprary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fastq_file_info_href;
    my $FILEHANDLE;
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
        FILEHANDLE    => { store => \$FILEHANDLE },
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
    say {$FILEHANDLE} q{splitted_files=(}
      . catfile( $temp_directory_quoted, q{*_splitted_*} )
      . q{)}, $NEWLINE;

    ## Iterate through array using a counter
    say {$FILEHANDLE}
      q?for ((file_counter=0; file_counter<${#splitted_files[@]}; file_counter++)); do ?;

    ## Rename each element of array to include splitted suffix in flowcell id
    print {$FILEHANDLE} $TAB . q?mv "${splitted_files[$file_counter]}" ?;
    print {$FILEHANDLE} catfile( $temp_directory_quoted, $EMPTY_STR );
    print {$FILEHANDLE} $fastq_file_info_href->{lane} . $UNDERSCORE;
    print {$FILEHANDLE} $fastq_file_info_href->{date} . $UNDERSCORE;
    print {$FILEHANDLE} $fastq_file_info_href->{flowcell} . q?-SP"$file_counter"?;
    print {$FILEHANDLE} $UNDERSCORE . $fastq_file_info_href->{sample_id} . $UNDERSCORE;
    print {$FILEHANDLE} $fastq_file_info_href->{index} . $UNDERSCORE;
    print {$FILEHANDLE} $fastq_file_info_href->{direction} . $infile_suffix;
    say   {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} $TAB . q?echo "${splitted_files[$file_counter]}" ?;
    say {$FILEHANDLE} q{done};
    return;
}
1;
