package MIP::Recipes::Build::Human_genome_prerequisites;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $DOT $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.07;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ build_human_genome_prerequisites };

}

sub build_human_genome_prerequisites {

## Function : Creates the human genome prerequisites using active_parameters{human_genome_reference} as reference.
## Returns  :
## Arguments: $active_parameter_href        => Active parameters for this analysis hash {REF}
##          : $case_id                      => Family ID
##          : $filehandle                   => Filehandle to write to. A new sbatch script will be generated if $filehandle is lacking, else write to exising $filehandle {Optional}
##          : $file_info_href               => File info hash {REF}
##          : $human_genome_reference       => Human genome reference
##          : $job_id_href                  => Job id hash {REF}
##          : $log                          => Log object
##          : $parameter_build_suffixes_ref => The human genome reference associated file endings {REF}
##          : $parameter_href               => Parameter hash {REF}
##          : $profile_base_command         => Submission profile base command
##          : $recipe_name                  => Program under evaluation
##          : $random_integer               => The random integer to create temporary file name
##          : $reference_dir                => MIP reference directory
##          : $sample_info_href             => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $filehandle;
    my $file_info_href;
    my $job_id_href;
    my $log;
    my $parameter_build_suffixes_ref;
    my $parameter_href;
    my $recipe_name;
    my $random_integer;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $human_genome_reference;
    my $profile_base_command;
    my $reference_dir;

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
        filehandle     => { store => \$filehandle, },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        human_genome_reference => {
            default     => $arg_href->{active_parameter_href}{human_genome_reference},
            store       => \$human_genome_reference,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        log                          => { store => \$log, },
        parameter_build_suffixes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$parameter_build_suffixes_ref,
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
        random_integer => { store => \$random_integer, strict_type => 1, },
        reference_dir  => {
            default     => $arg_href->{active_parameter_href}{reference_dir},
            store       => \$reference_dir,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::Parameter qw{ get_recipe_resources };
    use MIP::Program::Gnu::Coreutils qw{ gnu_rm gnu_ln };
    use MIP::Language::Java qw{ java_core };
    use MIP::Language::Shell qw{ check_exist_and_move_file };
    use MIP::Program::Gzip qw{ gzip };
    use MIP::Program::Samtools qw{ samtools_faidx };
    use MIP::Program::Picardtools qw{ picardtools_createsequencedictionary };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Recipes::Build::Capture_file_prerequisites
      qw{ build_capture_file_prerequisites };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $MAX_RANDOM_NUMBER => 100_00;

    my $recipe_file_path;
    my $submit_switch;

    ## Unpack parameters
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => q{mip},
        }
    );

    ## No supplied filehandle i.e. create new sbatch script
    if ( not defined $filehandle ) {

        $submit_switch = 1;

        ## Create anonymous filehandle
        $filehandle = IO::Handle->new();

        ## Generate a random integer between 0-10,000.
        $random_integer = int rand $MAX_RANDOM_NUMBER;

        ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
        ($recipe_file_path) = setup_script(
            {
                active_parameter_href           => $active_parameter_href,
                job_id_href                     => $job_id_href,
                filehandle                      => $filehandle,
                directory_id                    => $case_id,
                log                             => $log,
                recipe_name                     => $recipe_name,
                recipe_directory                => $recipe_name,
                source_environment_commands_ref => $recipe_resource{load_env_ref},
            }
        );
    }

    ## Check for compressed files
    if ( $file_info_href->{human_genome_compressed} ) {

        $log->warn( q{Will try to decompress }
              . $human_genome_reference
              . q{ before executing }
              . $recipe_name );

        ## Perl wrapper for writing gzip recipe to $filehandle
        gzip(
            {
                decompress       => 1,
                infile_paths_ref => [$human_genome_reference],
                filehandle       => $filehandle,
            }
        );
        say {$filehandle} $NEWLINE;

        ## Replace the .fasta.gz ending with .fasta since this will execute before the analysis
        ## Hence, changing the original file name ending from ".fastq" to ".fastq.gz".
        $human_genome_reference =~ s/.fasta.gz/.fasta/xsmg;

        $log->info( q{Set human_genome_reference to: } . $human_genome_reference );
        $file_info_href->{human_genome_compressed} = 0;
    }

    if ( exists $parameter_href->{exome_target_bed}{build_file}
        and $parameter_href->{exome_target_bed}{build_file} == 1 )
    {

        build_capture_file_prerequisites(
            {
                active_parameter_href        => $active_parameter_href,
                filehandle                   => $filehandle,
                file_info_href               => $file_info_href,
                job_id_href                  => $job_id_href,
                log                          => $log,
                parameter_build_suffixes_ref => \@{ $file_info_href->{exome_target_bed} },
                parameter_href               => $parameter_href,
                recipe_name                  => $recipe_name,
                sample_info_href             => $sample_info_href,
            }
        );
    }
    if ( $parameter_href->{human_genome_reference_file_endings}{build_file} == 1 ) {

      FILE_ENDING:
        foreach my $file_ending ( @{$parameter_build_suffixes_ref} ) {

            if ( $file_ending eq $DOT . q{dict} ) {

                $log->warn( q{Will try to create }
                      . $file_ending
                      . q{ file for }
                      . $human_genome_reference
                      . q{ before executing }
                      . $recipe_name );

                my $filename_prefix = catfile( $reference_dir,
                    $file_info_href->{human_genome_reference_name_prefix} );

                say {$filehandle} q{#CreateSequenceDictionary from reference};

                picardtools_createsequencedictionary(
                    {
                        filehandle => $filehandle,
                        java_jar   => catfile(
                            $active_parameter_href->{picardtools_path},
                            q{picard.jar}
                        ),
                        java_use_large_pages =>
                          $active_parameter_href->{java_use_large_pages},
                        memory_allocation => q{Xmx2g},
                        outfile_path      => $filename_prefix
                          . $UNDERSCORE
                          . $random_integer
                          . $file_ending,
                        referencefile_path => $human_genome_reference,
                        temp_directory     => $active_parameter_href->{temp_directory},
                    }
                );
                say {$filehandle} $NEWLINE;

                my $intended_file_path = $filename_prefix . $file_ending;
                my $temporary_file_path =
                  $filename_prefix . $UNDERSCORE . $random_integer . $file_ending;

                ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
                check_exist_and_move_file(
                    {
                        filehandle          => $filehandle,
                        intended_file_path  => $intended_file_path,
                        temporary_file_path => $temporary_file_path,
                    }
                );
            }
            if ( $file_ending eq $DOT . q{fai} ) {

                $log->warn( q{Will try to create }
                      . $file_ending
                      . q{ file for }
                      . $human_genome_reference
                      . q{ before executing }
                      . $recipe_name );

                my $human_genome_reference_temp_file =
                  $human_genome_reference . $UNDERSCORE . $random_integer;

                say {$filehandle} q{## Fai file from reference};
                gnu_ln(
                    {
                        filehandle  => $filehandle,
                        force       => 1,
                        link_path   => $human_genome_reference_temp_file,
                        symbolic    => 1,
                        target_path => $human_genome_reference,
                    }
                );
                say {$filehandle} $NEWLINE;

                samtools_faidx(
                    {
                        filehandle  => $filehandle,
                        infile_path => $human_genome_reference_temp_file,
                    }
                );
                say {$filehandle} $NEWLINE;

                my $intended_file_path = $human_genome_reference . $file_ending;
                my $temporary_file_path =
                  $human_genome_reference_temp_file . $file_ending;

                ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
                check_exist_and_move_file(
                    {
                        filehandle          => $filehandle,
                        intended_file_path  => $intended_file_path,
                        temporary_file_path => $temporary_file_path,
                    }
                );

                ## Remove soft link
                gnu_rm(
                    {
                        filehandle  => $filehandle,
                        force       => 1,
                        infile_path => $human_genome_reference_temp_file,
                    }
                );
                say {$filehandle} $NEWLINE;
            }
        }

        ## Only create once
        $parameter_href->{human_genome_reference_file_endings}{build_file} = 0;
    }

    ## Unless filehandle was supplied close it and submit
    if ($submit_switch) {

        close $filehandle;

        if ( $recipe_mode == 1 ) {

            submit_recipe(
                {
                    base_command       => $profile_base_command,
                    dependency_method  => q{island_to_samples},
                    case_id            => $case_id,
                    job_id_href        => $job_id_href,
                    log                => $log,
                    job_id_chain       => q{MAIN},
                    recipe_file_path   => $recipe_file_path,
                    sample_ids_ref     => \@{ $active_parameter_href->{sample_ids} },
                    submission_profile => $active_parameter_href->{submission_profile},
                }
            );
        }
    }
    return 1;
}

1;
