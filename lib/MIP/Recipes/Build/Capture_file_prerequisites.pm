package MIP::Recipes::Build::Capture_file_prerequisites;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $DOT $NEWLINE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.08;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ build_capture_file_prerequisites };

}

sub build_capture_file_prerequisites {

## Function : Creates the target "interval_list" and  "padded.interval_list" files.
## Returns  :
## Arguments: $active_parameter_href        => Active parameters for this analysis hash {REF}
##          : $case_id                      => Family ID
##          : $filehandle                   => Filehandle to write to
##          : $file_info_href               => File info hash {REF}
##          : $job_id_href                  => Job id hash {REF}
##          : $log                          => Log object
##          : $parameter_build_suffixes_ref => Exome target bed associated file endings
##          : $parameter_href               => Parameter hash {REF}
##          : $profile_base_command         => Submission profile base command
##          : $recipe_name                  => Program name
##          : $sample_info_href             => Info on samples and case hash {REF}
##          : $temp_directory               => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $filehandle;
    my $file_info_href;
    my $job_id_href;
    my $log;
    my $parameter_build_suffixes_ref;
    my $parameter_href;
    my $profile_base_command;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
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
        filehandle     => { store => \$filehandle, },
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

    use MIP::Get::Parameter qw{ get_recipe_resources };
    use MIP::Program::Gnu::Coreutils qw{ gnu_rm gnu_cat gnu_ln };
    use MIP::Language::Shell qw{ check_exist_and_move_file };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Picardtools qw{ picardtools_createsequencedictionary };
    use MIP::Program::Picardtools qw{ picardtools_intervallisttools };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $NR_OF_BASES_PADDING => 100;
    Readonly my $MAX_RANDOM_NUMBER   => 100_00;

    my $recipe_file_path;
    my $submit_switch;

    ## Unpack parameters
    my $interval_list_suffix        = $parameter_build_suffixes_ref->[0];
    my $padded_interval_list_suffix = $parameter_build_suffixes_ref->[1];
    my $recipe_mode                 = $active_parameter_href->{$recipe_name};
    my %recipe_resource             = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => q{mip},
        }
    );
    my $referencefile_path = $active_parameter_href->{human_genome_reference};

    ## No supplied filehandle i.e. create new sbatch script
    if ( not defined $filehandle ) {

        $submit_switch = 1;

        ## Create anonymous filehandle
        $filehandle = IO::Handle->new();

        ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
        ($recipe_file_path) = setup_script(
            {
                active_parameter_href           => $active_parameter_href,
                filehandle                      => $filehandle,
                directory_id                    => $case_id,
                job_id_href                     => $job_id_href,
                log                             => $log,
                recipe_directory                => $recipe_name,
                recipe_name                     => $recipe_name,
                source_environment_commands_ref => $recipe_resource{load_env_ref},
            }
        );
    }

    ## Generate a random integer.
    my $random_integer = int rand $MAX_RANDOM_NUMBER;

  BED_FILE:
    foreach
      my $exome_target_bed_file ( keys %{ $active_parameter_href->{exome_target_bed} } )
    {

        $log->warn( q{Will try to create required }
              . $exome_target_bed_file
              . q{ associated file(s) before executing }
              . $recipe_name );

        ## Add random integer
        my $exome_target_bed_file_random =
          $exome_target_bed_file . $UNDERSCORE . $random_integer;

        say {$filehandle} q{## CreateSequenceDictionary from reference};

        picardtools_createsequencedictionary(
            {
                filehandle => $filehandle,
                java_jar =>
                  catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                memory_allocation    => q{Xmx2g},
                outfile_path         => $exome_target_bed_file_random . $DOT . q{dict},
                referencefile_path   => $referencefile_path,
                temp_directory       => $temp_directory,
            }
        );
        say {$filehandle} $NEWLINE;

        say {$filehandle} q{## Add target file to headers from sequence dictionary};
        gnu_cat(
            {
                filehandle       => $filehandle,
                infile_paths_ref => [
                    $exome_target_bed_file_random . $DOT . q{dict},
                    $exome_target_bed_file
                ],
                stdoutfile_path => $exome_target_bed_file_random . $DOT . q{dict_body},
            }
        );
        say {$filehandle} $NEWLINE;

        ## Remove track and browser header lines and reformat columns in file
        _reformat_capture_file(
            {
                exome_target_bed_file_random => $exome_target_bed_file_random,
                filehandle                   => $filehandle,
            }
        );

        say {$filehandle} q{## Create} . $interval_list_suffix;

        my @infile_paths_ref =
          ( $exome_target_bed_file_random . $DOT . q{dict_body_col_5.interval_list} );
        my $interval_list_outfile_path =
            $exome_target_bed_file_random
          . $DOT
          . q{dict_body_col_5}
          . $UNDERSCORE
          . $interval_list_suffix;
        picardtools_intervallisttools(
            {
                filehandle       => $filehandle,
                infile_paths_ref => \@infile_paths_ref,
                java_jar =>
                  catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                memory_allocation    => q{Xmx2g},
                outfile_path         => $interval_list_outfile_path,
                referencefile_path   => $referencefile_path,
                temp_directory       => $temp_directory,
            }
        );
        say {$filehandle} $NEWLINE;

        my $intended_file_path = $exome_target_bed_file . $interval_list_suffix;
        my $temporary_file_path =
            $exome_target_bed_file_random
          . $DOT
          . q{dict_body_col_5}
          . $UNDERSCORE
          . $interval_list_suffix;

        ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
        check_exist_and_move_file(
            {
                filehandle          => $filehandle,
                intended_file_path  => $intended_file_path,
                temporary_file_path => $temporary_file_path,
            }
        );

        say {$filehandle} q{#Create} . $padded_interval_list_suffix;

        my $padded_interval_list_outfile_path =
            $exome_target_bed_file_random
          . $DOT
          . q{dict_body_col_5}
          . $padded_interval_list_suffix;
        picardtools_intervallisttools(
            {
                filehandle       => $filehandle,
                infile_paths_ref => \@infile_paths_ref,
                java_jar =>
                  catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                memory_allocation    => q{Xmx2g},
                outfile_path         => $padded_interval_list_outfile_path,
                padding              => $NR_OF_BASES_PADDING,
                referencefile_path   => $referencefile_path,
                temp_directory       => $active_parameter_href->{temp_directory},
            }
        );
        say {$filehandle} $NEWLINE;

        $intended_file_path = $exome_target_bed_file . $padded_interval_list_suffix;
        $temporary_file_path =
            $exome_target_bed_file_random
          . $DOT
          . q{dict_body_col_5}
          . $padded_interval_list_suffix;

        ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
        check_exist_and_move_file(
            {
                filehandle          => $filehandle,
                intended_file_path  => $intended_file_path,
                temporary_file_path => $temporary_file_path,
            }
        );

        ## Remove temporary files
        say {$filehandle} q{#Remove temporary files};

        my @temp_files = (
            $exome_target_bed_file_random . $DOT . q{dict_body_col_5.interval_list},
            $exome_target_bed_file_random . $DOT . q{dict_body},
            $exome_target_bed_file_random . $DOT . q{dict},
        );
      FILE_TO_REMOVE:
        foreach my $file (@temp_files) {

            gnu_rm(
                {
                    filehandle  => $filehandle,
                    force       => 1,
                    infile_path => $file,
                }
            );
            say {$filehandle} $NEWLINE;
        }
    }

    ## Only create once
    $parameter_href->{exome_target_bed}{build_file} = 0;

    ## Unless filehandle was supplied close filehandle and submit
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

sub _reformat_capture_file {

## Function : Remove track and browser header lines and reformat columns in file
## Returns  :
## Arguments: $exome_target_bed_file_random => Exome target bed file
##          : $filehandle                   => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $exome_target_bed_file_random;
    my $filehandle;

    my $tmpl = {
        exome_target_bed_file_random => {
            defined     => 1,
            required    => 1,
            store       => \$exome_target_bed_file_random,
            strict_type => 1,
        },
        filehandle => { store => \$filehandle, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Remove unnecessary info and reformat
    say {$filehandle}
      q{#Remove target annotations, 'track', 'browse' and keep only 5 columns};

    ## Execute perl
    print {$filehandle} q?perl -nae '?;

    ## If header line
    print {$filehandle} q?if ($_=~/@/) { ?;

    ## Write to stdout
    print {$filehandle} q?print $_;} ?;

    ## If track - do nothing
    print {$filehandle} q?elsif ($_=~/^track/) {} ?;

    ## If browser - do nothing
    print {$filehandle} q?elsif ($_=~/^browser/) {} ?;

    ## Else print reformated line to stdout
    print {$filehandle}
q?else {print @F[0], "\t", (@F[1] + 1), "\t", @F[2], "\t", "+", "\t", "-", "\n";}' ?;

    ## Infile
    print {$filehandle} $exome_target_bed_file_random . $DOT . q{dict_body} . $SPACE;

    ## Write to
    print {$filehandle} q{>} . $SPACE;
    say   {$filehandle} $exome_target_bed_file_random . $DOT
      . q{dict_body_col_5.interval_list}, $NEWLINE;

    return;
}

1;
