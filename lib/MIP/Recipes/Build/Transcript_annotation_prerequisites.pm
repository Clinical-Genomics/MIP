package MIP::Recipes::Build::Transcript_annotation_prerequisites;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $DASH $DOT $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ build_transcript_annotation_prerequisites };

}

sub build_transcript_annotation_prerequisites {

## Function : Creates the transcript annotatin refFlat file.
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

    use MIP::Language::Shell qw{ check_exist_and_move_file };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Ucsc qw{ ucsc_gtf_to_genepred };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $MAX_RANDOM_NUMBER => 100_00;

    my $recipe_file_path;
    my $submit_switch;

    ## Unpack parameters

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Generate a random integer.
    my $random_integer              = int rand $MAX_RANDOM_NUMBER;
    my $annotation_file_path        = $active_parameter_href->{transcript_annotation};
    my $annotation_file_path_random = $annotation_file_path . $UNDERSCORE . $random_integer;

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
                recipe_directory                => $recipe_name,
                recipe_name                     => $recipe_name,
                source_environment_commands_ref => $recipe{load_env_ref},
            }
        );
    }

    $log->warn( q{Will try to create required }
          . $annotation_file_path
          . q{ associated file(s) before executing }
          . $recipe_name );

    my %build_transcript_annotation = (
        q{.bed} => {
            extra_arg_href => {},
            method         => \&_build_bed,
        },
        q{.refflat} => {
            extra_arg_href => {
                annotation_file_path_random => $annotation_file_path_random,
            },
            method => \&_build_refflat,
        },
        q{.rrna.interval_list} => {
            extra_arg_href => {
                active_parameter_href       => $active_parameter_href,
                annotation_file_path_random => $annotation_file_path_random,
            },
            method => \&_build_rrna_interval_list,
        },
    );

  ANNOTATION_SUFFIX:
    foreach my $annotation_suffix ( @{$parameter_build_suffixes_ref} ) {

        my $intended_file_path = $annotation_file_path . $annotation_suffix;
        my $temp_file_path     = $annotation_file_path_random . $annotation_suffix;

        ## Build annotaion
        $build_transcript_annotation{$annotation_suffix}{method}->(
            {
                %{ $build_transcript_annotation{$annotation_suffix}{extra_arg_href} },
                annotation_file_path => $annotation_file_path,
                filehandle           => $filehandle,
                temp_file_path       => $temp_file_path,
            }
        );

        ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
        check_exist_and_move_file(
            {
                filehandle          => $filehandle,
                intended_file_path  => $intended_file_path,
                temporary_file_path => $temp_file_path,
            }
        );
    }

    ## Only create once
    $parameter_href->{transcript_annotation_file_endings}{build_file} = 0;

    ## Unless filehandle was supplied close filehandle and submit
    if ($submit_switch) {

        close $filehandle;

        if ( $recipe{mode} == 1 ) {

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

sub _build_bed {

## Function : Creates the transcript annotation bed file
## Returns  :
## Arguments: $annotation_file_path => Annotation file path
##          : $filehandle           => Filehandle to write to
##          : $temp_file_path       => Temp file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $annotation_file_path;
    my $filehandle;
    my $temp_file_path;

    my $tmpl = {
        annotation_file_path => {
            required    => 1,
            store       => \$annotation_file_path,
            strict_type => 1,
        },
        filehandle     => { store => \$filehandle, },
        temp_file_path => {
            required    => 1,
            store       => \$temp_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Gtf2bed qw{ gtf2bed };

    say {$filehandle} q{## Converting to bed format};
    gtf2bed(
        {
            filehandle      => $filehandle,
            infile_path     => $annotation_file_path,
            stdoutfile_path => $temp_file_path,
        }
    );
    say {$filehandle} $NEWLINE;

    return;
}

sub _build_refflat {

## Function : Creates the transcript annotation refflat file.
## Returns  :
## Arguments: $annotation_file_path        => Annotation file path
##          : $annotation_file_path_random => Annotation suffix
##          : $filehandle                  => Filehandle to write to
##          : $temp_file_path              => Temp file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $annotation_file_path;
    my $annotation_file_path_random;
    my $filehandle;
    my $temp_file_path;

    my $tmpl = {
        annotation_file_path => {
            required    => 1,
            store       => \$annotation_file_path,
            strict_type => 1,
        },
        annotation_file_path_random => {
            required    => 1,
            store       => \$annotation_file_path_random,
            strict_type => 1,
        },
        filehandle     => { store => \$filehandle, },
        temp_file_path => {
            required    => 1,
            store       => \$temp_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Perl qw{ perl_nae_oneliners };
    use MIP::Program::Gnu::Coreutils qw{ gnu_rm };
    use MIP::Program::Ucsc qw{ ucsc_gtf_to_genepred };

    ## Set file names
    my $temp_genepred_file_path = $annotation_file_path_random . $DOT . q{genePred};

    say {$filehandle} q{## Convert gtf to extended genePred };
    ucsc_gtf_to_genepred(
        {
            extended_genepred  => 1,
            filehandle         => $filehandle,
            gene_name_as_name2 => 1,
            infile_path        => $annotation_file_path,
            outfile_path       => $temp_genepred_file_path,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Convert genePred to refFlat};

    perl_nae_oneliners(
        {
            filehandle      => $filehandle,
            oneliner_name   => q{genepred_to_refflat},
            stdinfile_path  => $temp_genepred_file_path,
            stdoutfile_path => $temp_file_path,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Remove temporary files};
    gnu_rm(
        {
            filehandle  => $filehandle,
            force       => 1,
            infile_path => $temp_genepred_file_path,
        }
    );
    print {$filehandle} $NEWLINE;

    return;
}

sub _build_rrna_interval_list {

## Function : Creates the transcript annotation ribomal RNA interval_list
## Returns  :
## Arguments: $active_parameter_href       => Active parameter hash {REF}
##          : $annotation_file_path        => Annotation file path
##          : $annotation_file_path_random => Annotation suffix
##          : $filehandle                  => Filehandle to write to
##          : $temp_file_path              => Temp file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $annotation_file_path;
    my $annotation_file_path_random;
    my $filehandle;
    my $temp_file_path;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        annotation_file_path => {
            required    => 1,
            store       => \$annotation_file_path,
            strict_type => 1,
        },
        annotation_file_path_random => {
            required    => 1,
            store       => \$annotation_file_path_random,
            strict_type => 1,
        },
        filehandle     => { store => \$filehandle, },
        temp_file_path => {
            required    => 1,
            store       => \$temp_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Perl qw{ perl_nae_oneliners };
    use MIP::Program::Gtf2bed qw{ gtf2bed };
    use MIP::Program::Picardtools qw{
      picardtools_bedtointervallist
      picardtools_createsequencedictionary
    };
    use MIP::Program::Gnu::Coreutils qw{ gnu_rm };

    ## Set file names
    my $temp_rrna_bed_file_path = $annotation_file_path_random . $DOT . q{rrna.bed};
    my $temp_dict_file_path     = $annotation_file_path_random . $DOT . q{dict};

    say {$filehandle} q{## Getting rRNA transcripts and converting to bed format};
    perl_nae_oneliners(
        {
            filehandle     => $filehandle,
            oneliner_name  => q{get_rrna_transcripts},
            stdinfile_path => $annotation_file_path,
        }
    );
    print {$filehandle} $PIPE . $SPACE;

    gtf2bed(
        {
            filehandle      => $filehandle,
            infile_path     => $DASH,
            stdoutfile_path => $temp_rrna_bed_file_path,
        }
    );
    say {$filehandle} $NEWLINE;

    picardtools_createsequencedictionary(
        {
            filehandle => $filehandle,
            java_jar   => catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
            java_use_large_pages => $active_parameter_href->{java_use_large_pages},
            memory_allocation    => q{Xmx2g},
            outfile_path         => $temp_dict_file_path,
            referencefile_path   => $active_parameter_href->{human_genome_reference},
            temp_directory       => $active_parameter_href->{temp_directory},
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Convert bed to interval_list format};
    picardtools_bedtointervallist(
        {
            filehandle  => $filehandle,
            infile_path => $temp_rrna_bed_file_path,
            java_jar    => catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
            java_use_large_pages => $active_parameter_href->{java_use_large_pages},
            memory_allocation    => q{Xmx2g},
            outfile_path         => $temp_file_path,
            sequence_dictionary  => $temp_dict_file_path,
            temp_directory       => $active_parameter_href->{temp_directory},
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Remove temporary files};
    foreach my $temp_file ( $temp_dict_file_path, $temp_rrna_bed_file_path ) {
        gnu_rm(
            {
                filehandle  => $filehandle,
                force       => 1,
                infile_path => $temp_file,
            }
        );
        print {$filehandle} $NEWLINE;
    }

    return;
}

1;
