package MIP::Recipes::Build::Transcript_annotation_prerequisites;

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
    our $VERSION = 1.00;

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
##          : $infile_lane_prefix_href      => Infile(s) without the ".ending" {REF}
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
    my $infile_lane_prefix_href;
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
    use MIP::Program::Gnu::Coreutils qw{ gnu_rm };
    use MIP::Language::Perl qw{ perl_nae_oneliners };
    use MIP::Language::Shell qw{ check_exist_and_move_file };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Ucsc qw{ ucsc_gtf_to_genepred };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $MAX_RANDOM_NUMBER => 100_00;

    my $recipe_file_path;
    my $submit_switch;

    ## Unpack parameters
    my $refflat_suffix  = $parameter_build_suffixes_ref->[0];
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => q{mip},
        }
    );
    my $annotation_file_path = $active_parameter_href->{transcript_annotation};

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

    $log->warn( q{Will try to create required }
          . $annotation_file_path
          . q{ associated file(s) before executing }
          . $recipe_name );

    ## Set file names
    my $annotation_file_path_random =
      $annotation_file_path . $UNDERSCORE . $random_integer;
    my $temp_genepred_file_path = $annotation_file_path_random . $DOT . q{genePred};
    my $temp_refflat_file_path  = $annotation_file_path_random . $refflat_suffix;
    my $intended_file_path      = $annotation_file_path . $refflat_suffix;

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
            stdoutfile_path => $temp_refflat_file_path,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
    check_exist_and_move_file(
        {
            filehandle          => $filehandle,
            intended_file_path  => $intended_file_path,
            temporary_file_path => $temp_refflat_file_path,
        }
    );

    ## Remove temporary files
    say {$filehandle} q{#Remove temporary files};
    gnu_rm(
        {
            filehandle  => $filehandle,
            force       => 1,
            infile_path => $temp_genepred_file_path,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Only create once
    $parameter_href->{transcript_annotation_file_endings}{build_file} = 0;

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

1;
