package MIP::Recipes::Analysis::Bamcalibrationblock;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
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
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_bamcalibrationblock };

}

## Constants
Readonly my $CLOSE_BRACKET => q{]};
Readonly my $NEWLINE       => qq{\n};
Readonly my $OPEN_BRACKET  => q{[};
Readonly my $TAB           => qq{\t};
Readonly my $UNDERSCORE    => q{_};

sub analysis_bamcalibrationblock {

## Function : Run consecutive analysis recipes on the same node
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $bamcal_ar_href          => Analysis recipes for bamcalibration block
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $log                     => Log object to write to
##          : $order_recipes_ref       => Order of recipes
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $bamcal_ar_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $log;
    my $order_recipes_ref;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        bamcal_ar_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$bamcal_ar_href,
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
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        order_recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_recipes_ref,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
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

    use MIP::Log::MIP_log4perl qw{ log_display_recipe_for_user };
    use MIP::Get::Parameter qw{ get_recipe_parameters };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $PROCESS_TIME => 80;
    Readonly my $CORE_NUMBER  => $active_parameter_href->{max_cores_per_node};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Broadcasting
  PROGRAM:
    foreach my $recipe ( @{$order_recipes_ref} ) {

        ## Only for active recipes
        next PROGRAM if ( not $active_parameter_href->{$recipe} );

        ## For displaying
        log_display_recipe_for_user(
            {
                indent_level => 1,
                log          => $log,
                recipe       => $recipe,
            }
        );
    }

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        my $xargs_file_counter = 0;
        my ( $c_n, $t, @source_environment_cmds ) = get_recipe_parameters(
            {
                active_parameter_href => $active_parameter_href,
                recipe_name           => $recipe_name,
            }
        );

        ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
        my ( $file_path, $recipe_info_path ) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                core_number           => $CORE_NUMBER,
                directory_id          => $sample_id,
                FILEHANDLE            => $FILEHANDLE,
                job_id_href           => $job_id_href,
                log                   => $log,
                recipe_directory      => $active_parameter_href->{outaligner_dir},
                recipe_name           => $recipe_name,
                process_time          => $PROCESS_TIME,
                source_environment_commands_ref => \@source_environment_cmds,
                temp_directory                  => $temp_directory,
            }
        );

      PROGRAM:
        foreach my $recipe ( @{$order_recipes_ref} ) {

            ## Only for active recipes
            next PROGRAM if ( not $active_parameter_href->{$recipe} );

            ($xargs_file_counter) = $bamcal_ar_href->{$recipe}->(
                {
                    active_parameter_href   => $active_parameter_href,
                    file_info_href          => $file_info_href,
                    file_path               => $file_path,
                    FILEHANDLE              => $FILEHANDLE,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    parameter_href          => $parameter_href,
                    recipe_info_path        => $recipe_info_path,
                    recipe_name             => $recipe,
                    sample_id               => $sample_id,
                    sample_info_href        => $sample_info_href,
                    xargs_file_counter      => $xargs_file_counter,
                }
            );
        }
    }

    return;
}

1;
