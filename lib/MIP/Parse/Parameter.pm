package MIP::Parse::Parameter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
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
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ parse_infiles parse_prioritize_variant_callers parse_start_with_recipe };

}

sub parse_infiles {

## Function : Collects the ".fastq(.gz)" files from the supplied infiles directory. Checks if any files exist.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $file_info_href        => File info hash {REF}
##          : $log                   => Log object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $log;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Parameter qw{ check_infiles };
    use MIP::Get::File qw{ get_files get_matching_values_key };
    use MIP::Set::File qw{ set_infiles };

    ## Collects inputfiles governed by sample_ids
  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Return the key if the hash value and query match
        my $infile_directory = get_matching_values_key(
            {
                active_parameter_href => $active_parameter_href,
                parameter_name        => q{infile_dirs},
                query_value           => $sample_id,
            }
        );

        ## Get the file(s) from filesystem
        my @infiles = get_files(
            {
                file_directory   => $infile_directory,
                rule_name        => q{*.fastq*},
                rule_skip_subdir => q{original_fastq_files},
            }
        );

        ## Check infiles found and that they contain sample_id
        check_infiles(
            {
                infiles_ref      => \@infiles,
                infile_directory => $infile_directory,
                log              => $log,
                sample_id        => $sample_id,
            }
        );

        ## Set the infile features i.e. dir and files
        set_infiles(
            {
                file_info_href   => $file_info_href,
                infiles_ref      => \@infiles,
                infile_directory => $infile_directory,
                sample_id        => $sample_id,
            }
        );

        ## Broadcast to user
        $log->info(q{Reads from platform:});
        $log->info( q{Sample id: } . $sample_id );
        $log->info(qq{\tInputfiles:});

        ## Log each file from platform
      FILE:
        foreach my $file (@infiles) {

            # Indent for visability
            $log->info( qq{\t\t}, $file );
        }
    }
    return 1;
}

sub parse_prioritize_variant_callers {

## Function : Check that all active variant callers have a prioritization order and that the prioritization elements match a supported variant caller
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $log                   => Log object
##          : $parameter_href        => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Parameter qw{ check_prioritize_variant_callers };

    my %priority_call_parameter = (
        variant_callers            => q{gatk_combinevariants_prioritize_caller},
        structural_variant_callers => q{sv_svdb_merge_prioritize},
    );

    while ( my ( $variant_caller_type, $prioritize_parameter_name ) =
        each %priority_call_parameter )
    {

        ## Check if we have any active callers
        my $activate_caller_tracker = 0;

      CALLER:
        foreach my $variant_caller (
            @{ $parameter_href->{dynamic_parameter}{$variant_caller_type} } )
        {

            if ( $active_parameter_href->{$variant_caller} ) {

                $activate_caller_tracker++;
            }
        }
        if ($activate_caller_tracker) {

            check_prioritize_variant_callers(
                {
                    active_parameter_href => $active_parameter_href,
                    log                   => $log,
                    parameter_href        => $parameter_href,
                    parameter_name        => $prioritize_parameter_name,
                    variant_callers_ref =>
                      \@{ $parameter_href->{dynamic_parameter}{$variant_caller_type} },
                }
            );
        }
        else {
            ## No active callers at all
            return;
        }
    }
    return 1;
}

sub parse_start_with_recipe {

## Function : Get initiation recipe, downstream dependencies and update recipe modes fo start_with_recipe parameter
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $initiation_file       => Initiation file for pipeline
##          : $log                   => Log object
##          : $parameter_href        => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $initiation_file;
    my $log;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        initiation_file => {
            defined     => 1,
            required    => 1,
            store       => \$initiation_file,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Parameter qw{ check_recipe_exists_in_hash };
    use MIP::File::Format::Yaml qw{ load_yaml };
    use MIP::Get::Analysis qw{ get_dependency_tree };
    use MIP::Update::Recipes qw{  update_recipe_mode_with_start_with };

    return if ( not defined $active_parameter_href->{start_with_recipe} );

    my %dependency_tree = load_yaml(
        {
            yaml_file => $initiation_file,
        }
    );

    check_recipe_exists_in_hash(
        {
            log            => $log,
            parameter_name => $active_parameter_href->{start_with_recipe},
            query_ref      => \$active_parameter_href->{start_with_recipe},
            truth_href     => $parameter_href,
        }
    );

    my @start_with_recipes;
    my $is_recipe_found = 0;
    my $is_chain_found  = 0;

    ## Collects all downstream recipes from initation point
    get_dependency_tree(
        {
            dependency_tree_href   => \%dependency_tree,
            is_recipe_found_ref    => \$is_recipe_found,
            is_chain_found_ref     => \$is_chain_found,
            recipe                 => $active_parameter_href->{start_with_recipe},
            start_with_recipes_ref => \@start_with_recipes,
        }
    );

    ## Update recipe mode depending on start with flag
    update_recipe_mode_with_start_with(
        {
            active_parameter_href  => $active_parameter_href,
            recipes_ref            => \@{ $parameter_href->{dynamic_parameter}{recipe} },
            start_with_recipes_ref => \@start_with_recipes,
        }
    );
    return 1;
}

1;
