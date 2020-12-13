package MIP::Io::Recipes;

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

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ build_file_prefix_tag set_file_prefix_tag };
}

sub build_file_prefix_tag {

## Function : Build the file tags depending on which recipes are used by the user to relevant chain.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $case_id               => Case id
##          : $file_info_href        => Info on files hash {REF}
##          : $order_recipes_ref     => Order of addition to parameter array {REF}
##          : $parameter_href        => Parameter hash {REF}
##          : $sample_ids_ref        => Sample ids {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $case_id;
    my $file_info_href;
    my $order_recipes_ref;
    my $parameter_href;
    my $sample_ids_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id => {
            defined     => 1,
            required    => 1,
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
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Io::Recipes qw{ set_file_prefix_tag };
    use MIP::Parameter qw{ get_parameter_attribute };

    ## Used to enable seqential build-up of file_tags between analysis recipes
    my %temp_file_ending;

  PARAMETER:
    foreach my $recipe_name ( @{$order_recipes_ref} ) {

        ## Only active parameters
        next PARAMETER
          if ( not defined $active_parameter_href->{$recipe_name} );

        ## Get parameter recipe attributes and unpack
        my %recipe_attribute = get_parameter_attribute(
            {
                parameter_href => $parameter_href,
                parameter_name => $recipe_name,
            }
        );
        my $current_chain = $recipe_attribute{chain};
        my $file_tag      = $recipe_attribute{file_tag};

        ## Skip parameters with no file tag
        next PARAMETER if ( $file_tag eq q{nofile_tag} );

        ### Per sample_id
      SAMPLE_ID:
        foreach my $sample_id ( @{$sample_ids_ref} ) {

            ## Set file tag for recipe and update sequential build-up of fileending
            $temp_file_ending{$current_chain}{$sample_id} = set_file_prefix_tag(
                {
                    current_chain         => $current_chain,
                    file_tag              => $file_tag,
                    file_info_href        => $file_info_href,
                    id                    => $sample_id,
                    is_active_recipe      => $active_parameter_href->{$recipe_name},
                    recipe_name           => $recipe_name,
                    temp_file_ending_href => \%temp_file_ending,
                }
            );
        }

        ### Per case_id
        ## Set file tag for recipe and update sequential build-up of fileending
        $temp_file_ending{$current_chain}{$case_id} = set_file_prefix_tag(
            {
                current_chain         => $current_chain,
                file_tag              => $file_tag,
                file_info_href        => $file_info_href,
                id                    => $case_id,
                is_active_recipe      => $active_parameter_href->{$recipe_name},
                recipe_name           => $recipe_name,
                temp_file_ending_href => \%temp_file_ending,
            }
        );
    }
    return;
}

sub set_file_prefix_tag {

## Function : Set the file tag depending on active recipes
## Returns  : $file_tag_to_set
## Arguments: $current_chain         => Name of current chain
##          : $file_info_href        => Info on files hash {REF}
##          : $file_tag              => File tag to set
##          : $id                    => To change id for case or sample
##          : $is_active_recipe      => Active recipe for this analysis
##          : $recipe_name           => Recipe to add file tag for
##          : $temp_file_ending_href => Store sequential build of file tag

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $current_chain;
    my $file_info_href;
    my $file_tag;
    my $id;
    my $is_active_recipe;
    my $recipe_name;
    my $temp_file_ending_href;

    my $tmpl = {
        current_chain => {
            defined     => 1,
            required    => 1,
            store       => \$current_chain,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_tag => {
            defined     => 1,
            required    => 1,
            store       => \$file_tag,
            strict_type => 1,
        },
        id => {
            defined     => 1,
            required    => 1,
            store       => \$id,
            strict_type => 1,
        },
        is_active_recipe => {
            defined     => 1,
            required    => 1,
            store       => \$is_active_recipe,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        temp_file_ending_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$temp_file_ending_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info qw{ set_file_tag };

    ## Default from previous recipes
    my $file_tag_to_set = $temp_file_ending_href->{$current_chain}{$id};

    ## File_ending should be added for this recipe
    if ($is_active_recipe) {

        ## Inherit from MAIN if new branch and first recipe
        $file_tag_to_set = _inherit_chain_main(
            {
                current_chain         => $current_chain,
                id                    => $id,
                temp_file_ending_href => $temp_file_ending_href,
            }
        );

        if ( defined $file_tag_to_set ) {

            ## Add new file tag to sequential build-up
            $file_tag_to_set .= $file_tag;
        }
        else {
            ## First recipe that have a filending

            $file_tag_to_set = $file_tag;
        }
    }

    ## Set the file tag
    set_file_tag(
        {
            file_info_href => $file_info_href,
            file_tag       => $file_tag_to_set,
            id             => $id,
            recipe_name    => $recipe_name,
        }
    );

    return $file_tag_to_set;
}

sub _inherit_chain_main {

## Function : Inherit file tags from MAIN chain
## Returns  :
## Arguments: $current_chain         => Name of current chain
##          : $id                    => To change id for case or sample
##          : $temp_file_ending_href => Store sequential build of file tag

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $current_chain;
    my $id;
    my $temp_file_ending_href;

    my $tmpl = {
        current_chain => {
            defined     => 1,
            required    => 1,
            store       => \$current_chain,
            strict_type => 1,
        },
        id => {
            defined     => 1,
            required    => 1,
            store       => \$id,
            strict_type => 1,
        },
        temp_file_ending_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$temp_file_ending_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Unpack
    my $branch_file_tag = $temp_file_ending_href->{$current_chain}{$id};
    my $main_file_tag   = $temp_file_ending_href->{MAIN}{$id};

    ## Return if already on MAIN branch
    return $branch_file_tag if ( $current_chain eq q{MAIN} );

    ## Check if first recipe on other branch
    return $branch_file_tag if ( defined $branch_file_tag );

    ## First recipe on branch other than MAIN
    ## Inherit current MAIN chain file tag
    $temp_file_ending_href->{$current_chain}{$id} = $main_file_tag;
    return $main_file_tag;
}

1;
