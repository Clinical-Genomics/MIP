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

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ set_file_prefix_tag };
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
