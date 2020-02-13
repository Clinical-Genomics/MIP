package MIP::Dependency_tree;

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
use List::MoreUtils qw{ any firstidx };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_dependency_tree
      get_dependency_tree_chain
      get_dependency_subtree
      get_recipes_for_dependency_tree_chain
      get_recipe_dependency_tree_chain
      set_dependency_tree_order
    };
}

sub get_dependency_tree {

## Function  : Collects all downstream recipes from initation point.
## Returns   :
## Arguments : $current_chain          => Current chain
##           : $is_recipe_found_ref    => Found initiation recipe {REF}
##           : $is_chain_found_ref     => Found recipe chain
##           : $recipe                 => Initiation point
##           : $start_with_recipes_ref => Store recipes
##           : $dependency_tree_href   => Dependency hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $current_chain;
    my $is_recipe_found_ref;
    my $is_chain_found_ref;
    my $recipe;
    my $start_with_recipes_ref;
    my $dependency_tree_href;

    my $tmpl = {
        current_chain => {
            store       => \$current_chain,
            strict_type => 1,
        },
        is_recipe_found_ref => {
            default     => \$$,
            store       => \$is_recipe_found_ref,
            strict_type => 1,
        },
        is_chain_found_ref => {
            default     => \$$,
            store       => \$is_chain_found_ref,
            strict_type => 1,
        },
        recipe => {
            required    => 1,
            store       => \$recipe,
            strict_type => 1,
        },
        start_with_recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$start_with_recipes_ref,
            strict_type => 1,
        },
        dependency_tree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$dependency_tree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Copy hash to enable recursive removal of keys
    my %tree = %{$dependency_tree_href};

  KEY_VALUE_PAIR:
    while ( my ( $key, $value ) = each %tree ) {

        ## Do not enter into more chains than one if recipe and chain is found
        next KEY_VALUE_PAIR
          if ( $key =~ /CHAIN_/sxm
            && ${$is_recipe_found_ref}
            && ${$is_chain_found_ref} ne q{CHAIN_MAIN} );

        ## Do not add recipe name or PARALLEL
        if ( $key =~ /CHAIN_/sxm ) {

            $current_chain = $key;
        }

        if ( ref $value eq q{ARRAY} ) {
            ## Inspect element

          ELEMENT:
            foreach my $element ( @{$value} ) {

                ## Call recursive
                if ( ref $element eq q{HASH} ) {

                    get_dependency_tree(
                        {
                            current_chain          => $current_chain,
                            dependency_tree_href   => $element,
                            is_recipe_found_ref    => $is_recipe_found_ref,
                            is_chain_found_ref     => $is_chain_found_ref,
                            recipe                 => $recipe,
                            start_with_recipes_ref => $start_with_recipes_ref,
                        }
                    );
                }
                ## Found initiator recipe
                if ( ref $element ne q{HASH}
                    && $element eq $recipe )
                {

                    ## Start collecting recipes downstream
                    ${$is_recipe_found_ref} = 1;

                    ## Found chain that recipe belongs to
                    # Set is part of chain signal
                    ${$is_chain_found_ref} = $current_chain;

                }

                ## Special case for parallel section
                if ( $key eq q{PARALLEL}
                    && ${$is_recipe_found_ref} )
                {

                    if ( any { $_ eq $recipe } @{$value} ) {

                        ## Add only start_with recipe from parallel section
                        push @{$start_with_recipes_ref}, $recipe;

                        ## Skip any remaining hash_ref or element
                        last ELEMENT;
                    }
                }

                ## Add downstream recipes
                if ( ref $element ne q{HASH}
                    && ${$is_recipe_found_ref} )
                {

                    push @{$start_with_recipes_ref}, $element;
                }
            }
        }

        ## Remove identifier
        delete $tree{$key};
    }
    return;
}

sub get_dependency_tree_chain {

## Function  : Sets chain id to parameters hash from the dependency tree
## Returns   :
## Arguments : $current_chain        => Current chain
##           : $dependency_tree_href => Dependency hash {REF}
##           : $parameter_href       => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $current_chain;
    my $dependency_tree_href;
    my $parameter_href;

    my $tmpl = {
        current_chain => {
            store       => \$current_chain,
            strict_type => 1,
        },
        dependency_tree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$dependency_tree_href,
            strict_type => 1,
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

    ## Copy hash to enable recursive removal of keys
    my %tree = %{$dependency_tree_href};

  KEY_VALUE_PAIR:
    while ( my ( $key, $value ) = each %tree ) {

        ## Add ID of chain
        my ($chain_id) = $key =~ /CHAIN_(\S+)/sxm;

        ## If chain_id is found
        if ( defined $chain_id ) {

            ## Set current chain
            $current_chain = $chain_id;
        }

        ## Call recursive
        if ( ref $value eq q{HASH} ) {

            get_dependency_tree_chain(
                {
                    current_chain        => $current_chain,
                    dependency_tree_href => $value,
                    parameter_href       => $parameter_href,
                }
            );
        }
        elsif ( ref $value eq q{ARRAY} ) {
            ## Inspect element

          ELEMENT:
            foreach my $element ( @{$value} ) {

                ## Call recursive
                if ( ref $element eq q{HASH} ) {

                    get_dependency_tree_chain(
                        {
                            current_chain        => $current_chain,
                            dependency_tree_href => $element,
                            parameter_href       => $parameter_href,
                        }
                    );
                }

                ## Found recipes
                if ( ref $element ne q{HASH} ) {

                    $parameter_href->{$element}{chain} = $current_chain;
                }

                if ( $key eq q{PARALLEL} ) {

                    $parameter_href->{$element}{chain} = uc $element;
                }

                ## Hash in PARALLEL section create anonymous chain ID
                ## E.g. haplotypecaller->genotypegvcfs
                if ( $key eq uc $element ) {

                  RECIPE:
                    foreach my $recipe ( @{$value} ) {

                        $parameter_href->{$recipe}{chain} = uc $element;
                    }
                    last ELEMENT;
                }
            }
        }

        ## Remove identifier
        delete $tree{$key};
    }
    return;
}

sub get_dependency_subtree {

## Function : Get part of dependency tree
## Returns  : %dependency_tree
## Arguments: $chain_initiation_point  => Chain to operate on
##          : $dependency_tree_href    => Dependency hash {REF}
##          : $dependency_subtree_href => Dependency sub hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $chain_initiation_point;
    my $dependency_tree_href;
    my $dependency_subtree_href;

    my $tmpl = {
        chain_initiation_point => {
            store       => \$chain_initiation_point,
            strict_type => 1,
        },
        dependency_tree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$dependency_tree_href,
            strict_type => 1,
        },
        dependency_subtree_href => {
            default     => {},
            required    => 1,
            store       => \$dependency_subtree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Return if tree is found
    return if ( defined $dependency_subtree_href->{$chain_initiation_point} );

    ## Copy hash to enable recursive removal of keys
    my %tree = %{$dependency_tree_href};

  KEY_VALUE_PAIR:
    while ( my ( $key, $value ) = each %tree ) {

        ## Save subtree if it matches chain
        if ( $key eq $chain_initiation_point ) {
            $dependency_subtree_href->{$chain_initiation_point} = $value;
        }

        ## Inspect element
        if ( ref $value eq q{ARRAY} ) {

          ELEMENT:
            foreach my $element ( @{$value} ) {

                ## Call recursive
                if ( ref $element eq q{HASH} ) {

                    get_dependency_subtree(
                        {
                            dependency_tree_href    => $element,
                            dependency_subtree_href => $dependency_subtree_href,
                            chain_initiation_point  => $chain_initiation_point,
                        }
                    );
                }
            }
            delete $tree{$key};
        }
    }
    return;
}

sub set_dependency_tree_order {

## Function  : Collects order of all recipes from initiation.
## Returns   :
## Arguments : $recipes_ref          => Recipes {REF}
##           : $dependency_tree_href => Dependency hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $recipes_ref;
    my $dependency_tree_href;

    my $tmpl = {
        recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$recipes_ref,
            strict_type => 1,
        },
        dependency_tree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$dependency_tree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Copy hash to enable recursive removal of keys
    my %tree = %{$dependency_tree_href};

  KEY_VALUE_PAIR:
    while ( my ( $key, $value ) = each %tree ) {

        ## Call recursive
        if ( ref $value eq q{HASH} ) {

            set_dependency_tree_order(
                {
                    dependency_tree_href => $value,
                    recipes_ref          => $recipes_ref,
                }
            );
        }
        elsif ( ref $value eq q{ARRAY} ) {
            ## Inspect element

          ELEMENT:
            foreach my $element ( @{$value} ) {

                ## Call recursive
                if ( ref $element eq q{HASH} ) {

                    set_dependency_tree_order(
                        {
                            dependency_tree_href => $element,
                            recipes_ref          => $recipes_ref,
                        }
                    );
                }
                ## Found recipe
                if ( ref $element ne q{HASH} ) {

                    ## Add to order
                    push @{$recipes_ref}, $element;
                }
            }
        }

        ## Remove identifier
        delete $tree{$key};
    }
    return;
}

sub get_recipes_for_dependency_tree_chain {

## Function  : Collects all recipes downstream of initation point
## Returns   : @chain_recipes
## Arguments : $chain_initiation_point  => Chain to operate on
##           : $dependency_tree_href    => Dependency hash {REF}
##           : $recipe_initiation_point => Recipe to start with

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $chain_initiation_point;
    my $dependency_tree_href;
    my $recipe_initiation_point;

    my $tmpl = {
        chain_initiation_point => {
            defined     => 1,
            required    => 1,
            store       => \$chain_initiation_point,
            strict_type => 1,
        },
        dependency_tree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$dependency_tree_href,
            strict_type => 1,
        },
        recipe_initiation_point => {
            defined     => 1,
            store       => \$recipe_initiation_point,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Get the dependency subtree
    my $dependency_subtree_href = {};
    get_dependency_subtree(
        {
            dependency_subtree_href => $dependency_subtree_href,
            dependency_tree_href    => $dependency_tree_href,
            chain_initiation_point  => $chain_initiation_point,
        }
    );

    ## Get the recipes
    my @recipes;
    set_dependency_tree_order(
        {
            dependency_tree_href => $dependency_subtree_href,
            recipes_ref          => \@recipes,
        }
    );

    ## Slice if $recipe_initiation_point is defined
    if ($recipe_initiation_point) {
        my $initiation_idx = firstidx { $_ eq $recipe_initiation_point } @recipes;
        @recipes = @recipes[ $initiation_idx .. $#recipes ];

    }
    return @recipes;
}

sub get_recipe_dependency_tree_chain {

## Function  : Get the chain to which a recipe belongs
## Returns   :
## Arguments : $chain_id_ref         => Chain found {REF}
##           : $current_chain        => Current chain
##           : $dependency_tree_href => Dependency hash {REF}
##           : $recipe               => Initiation point

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $chain_id_ref;
    my $current_chain;
    my $dependency_tree_href;
    my $recipe;

    my $tmpl = {
        chain_id_ref => {
            default     => \$$,
            store       => \$chain_id_ref,
            strict_type => 1,
        },
        current_chain => {
            store       => \$current_chain,
            strict_type => 1,
        },
        dependency_tree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$dependency_tree_href,
            strict_type => 1,
        },
        recipe => {
            required    => 1,
            store       => \$recipe,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Return if chain has been found
    return if ( ${$chain_id_ref} );

    ## Copy hash to enable recursive removal of keys
    my %tree = %{$dependency_tree_href};

  KEY_VALUE_PAIR:
    while ( my ( $key, $value ) = each %tree ) {

        ## Don't store PARALLEL as the current chain
        if ( $key =~ /CHAIN_/sxm ) {

            $current_chain = $key;
        }

        ## Call recursive
        if ( ref $value eq q{HASH} ) {

            get_recipe_dependency_tree_chain(
                {
                    chain_id_ref         => $chain_id_ref,
                    current_chain        => $current_chain,
                    dependency_tree_href => $value,
                    recipe               => $recipe,
                }
            );
        }
        ## Inspect element
        if ( ref $value eq q{ARRAY} ) {

          ELEMENT:
            foreach my $element ( @{$value} ) {

                ## Call recursive
                if ( ref $element eq q{HASH} ) {

                    get_recipe_dependency_tree_chain(
                        {
                            chain_id_ref         => $chain_id_ref,
                            current_chain        => $current_chain,
                            dependency_tree_href => $element,
                            recipe               => $recipe,
                        }
                    );
                }
                ## Found recipe
                if ( ( ref $element ne q{HASH} ) && ( $element eq $recipe ) ) {

                    ## Save current chain
                    ${$chain_id_ref} = $current_chain;

                    last ELEMENT;
                }
            }
            delete $tree{$key};
        }
    }
    return;
}

1;
