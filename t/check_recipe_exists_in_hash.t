#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Check} => [qw{ check_recipe_exists_in_hash }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Check qw{ check_recipe_exists_in_hash };

diag(   q{Test check_recipe_exists_in_hash from Check.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( {} );

## Given recipe names
my %parameter = ( glnexus_merge => 1, );

my %active_parameter = (
    recipe_time => {
        bwa_mem       => 1,
        glnexus_merge => 1,
    },
    associated_recipe => [ qw{ fastqc_ar }, ],
);
## When one does not exist in truth hash
trap {
    check_recipe_exists_in_hash(
        {
            parameter_name => q{recipe_time},
            query_ref      => \%{ $active_parameter{recipe_time} },
            truth_href     => \%parameter,
        }
    )
};

## Then exist and throw error
ok( $trap->exit, q{Exit if recipe key does not exist} );
like( $trap->stderr, qr/FATAL/xms, q{Throw FATAL log message} );

## Given recipe names
%parameter = (
    glnexus_merge => 1,
    bwa_mem       => 1,
);

## When all exists in truth hash
my $return = check_recipe_exists_in_hash(
    {
        parameter_name => q{recipe_time},
        query_ref      => \%{ $active_parameter{recipe_time} },
        truth_href     => \%parameter,
    }
);
is( $return, undef, q{All recipe keys exists in truth hash} );

## Given recipe names, when none exists in truth hash
trap {
    check_recipe_exists_in_hash(
        {
            parameter_name => q{associated_recipe},
            query_ref      => \@{ $active_parameter{associated_recipe} },
            truth_href     => \%parameter,
        }
    )
};

## Then exist and throw error
ok( $trap->exit, q{Exit if recipe element does not exist} );
like( $trap->stderr, qr/FATAL/xms, q{Throw FATAL log message} );

## Given recipe names
%parameter = (
    glnexus_merge => 1,
    bwa_mem       => 1,
    fastqc_ar     => 1,
);

## When all exists in truth hash
$return = check_recipe_exists_in_hash(
    {
        parameter_name => q{associated_recipe},
        query_ref      => \@{ $active_parameter{associated_recipe} },
        truth_href     => \%parameter,
    }
);
is( $return, undef, q{All recipe element exists in truth hash} );

## Given a recipe name
my $recipe_name = q{this_recipe_does_not_exist};

%parameter = ( bwa_mem => 1, );

## When recipe does not exists in truth hash
trap {
    check_recipe_exists_in_hash(
        {
            parameter_name => $recipe_name,
            query_ref      => \$recipe_name,
            truth_href     => \%parameter,
        }
    )
};

## Then exist and throw error
ok( $trap->exit, q{Exit if recipe does not exist} );
like( $trap->stderr, qr/FATAL/xms, q{Throw FATAL log message} );

## Given a recipe name
$recipe_name = q{bwa_mem};

## When recipe exists in truth hash
$return = check_recipe_exists_in_hash(
    {
        parameter_name => $recipe_name,
        query_ref      => \$recipe_name,
        truth_href     => \%parameter,
    }
);
is( $return, undef, q{Recipe exists in truth hash} );

done_testing();
