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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Active_parameter} => [qw{ add_recipe_bind_paths }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ add_recipe_bind_paths };

diag(   q{Test add_recipe_bind_paths from Active_parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my @export_bind_paths = qw{ path1 path2 };
my %active_parameter  = (
    recipe_bind_path => {
        recipe_1 => [qw{ path3 }],
    },
);

## Given no recipe_bind_path
add_recipe_bind_paths(
    {
        active_parameter_href => \%active_parameter,
        export_bind_paths_ref => \@export_bind_paths,
        recipe_name           => q{recipe_2},
    }
);

## Then return export_bind_paths unchanged
my @expected_bind_paths = qw{ path1 path2 };
is_deeply( \@export_bind_paths, \@expected_bind_paths, q{No extra bind paths } );

## Given extra recipe_bind_path
add_recipe_bind_paths(
    {
        active_parameter_href => \%active_parameter,
        export_bind_paths_ref => \@export_bind_paths,
        recipe_name           => q{recipe_1},
    }
);

## Then return export_bind_paths unchanged
@expected_bind_paths = qw{ path1 path2 path3};
is_deeply( \@export_bind_paths, \@expected_bind_paths, q{Add extra bind paths } );

done_testing();
