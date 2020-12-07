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
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Parameter}      => [qw{ parse_parameter_recipe_names }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parameter qw{ parse_parameter_recipe_names };

diag(   q{Test parse_parameter_recipe_names from Parameter.pm v}
      . $MIP::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( {} );

## Given recipe names
my %parameter = (
    markduplicates_picardtools_opt_dup_dist => {
        associated_recipe => [ qw{ mip }, ],
        data_type         => q{SCALAR},
        default           => 1,
        type              => q{recipe_argument}
    },
    mip => {
        associated_recipe => [ qw{ mip }, ],
        data_type         => q{SCALAR},
        default           => 1,
        type              => q{mip}
    },
);

## When recipe name exist in truth hash
my $is_ok = parse_parameter_recipe_names(
    {
        parameter_href => \%parameter,
    }
);

## Then return true
ok( $is_ok, q{Parsed recipe names in parameter} );

done_testing();
