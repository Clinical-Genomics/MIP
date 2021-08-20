#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Parse} => [qw{ parse_initiation_recipe }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Parse qw{ parse_initiation_recipe };

diag(   q{Test parse_initiation_recipe from Parse.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( {} );

## Given no defined start_with_recipe or start_after_recipe parameter
my %active_parameter;
my %parameter = (
    bwa_mem        => { default => 0, },
    markduplicates => { default => 0, },
);
my %dependency_tree = test_mip_hashes( { mip_hash_name => q{dependency_tree_dna} } );
$parameter{dependency_tree_href} = \%dependency_tree;

my $return = parse_initiation_recipe(
    {
        active_parameter_href => \%active_parameter,
        parameter_href        => \%parameter,
    },
);

## Then skip parsing
is( $return, undef, q{Skip parsing} );

## Given start_with_recipe parameter, when defined
$active_parameter{start_with_recipe} = q{bwa_mem};

my $is_ok = parse_initiation_recipe(
    {
        active_parameter_href => \%active_parameter,
        parameter_href        => \%parameter,
    },
);

## Then return true for successful parsing
ok( $is_ok, q{Parsed programs from start_with_flag} );

## Given start_after_recipe parameter, when defined
$active_parameter{start_after_recipe} = q{markduplicates};
$active_parameter{start_with_recipe}  = undef;

$is_ok = parse_initiation_recipe(
    {
        active_parameter_href => \%active_parameter,
        parameter_href        => \%parameter,
    },
);

## Then return true for successful parsing
ok( $is_ok, q{Parsed programs from start_after_flag} );

done_testing();
