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
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Recipes::Analysis::Cadd} => [qw{ analysis_cadd analysis_cadd_gb_38 }],
        q{MIP::Set::Analysis}           => [qw{ set_recipe_cadd }],
        q{MIP::Test::Fixtures}          => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::Analysis qw{ set_recipe_cadd };
use MIP::Recipes::Analysis::Cadd qw{ analysis_cadd analysis_cadd_gb_38 };

diag(   q{Test set_recipe_cadd from Analysis.pm v}
      . $MIP::Set::Analysis::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $GENOME_BUILD_VERSION_38 => 38;
Readonly my $GENOME_BUILD_VERSION_37 => 37;

## Given genome build when version is lesser than 38
my %analysis_recipe;
set_recipe_cadd(
    {
        analysis_recipe_href           => \%analysis_recipe,
        human_genome_reference_version => $GENOME_BUILD_VERSION_37,
    }
);
my %expected_analysis_recipe = ( cadd_ar => \&analysis_cadd, );

## Then set cadd recipe
is_deeply( \%analysis_recipe, \%expected_analysis_recipe,
    q{Set cadd recipe for } . $GENOME_BUILD_VERSION_37 );

## Given genome build when version is 38
set_recipe_cadd(
    {
        analysis_recipe_href           => \%analysis_recipe,
        human_genome_reference_version => $GENOME_BUILD_VERSION_38,
    }
);
$expected_analysis_recipe{cadd_ar} = \&analysis_cadd_gb_38;

## Then set cadd_gb_38 recipe
is_deeply( \%analysis_recipe, \%expected_analysis_recipe,
    q{Set cadd_gb_38 recipe for } . $GENOME_BUILD_VERSION_38 );

done_testing();
