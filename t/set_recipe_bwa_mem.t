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
        q{MIP::Recipes::Analysis::Bwa_mem} =>
          [qw{ analysis_bwa_mem analysis_run_bwa_mem }],
        q{MIP::Set::Analysis}  => [qw{ set_recipe_bwa_mem }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::Analysis qw{ set_recipe_bwa_mem };
use MIP::Recipes::Analysis::Bwa_mem qw{ analysis_bwa_mem analysis_run_bwa_mem };

diag(   q{Test set_recipe_bwa_mem from Analysis.pm v}
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
Readonly my $GENOME_BUILD_VERSION_20 => 20;
Readonly my $GENOME_BUILD_VERSION_19 => 19;

## Given hg genome source when genome build version 19
my $reference_source = q{hg};
my %analysis_recipe;
set_recipe_bwa_mem(
    {
        analysis_recipe_href           => \%analysis_recipe,
        human_genome_reference_source  => $reference_source,
        human_genome_reference_version => $GENOME_BUILD_VERSION_19,
    }
);
my %expected_analysis_recipe = ( bwa_mem => \&analysis_bwa_mem, );
## Then set bwa mem recipe
is_deeply( \%analysis_recipe, \%expected_analysis_recipe,
    q{Set bwa mem recipe for } . $reference_source . $GENOME_BUILD_VERSION_19 );

## Given hg genome source when genome build version 20
set_recipe_bwa_mem(
    {
        analysis_recipe_href           => \%analysis_recipe,
        human_genome_reference_source  => $reference_source,
        human_genome_reference_version => $GENOME_BUILD_VERSION_20,
    }
);
$expected_analysis_recipe{bwa_mem} = \&analysis_run_bwa_mem;

## Then set run-bwa mem recipe
is_deeply( \%analysis_recipe, \%expected_analysis_recipe,
    q{Set run-bwa mem recipe for } . $reference_source . $GENOME_BUILD_VERSION_20 );

## Given grch genome source when genome build version 37
$reference_source = q{grch};
set_recipe_bwa_mem(
    {
        analysis_recipe_href           => \%analysis_recipe,
        human_genome_reference_source  => $reference_source,
        human_genome_reference_version => $GENOME_BUILD_VERSION_37,
    }
);
$expected_analysis_recipe{bwa_mem} = \&analysis_bwa_mem;

## Then set bwa mem recipe
is_deeply( \%analysis_recipe, \%expected_analysis_recipe,
    q{Set bwa mem recipe for } . $reference_source . $GENOME_BUILD_VERSION_37 );

## Given grch genome source when genome build version 38
$reference_source = q{grch};
set_recipe_bwa_mem(
    {
        analysis_recipe_href           => \%analysis_recipe,
        human_genome_reference_source  => $reference_source,
        human_genome_reference_version => $GENOME_BUILD_VERSION_38,
    }
);
$expected_analysis_recipe{bwa_mem} = \&analysis_run_bwa_mem;

## Then set run-bwa mem recipe
is_deeply( \%analysis_recipe, \%expected_analysis_recipe,
    q{Set run-bwa mem recipe for } . $reference_source . $GENOME_BUILD_VERSION_38 );

done_testing();
