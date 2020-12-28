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

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Analysis::Bwa_mem} => [qw{ analysis_bwa_mem analysis_run_bwa_mem }],
        q{MIP::Set::Analysis}              => [qw{ set_recipe_bwa_mem }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::Analysis qw{ set_recipe_bwa_mem };
use MIP::Recipes::Analysis::Bwa_mem qw{ analysis_bwa_mem analysis_run_bwa_mem };

diag(   q{Test set_recipe_bwa_mem from Analysis.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $GENOME_BUILD_VERSION_38 => 38;
Readonly my $GENOME_BUILD_VERSION_37 => 37;
Readonly my $GENOME_BUILD_VERSION_19 => 19;

## Given genome build version 19
my %analysis_recipe;
set_recipe_bwa_mem(
    {
        analysis_recipe_href           => \%analysis_recipe,
        human_genome_reference_version => $GENOME_BUILD_VERSION_19,
        run_bwakit                     => 1,
    }
);
my %expected_analysis_recipe = ( bwa_mem => \&analysis_bwa_mem, );

## Then set bwa mem recipe
is_deeply( \%analysis_recipe, \%expected_analysis_recipe,
    q{Set bwa mem recipe for hg} . $GENOME_BUILD_VERSION_19 );

## Given genome build version 37
set_recipe_bwa_mem(
    {
        analysis_recipe_href           => \%analysis_recipe,
        human_genome_reference_version => $GENOME_BUILD_VERSION_37,
        run_bwakit                     => 1,
    }
);
$expected_analysis_recipe{bwa_mem} = \&analysis_bwa_mem;

## Then set bwa mem recipe
is_deeply( \%analysis_recipe, \%expected_analysis_recipe,
    q{Set bwa mem recipe for grch} . $GENOME_BUILD_VERSION_37 );

## Given genome build version 38
set_recipe_bwa_mem(
    {
        analysis_recipe_href           => \%analysis_recipe,
        human_genome_reference_version => $GENOME_BUILD_VERSION_38,
        run_bwakit                     => 1,
    }
);
$expected_analysis_recipe{bwa_mem} = \&analysis_run_bwa_mem;

## Then set run-bwa mem recipe
is_deeply( \%analysis_recipe, \%expected_analysis_recipe,
    q{Set run-bwa mem recipe for grch} . $GENOME_BUILD_VERSION_38 );

## Given genome build version 38 and bwakit turned off
set_recipe_bwa_mem(
    {
        analysis_recipe_href           => \%analysis_recipe,
        human_genome_reference_version => $GENOME_BUILD_VERSION_38,
        run_bwakit                     => 0,
    }
);
$expected_analysis_recipe{bwa_mem} = \&analysis_bwa_mem;

## Then set run-bwa mem recipe
is_deeply( \%analysis_recipe, \%expected_analysis_recipe,
    q{Set bwa mem recipe for grch} . $GENOME_BUILD_VERSION_38 );
done_testing();
