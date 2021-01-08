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
        q{MIP::Recipes::Analysis::Deepvariant} => [qw{ analysis_deepvariant }],
        q{MIP::Recipes::Analysis::Deeptrio}    => [qw{ analysis_deeptrio }],
        q{MIP::Set::Analysis}                  => [qw{ set_recipe_deepvariant }],
        q{MIP::Test::Fixtures}                 => [qw{ test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::Analysis qw{ set_recipe_deepvariant };
use MIP::Recipes::Analysis::Deeptrio qw{ analysis_deeptrio };
use MIP::Recipes::Analysis::Deepvariant qw{ analysis_deepvariant };
use MIP::Test::Fixtures qw{ test_mip_hashes };

diag(   q{Test set_recipe_deepvariant from Analysis.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %sample_info = test_mip_hashes( { mip_hash_name => q{qc_sample_info} } );
my %analysis_recipe;

## When constellation state is a duo
$sample_info{has_duo}  = 1;
$sample_info{has_trio} = 0;

set_recipe_deepvariant(
    {
        analysis_recipe_href => \%analysis_recipe,
        sample_info_href     => \%sample_info,
    }
);
my %expected_analysis_recipe = (
    deepvariant => undef,
    deeptrio    => \&analysis_deeptrio,
);

## Then set deepvariant recipe
is_deeply( \%analysis_recipe, \%expected_analysis_recipe, q{Set deepvariant recipe for a duo} );

## When constellation state is a trio
$sample_info{has_duo}  = 0;
$sample_info{has_trio} = 1;
set_recipe_deepvariant(
    {
        analysis_recipe_href => \%analysis_recipe,
        sample_info_href     => \%sample_info,
    }
);
$expected_analysis_recipe{deeptrio}    = \&analysis_deeptrio;
$expected_analysis_recipe{deepvariant} = undef;

## Then set deepvariant recipe
is_deeply( \%analysis_recipe, \%expected_analysis_recipe, q{Set deepvariant recipe for a trio} );

## When constellation state is neither a duo or trio
$sample_info{has_duo}  = 0;
$sample_info{has_trio} = 0;
set_recipe_deepvariant(
    {
        analysis_recipe_href => \%analysis_recipe,
        sample_info_href     => \%sample_info,
    }
);
$expected_analysis_recipe{deeptrio}    = undef;
$expected_analysis_recipe{deepvariant} = \&analysis_deepvariant;

## Then set bwa mem recipe
is_deeply( \%analysis_recipe, \%expected_analysis_recipe,
    q{Set deepvariant recipe when the case is neither a duo or a trio} );

done_testing();
