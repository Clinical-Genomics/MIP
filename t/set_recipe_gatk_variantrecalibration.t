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
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Set::Analysis}  => [qw{ set_recipe_gatk_variantrecalibration }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::Analysis qw{ set_recipe_gatk_variantrecalibration };
use MIP::Recipes::Analysis::Gatk_cnnscorevariants qw{ analysis_gatk_cnnscorevariants };

diag(   q{Test set_recipe_gatk_variantrecalibration from Analysis.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 1, } );

## Given more than a single sample
my %analysis_recipe;
my $gatk_cnnscorevariants = 1;
my @sample_ids            = qw{ sample_1 sample_2 };

set_recipe_gatk_variantrecalibration(
    {
        analysis_recipe_href => \%analysis_recipe,
        log                  => $log,
        sample_ids_ref       => \@sample_ids,
        use_cnnscorevariants => $gatk_cnnscorevariants,
    }
);

## Then skip setting cnnscorevariants recipe in hash
is( $analysis_recipe{gatk_variantrecalibration},
    undef, q{Skip when more than single sample} );

## Given single sample, when not to use cnnscorevariants
$gatk_cnnscorevariants = 0;
pop @sample_ids;

set_recipe_gatk_variantrecalibration(
    {
        analysis_recipe_href => \%analysis_recipe,
        log                  => $log,
        sample_ids_ref       => \@sample_ids,
        use_cnnscorevariants => $gatk_cnnscorevariants,
    }
);

## Then skip setting cnnscorevariants recipe in hash
is( $analysis_recipe{gatk_variantrecalibration},
    undef, q{Skip when not using cnnscorevariants} );

## Given single sample, when turned on cnnscorevariants
$gatk_cnnscorevariants = 1;

set_recipe_gatk_variantrecalibration(
    {
        analysis_recipe_href => \%analysis_recipe,
        log                  => $log,
        sample_ids_ref       => \@sample_ids,
        use_cnnscorevariants => $gatk_cnnscorevariants,
    }
);

## Then set cnnscorevariants recipe in hash
my %expected_recipe = ( gatk_variantrecalibration => \&analysis_gatk_cnnscorevariants, );
is_deeply( \%analysis_recipe, \%expected_recipe, q{Set cnnscorevariants recipe in hash} );

done_testing();
