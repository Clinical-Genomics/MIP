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


## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Set::Analysis}  => [qw{ set_recipe_on_analysis_type }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::Analysis qw{ set_recipe_on_analysis_type };
use MIP::Recipes::Analysis::Mip_vcfparser
  qw{ analysis_mip_vcfparser_sv_wes analysis_mip_vcfparser_sv_wgs };

diag(   q{Test set_recipe_on_analysis_type from Analysis.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given no recipes and wes
my $consensus_analysis_type = q{wes};
my %analysis_recipe;
my %expected_recipe;

## Update which recipe to use depending on consensus analysis type
set_recipe_on_analysis_type(
    {
        consensus_analysis_type => $consensus_analysis_type,
        analysis_recipe_href    => \%analysis_recipe,
    }
);

## Then set wgs recipe
is_deeply( \%analysis_recipe, \%expected_recipe, q{Set no recipe when no program} );

## Given program name and wes
$analysis_recipe{sv_vcfparser} = undef;
my %expected_wes_recipe = ( sv_vcfparser => \&analysis_mip_vcfparser_sv_wes, );
## Update which recipe to use depending on consensus analysis type
set_recipe_on_analysis_type(
    {
        consensus_analysis_type => $consensus_analysis_type,
        analysis_recipe_href    => \%analysis_recipe,
    }
);

## Then set wgs recipe
is(
    $analysis_recipe{sv_vcfparser},
    $expected_wes_recipe{sv_vcfparser},
    q{Set wes recipe}
);

## Given an existing wes recipes and wgs
$consensus_analysis_type = q{wgs};
my %expected_wgs_recipe = ( sv_vcfparser => \&analysis_mip_vcfparser_sv_wgs, );
## Update which recipe to use depending on consensus analysis type
set_recipe_on_analysis_type(
    {
        consensus_analysis_type => $consensus_analysis_type,
        analysis_recipe_href    => \%analysis_recipe,
    }
);

## Then wgs recipe should be set instead of wes
is(
    $analysis_recipe{sv_vcfparser},
    $expected_wgs_recipe{sv_vcfparser},
    q{Set wgs recipe}
);

## Given a mixed consensus analysis type
$consensus_analysis_type = q{mixed};
$analysis_recipe{sv_vcfparser} = undef;
my %expected_default_recipe = ( sv_vcfparser => \&analysis_mip_vcfparser_sv_wgs, );

## Update which recipe to use depending on consensus analysis type
set_recipe_on_analysis_type(
    {
        consensus_analysis_type => $consensus_analysis_type,
        analysis_recipe_href    => \%analysis_recipe,
    }
);

## Then wgs recipe should be set as a default
is(
    $analysis_recipe{sv_vcfparser},
    $expected_default_recipe{sv_vcfparser},
    q{Set default recipe}
);

done_testing();
