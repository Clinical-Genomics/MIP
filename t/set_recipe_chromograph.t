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
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

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
        q{MIP::Recipes::Analysis::Chromograph} =>
          [qw{ analysis_chromograph analysis_chromograph_proband }],
        q{MIP::Set::Analysis}  => [qw{ set_recipe_chromograph }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::Analysis qw{ set_recipe_chromograph };
use MIP::Recipes::Analysis::Chromograph
  qw{ analysis_chromograph analysis_chromograph_proband };

diag(   q{Test set_recipe_chromograph from Analysis.pm v}
      . $MIP::Set::Analysis::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given
my %analysis_recipe;
my %sample_info    = test_mip_hashes( { mip_hash_name => q{qc_sample_info}, } );
my $proband_id     = q{ADM1059A1};
my $not_proband_id = q{ADM1059A2};

## Given a non trio
$sample_info{has_trio} = 0;

## Set chromograph recipe depending on pedigree and proband
set_recipe_chromograph(
    {
        analysis_recipe_href => \%analysis_recipe,
        sample_id            => $proband_id,
        sample_info_href     => \%sample_info,
    }
);
my %expected_single_sample_recipe = ( chromograph_ar => \&analysis_chromograph, );

## Then use not proband and not trio recipe
is_deeply(
    \%analysis_recipe,
    \%expected_single_sample_recipe,
    q{Set single sample recipe for not trio}
);

## Given a trio when not proband
$sample_info{has_trio} = 1;

## Set chromograph recipe depending on pedigree and proband
set_recipe_chromograph(
    {
        analysis_recipe_href => \%analysis_recipe,
        sample_id            => $not_proband_id,
        sample_info_href     => \%sample_info,
    }
);
my %expected_trio_not_proband_recipe = ( chromograph_ar => \&analysis_chromograph, );

## Then use not proband and not trio recipe
is_deeply(
    \%analysis_recipe,
    \%expected_trio_not_proband_recipe,
    q{Set sample recipe when trio but not proband}
);

## Given a trio when proband

## Set chromograph recipe depending on pedigree and proband
set_recipe_chromograph(
    {
        analysis_recipe_href => \%analysis_recipe,
        sample_id            => $proband_id,
        sample_info_href     => \%sample_info,
    }
);
my %expected_trio_and_proband_recipe =
  ( chromograph_ar => \&analysis_chromograph_proband, );

## Then use proband and trio recipe
is_deeply(
    \%analysis_recipe,
    \%expected_trio_and_proband_recipe,
    q{Set sample recipe when trio and proband}
);

done_testing();
