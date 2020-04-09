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
        q{MIP::Qccollect}      => [qw{ evaluate_case_qc_parameters }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qccollect qw{ evaluate_case_qc_parameters };

diag(   q{Test evaluate_case_qc_parameters from Qccollect.pm v}
      . $MIP::Qccollect::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

Readonly my $FRACTION_OF_COMMON_VARIANTS      => 0.3;
Readonly my $FRACTION_OF_COMMON_VARIANTS_EVAL => 0.55;

## Given
my $metric_lt = q{fraction_of_common_variants};
my $recipe_lt = q{variant_integrity_ar_father};
my %qc_data =
  ( recipe => { $recipe_lt => { $metric_lt => $FRACTION_OF_COMMON_VARIANTS }, }, );

my %evaluate_metric = (
    variant_integrity_ar_father => {
        fraction_of_common_variants => {
            lt => $FRACTION_OF_COMMON_VARIANTS,
        },
    },
);

my %sample_info = test_mip_hashes(
    {
        mip_hash_name => q{qc_sample_info},
        recipe_name   => $recipe_lt,
    }
);

my $is_ok = evaluate_case_qc_parameters(
    {
        evaluate_metric_href => \%evaluate_metric,
        qc_data_href         => \%qc_data,
    }
);

## Then
ok( $is_ok, q{Evaluated case qc parameter} );

done_testing();
