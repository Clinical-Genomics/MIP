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
our $VERSION = 1.02;

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

Readonly my $PCT_PF_READS_ALIGNED => 0.9;
Readonly my $PCT_ADAPTER          => 0.0005;

## Given qc_data with metrics from an analysis
my $metric_gt = q{PCT_ADAPTER};
my $metric_lt = q{PCT_PF_READS_ALIGNED};
my $recipe    = q{collectmultiplemetrics};
my %qc_data   = (
    recipe => {
        $recipe => { $metric_lt => $PCT_PF_READS_ALIGNED },
        $recipe => { $metric_gt => $PCT_ADAPTER },
    },
);

## Given evaluation thresholds
my %evaluate_metric = (
    collectmultiplemetrics => {
        PCT_PF_READS_ALIGNED => {
            lt => $PCT_PF_READS_ALIGNED,
        },
        PCT_ADAPTER => {
            gt => $PCT_ADAPTER,
        },
    },
);

## Given a the recipe in qc_sample_info
test_mip_hashes(
    {
        mip_hash_name => q{qc_sample_info},
        recipe_name   => $recipe,
    }
);

## When evaluating the recipe and metrics
my $is_ok = evaluate_case_qc_parameters(
    {
        evaluate_metric_href => \%evaluate_metric,
        qc_data_href         => \%qc_data,
    }
);

## Then metrics should be evaluated and return true
ok( $is_ok, q{Evaluated case qc parameter} );

done_testing();
