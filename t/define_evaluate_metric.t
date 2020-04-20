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
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

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
        q{MIP::Qccollect}      => [qw{ define_evaluate_metric }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qccollect qw{ define_evaluate_metric };

diag(   q{Test define_evaluate_metric from Qccollect.pm v}
      . $MIP::Qccollect::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

Readonly my $PCT_PF_READS_ALIGNED        => 0.95;
Readonly my $PCT_ADAPTER                 => 0.0005;
Readonly my $VARIANT_INTEGRITY_AR_MENDEL => 0.06;

my $log = test_log( { no_screen => 1, } );

my %sample_info = test_mip_hashes(
    {
        mip_hash_name => q{qc_sample_info},
    }
);

## Given a file with evaluation metrics
my $eval_metric_file =
  catfile( dirname($Bin), qw{ t data references qc_eval_metric_-v1.0-.yaml} );

## Then set the relevant evaluation metrics for the analysis
my %evaluate_metric = define_evaluate_metric(
    {
        eval_metric_file => $eval_metric_file,
        sample_info_href => \%sample_info,
    }
);

my %expected = (
    ADM1059A1 => {
        collecthsmetrics => {
            MEAN_TARGET_COVERAGE => {
                lt => 150,
            },
        },
        collectmultiplemetrics => {
            PCT_PF_READS_ALIGNED => {
                lt => $PCT_PF_READS_ALIGNED,
            },
            PCT_ADAPTER => {
                gt => $PCT_ADAPTER,
            },
        },
    },
    ADM1059A2 => {
        collecthsmetrics => {
            MEAN_TARGET_COVERAGE => {
                lt => 150,
            },
        },
        collectmultiplemetrics => {
            PCT_PF_READS_ALIGNED => {
                lt => $PCT_PF_READS_ALIGNED,
            },
            PCT_ADAPTER => {
                gt => $PCT_ADAPTER,
            },
        },
    },
    variant_integrity_ar_mendel => {
        fraction_of_errors => {
            gt => $VARIANT_INTEGRITY_AR_MENDEL,
        },
    },
    ADM1059A3 => {
        collecthsmetrics => {
            MEAN_TARGET_COVERAGE => {
                lt => 150,
            },
        },
    },
);
is_deeply( \%evaluate_metric, \%expected, q{Define analysis eval metrics} );

done_testing();