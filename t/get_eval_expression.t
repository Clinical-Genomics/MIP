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
        q{MIP::Qccollect}      => [qw{ get_eval_expression }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qccollect qw{ get_eval_expression };

diag(   q{Test get_eval_expression from Qccollect.pm v}
      . $MIP::Qccollect::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

Readonly my $PCT_PF_READS_ALIGNED            => 0.95;
Readonly my $PCT_ADAPTER                     => 0.0005;
Readonly my $VARIANT_INTEGRITY_AR_MENDEL_WGS => 0.06;

my %analysis_eval_metric = (
    ADM1059A1 => {
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
            gt => $VARIANT_INTEGRITY_AR_MENDEL_WGS,
        },
    },
);

## Given a sample id and recipe outfile key
## Then get the relevant expression
my %eval_expression = get_eval_expression(
    {
        eval_metric_href => \%analysis_eval_metric,
        recipe           => q{collectmultiplemetrics},
        sample_id        => q{ADM1059A1},
    }
);
my %expected = (
    PCT_PF_READS_ALIGNED => {
        lt => $PCT_PF_READS_ALIGNED,
    },
    PCT_ADAPTER => {
        gt => $PCT_ADAPTER,
    },
);

is_deeply( \%eval_expression, \%expected, q{Get sample eval expression} );

## Given no sample_id and recipe outfile key
## Then get the relevant expression
%eval_expression = get_eval_expression(
    {
        eval_metric_href => \%analysis_eval_metric,
        recipe           => q{variant_integrity_ar_mendel},
    }
);
%expected = (
    fraction_of_errors => {
        gt => $VARIANT_INTEGRITY_AR_MENDEL_WGS,
    },
);
is_deeply( \%eval_expression, \%expected, q{Get case eval expression} );

done_testing();
