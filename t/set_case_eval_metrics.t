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
use MIP::Test::Fixtures qw{ test_log test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Qccollect}      => [qw{ set_case_eval_metrics }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qccollect qw{ set_case_eval_metrics };

diag(   q{Test set_case_eval_metrics from Qccollect.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

Readonly my $PCT_PF_READS_ALIGNED => 0.95;
Readonly my $PCT_ADAPTER          => 0.0005;

my $log = test_log( { no_screen => 1, } );

my %sample_info = test_mip_hashes(
    {
        mip_hash_name => q{qc_sample_info},
    }
);

## Remove RNA sample from test to get a wgs consensus type
delete $sample_info{sample}{ADM1059A3};
delete $sample_info{analysis_type}{ADM1059A3};

## Given a hash of evaluation metrics
my %eval_metric = (
    wgs => {
        collectmultiplemetrics => {
            PCT_PF_READS_ALIGNED => {
                lt => $PCT_PF_READS_ALIGNED,
            },
            PCT_ADAPTER => {
                gt => $PCT_ADAPTER,
            },
        },
    },
);

## Then set the relevant evaluation metrics for the analysis
my %analysis_eval_metric;
set_case_eval_metrics(
    {
        analysis_eval_metric_href => \%analysis_eval_metric,
        eval_metric_href          => \%eval_metric,
        sample_info_href          => \%sample_info,
    }
);

my %expected = (
    collectmultiplemetrics => {
        PCT_PF_READS_ALIGNED => {
            lt => $PCT_PF_READS_ALIGNED,
        },
        PCT_ADAPTER => {
            gt => $PCT_ADAPTER,
        },
    },
);
is_deeply( \%analysis_eval_metric, \%expected, q{Set case eval metrics} );

done_testing();
