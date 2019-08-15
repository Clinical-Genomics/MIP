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
use MIP::Constants qw{ $COLON $COMMA $SPACE $UNDERSCORE};
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
        q{MIP::Qccollect}      => [qw{ check_qc_metric define_evaluate_metric }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qccollect qw{ check_qc_metric define_evaluate_metric };

diag(   q{Test check_qc_metric from Qccollect.pm v}
      . $MIP::Qccollect::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $PCT_ADAPTER_PASS             => 0.0001;
Readonly my $PCT_ADAPTER_FAIL             => 0.0010;
Readonly my $PERCENTAGE_MAPPED_READS_PASS => 99;
Readonly my $PERCENTAGE_MAPPED_READS_FAIL => 90;

## Given sample id and metric to evaluate when requirement lesser than fulfilled
my $sample_id = q{ADM1059A1};
my $metric_lt = q{percentage_mapped_reads};
my %qc_data;
my $recipe_lt = q{bamstats};

my %sample_info = test_mip_hashes(
    {
        mip_hash_name => q{qc_sample_info},
        recipe_name   => $recipe_lt,
    }
);

## Defines recipes, metrics and thresholds to evaluate
my %evaluate_metric = define_evaluate_metric(
    {
        sample_info_href => \%sample_info,
    }
);

check_qc_metric(
    {
        metric                => $metric_lt,
        qc_data_href          => \%qc_data,
        qc_metric_value       => $PERCENTAGE_MAPPED_READS_PASS,
        recipe                => $recipe_lt,
        reference_metric_href => $evaluate_metric{$sample_id}{$recipe_lt}{$metric_lt},
    }
);

## Then qc data should be undef if evaluation passed
is( $qc_data{evaluation}{$recipe_lt},
    undef, q{Passed evaluation for sample and metric operation lt } );

## Given sample id and metric to evaluate when requirement lesser than NOT fulfilled
check_qc_metric(
    {
        metric                => $metric_lt,
        qc_data_href          => \%qc_data,
        qc_metric_value       => $PERCENTAGE_MAPPED_READS_FAIL,
        recipe                => $recipe_lt,
        reference_metric_href => $evaluate_metric{$sample_id}{$recipe_lt}{$metric_lt},
    }
);

my $status =
    q{FAILED:}
  . $recipe_lt
  . $UNDERSCORE
  . $metric_lt
  . $COLON
  . $PERCENTAGE_MAPPED_READS_FAIL;
my %expected_evaluation = ( evaluation => { $recipe_lt => [$status], }, );

## Then qc data should have a failed entry status
is_deeply( \%qc_data, \%expected_evaluation,
    q{Failed evaluation for sample and metric operation lt } );

## Given sample id and metric to evaluate when requirement greater than fulfilled
my $metric_gt = q{PCT_ADAPTER};
my $recipe_gt = q{collectmultiplemetrics};

check_qc_metric(
    {
        metric                => $metric_gt,
        qc_data_href          => \%qc_data,
        qc_metric_value       => $PCT_ADAPTER_PASS,
        recipe                => $recipe_gt,
        reference_metric_href => $evaluate_metric{$sample_id}{$recipe_gt}{$metric_gt},
    }
);

## Then qc data should be undef if evaluation passed
is( $qc_data{evaluation}{$recipe_gt},
    undef, q{Passed evaluation for sample and metric operation gt } );

## Given sample id and metric to evaluate when requirement greater than NOT fulfilled
check_qc_metric(
    {
        metric                => $metric_gt,
        qc_data_href          => \%qc_data,
        qc_metric_value       => $PCT_ADAPTER_FAIL,
        recipe                => $recipe_gt,
        reference_metric_href => $evaluate_metric{$sample_id}{$recipe_gt}{$metric_gt},
    }
);
$status = q{FAILED:} . $recipe_gt . $UNDERSCORE . $metric_gt . $COLON . $PCT_ADAPTER_FAIL;
$expected_evaluation{evaluation}{$recipe_gt} = [$status];

## Then qc data should have a failed entry status
is_deeply( \%qc_data, \%expected_evaluation,
    q{Failed evaluation for sample and metric operation gt } );

done_testing();
