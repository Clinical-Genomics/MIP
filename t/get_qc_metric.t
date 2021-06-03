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

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Qc_data} => [qw{ get_qc_metric }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qc_data qw{ get_qc_metric };

diag(   q{Test get_qc_metric from Qc_data.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a header
my $header = q{a_header};

## Given an input source
my $input = q{an_infile};

## Given a metric
my $metric_name  = q{AT_DROPOUT};
my $metric_value = 1;

## Given a recipe
my $recipe_name = q{picardtools_collectmultiplemetrics};

## Given two sample_ids
my $sample_id   = q{sample_1};
my $sample_id_2 = q{sample_2};

## Given qc_data
my %qc_data = (
    metrics => [
        {
            header => $header,
            id     => $sample_id,
            input  => $input,
            name   => $metric_name,
            step   => $recipe_name,
            value  => $metric_value,
        },
        {
            header => $header,
            id     => $sample_id_2,
            input  => $input,
            name   => $metric_name,
            step   => $recipe_name,
            value  => $metric_value,
        },
    ],
);

## When getting metrics from hash
my @metrics = get_qc_metric(
    {
        header       => $header,
        id           => $sample_id,
        input        => $input,
        metric_name  => $metric_name,
        metric_value => $metric_value,
        qc_data_href => \%qc_data,
        recipe_name  => $recipe_name,
    }
);

my @expected_metrics = (
    {
        header => $header,
        id     => $sample_id,
        input  => $input,
        name   => $metric_name,
        step   => $recipe_name,
        value  => $metric_value,
    },
);

## Then metric info should be returned
is_deeply( \@metrics, \@expected_metrics, q{Got metrics} );

## Given a non matching sample_id
my $not_matching_sample_id = q{sample_is_not_present};

## When getting metric from hash
my @non_matching_metrics = get_qc_metric(
    {
        header       => $header,
        id           => $not_matching_sample_id,
        input        => $input,
        metric_name  => $metric_name,
        metric_value => $metric_value,
        qc_data_href => \%qc_data,
        recipe_name  => $recipe_name,
    }
);

is( scalar @non_matching_metrics, 0, q{Return zero for non matching metric} );

done_testing();
