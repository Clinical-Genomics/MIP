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
    my %perl_module = ( q{MIP::Qc_data} => [qw{ set_metrics_to_store }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qc_data qw{ set_metrics_to_store };

diag(   q{Test set_metrics_to_store from Qc_data.pm}
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

## Given qc_data
my %qc_data;

## Given a recipe
my $recipe_name = q{picardtools_collectmultiplemetrics};

## Given a sample_id
my $sample_id = q{sample_1};

## When setting metric to store
set_metrics_to_store(
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

my %expected_metric_info = (
    metrics => [
        {
            header => $header,
            id     => $sample_id,
            input  => $input,
            name   => $metric_name,
            step   => $recipe_name,
            value  => $metric_value,
        },
    ],
);

## Then metric info should be set
is_deeply( \%qc_data, \%expected_metric_info, q{Set metrics in hash} );

## Given a duplicate metric

## When setting metric to store
set_metrics_to_store(
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

## Then metric info should be set
is_deeply( \%qc_data, \%expected_metric_info, q{Skip duplicate metric} );

done_testing();
