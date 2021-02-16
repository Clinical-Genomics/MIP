#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
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
        q{MIP::Store}          => [qw{ store_metrics }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Store qw{ store_metrics };

diag(   q{Test store_metrics from Store.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( { no_screen => 1, } );

my $test_dir             = File::Temp->newdir();
my $metrics_outfile_path = catfile( $test_dir, q{case_metrics_deliverables.yaml} );

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
my $sample_id   = q{ADM1059A1};
my $sample_id_2 = q{ADM1059A2};

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

## Given sample ids in sample_info
my %sample_info = test_mip_hashes( { mip_hash_name => q{qc_sample_info} } );

## When getting metrics from hash
my @metrics = store_metrics(
    {
        qc_data_href               => \%qc_data,
        sample_info_href           => \%sample_info,
        store_metrics_outfile_path => $metrics_outfile_path,
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
    {
        header => $header,
        id     => $sample_id_2,
        input  => $input,
        name   => $metric_name,
        step   => $recipe_name,
        value  => $metric_value,
    },
);

## Then metric info should be returned
is_deeply( scalar @metrics, 2, q{Got metrics} );

done_testing();
