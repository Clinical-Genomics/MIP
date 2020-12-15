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
use MIP::Test::Fixtures qw{ test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Qccollect}      => [qw{ parse_sample_qc_metric }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qccollect qw{ parse_sample_qc_metric };

diag(   q{Test parse_sample_qc_metric from Qccollect.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $PCT_TARGET_BASES_10X         => 0.95;
Readonly my $PERCENTAGE_MAPPED_READS_EVAL => 95;
Readonly my $PERCENTAGE_MAPPED_READS_PASS => 99;

## Given sample qc data when metric lacks a header
my $infile    = q{an_infile};
my $metric    = q{percentage_mapped_reads};
my $header_recipe = q{collecthsmetrics};
my $header_metric = q{PCT_TARGET_BASES_10X};
my $recipe    = q{bamstats};
my $sample_id = q{ADM1059A1};
my %qc_data   = (
    sample => {
        $sample_id =>
          { $infile => { $recipe => { $metric => $PERCENTAGE_MAPPED_READS_PASS }, }, },
    }
);

## Alias
my $qc_data_recipe_href =
  \%{ $qc_data{sample}{$sample_id}{$infile}{$recipe} };

my %sample_info = test_mip_hashes(
    {
        mip_hash_name => q{qc_sample_info},
        recipe_name   => $recipe,
    }
);

## Defines recipes, metrics and thresholds to evaluate
my %evaluate_metric = (
    $sample_id => {
        $recipe => {
            $metric => {
                lt => $PERCENTAGE_MAPPED_READS_EVAL,
            },
        },
        $header_recipe => {
            $header_metric => {
                lt => $PCT_TARGET_BASES_10X,
            },
        },
    },
);

my $is_ok = parse_sample_qc_metric(
    {
        evaluate_metric_href => \%evaluate_metric,
        qc_data_href         => \%qc_data,
        qc_data_recipe_href  => $qc_data_recipe_href,
        recipe_name          => $recipe,
        sample_id            => $sample_id,
    }
);
ok( $is_ok, q{Parsed sample qc data metrics} );

## Given sample qc data when metric has a header
my $data_header   = q{first_of_pair};
$qc_data{sample}{$sample_id}{$infile}{$header_recipe}{header}{$data_header}{$header_metric} = $PCT_TARGET_BASES_10X;

## Alias
my $qc_data_header_recipe_href =
  \%{ $qc_data{sample}{$sample_id}{$infile}{$header_recipe} };

$is_ok = parse_sample_qc_metric(
    {
        evaluate_metric_href => \%evaluate_metric,
        qc_data_href         => \%qc_data,
        qc_data_recipe_href  => $qc_data_header_recipe_href,
        recipe_name          => $header_recipe,
        sample_id            => $sample_id,
    }
);
ok( $is_ok, q{Parsed sample qc data header metrics} );

done_testing();
