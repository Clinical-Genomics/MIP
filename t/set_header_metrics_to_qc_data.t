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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Qc_data}        => [qw{ set_header_metrics_to_qc_data }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qc_data qw{ set_header_metrics_to_qc_data };

diag(   q{Test set_header_metrics_to_qc_data from Qc_data.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $AT_DROPOUT             => 1.721239;
Readonly my $BAIT_DESIGN_EFFICIENCY => 0.580692;

## Given a data table when sample and infile
my $regexp_header_key = q{header};
my $regexp_key        = q{data};
my $infile            = q{an_infile};
my %qc_data;
my $key         = q{AT_DROPOUT};
my $recipe_name = q{collecthsmetrics};
my $sample_id   = q{sample_1};
my $value       = $AT_DROPOUT;

set_header_metrics_to_qc_data(
    {
        infile            => $infile,
        key               => $key,
        qc_data_href      => \%qc_data,
        regexp_header_key => $regexp_header_key,
        regexp_key        => $regexp_key,
        recipe_name       => $recipe_name,
        sample_id         => $sample_id,
        value             => $value,
    }
);

my %expected_qc_data = (
    sample => {
        $sample_id => {
            $infile => {
                $recipe_name =>
                  { $regexp_header_key => { $regexp_key => { $key => $value, }, }, },
            },
        },
    },
);

## Then data table should be added to qc data
is_deeply( \%qc_data, \%expected_qc_data,
    q{Set table metrics in qc_data on sample and infile level } );

## Given a data table when only recipe
$infile    = undef;
$sample_id = undef;
my $key_bait_design_efficiency = q{BAIT_DESIGN_EFFICIENCY};
set_header_metrics_to_qc_data(
    {
        infile            => $infile,
        key               => $key_bait_design_efficiency,
        qc_data_href      => \%qc_data,
        regexp_header_key => $regexp_header_key,
        regexp_key        => $regexp_key,
        recipe_name       => $recipe_name,
        sample_id         => $sample_id,
        value             => $BAIT_DESIGN_EFFICIENCY,
    }
);

$expected_qc_data{$recipe_name}{$regexp_header_key}{$regexp_key}
  {$key_bait_design_efficiency} = $BAIT_DESIGN_EFFICIENCY;

## Then data table should be added to qc data
is_deeply( \%qc_data, \%expected_qc_data,
    q{Set table metrics in qc_data on case level } );

done_testing();
