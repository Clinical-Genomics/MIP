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
        q{MIP::Qc_data}        => [qw{ parse_qc_recipe_table_data }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qc_data qw{ parse_qc_recipe_table_data };

diag(   q{Test parse_qc_recipe_table_data from Qc_data.pm v}
      . $MIP::Qc_data::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $AT_DROPOUT => 1.721239;

## Given a data table when sample and infile
my $regexp_header_key = q{header};
my $regexp_key        = q{data};
my $infile            = q{an_infile};
my %qc_data;
my $key            = q{AT_DROPOUT};
my $recipe_name    = q{collecthsmetrics};
my $sample_id      = q{sample_1};
my $value          = $AT_DROPOUT;
my %qc_header      = ( $recipe_name => { $regexp_header_key => [$key], } );
my %qc_recipe_data = ( $recipe_name => { $regexp_key => [$value] } );
my %regexp         = ( $recipe_name => { $regexp_key => q{a_regexp} } );

my $is_ok = parse_qc_recipe_table_data(
    {
        infile              => $infile,
        qc_data_href        => \%qc_data,
        qc_header_href      => \%qc_header,
        qc_recipe_data_href => \%qc_recipe_data,
        regexp_href         => \%regexp,
        recipe              => $recipe_name,
        sample_id           => $sample_id,
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

## Then sub should return true
ok( $is_ok, q{Parsed qc recipe data table metrics} );

## Then data table should be added to qc data
is_deeply( \%qc_data, \%expected_qc_data,
    q{Set table metrics in qc_data on sample and infile level } );

done_testing();
