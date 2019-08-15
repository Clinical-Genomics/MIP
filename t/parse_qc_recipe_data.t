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
        q{MIP::Qc_data}        => [qw{ parse_qc_recipe_data }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qc_data qw{ parse_qc_recipe_data };

diag(   q{Test parse_qc_recipe_data from Qc_data.pm v}
      . $MIP::Qc_data::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $AT_DROPOUT => 1.721239;

## Given data when header key is "data"
my $infile = q{an_infile};
my %qc_data;
my %qc_header;
my $recipe_name    = q{collecthsmetrics};
my $regexp_key     = q{data};
my $sample_id      = q{sample_1};
my $value          = $AT_DROPOUT;
my %qc_recipe_data = ( $recipe_name => { $regexp_key => [$value], }, );
my %regexp         = ( $recipe_name => { $regexp_key => [$value], } );
my %sample_info;

my $is_ok = parse_qc_recipe_data(
    {
        infile              => $infile,
        qc_data_href        => \%qc_data,
        qc_header_href      => \%qc_header,
        qc_recipe_data_href => \%qc_recipe_data,
        recipe              => $recipe_name,
        regexp_href         => \%regexp,
        sample_id           => $sample_id,
        sample_info_href    => \%sample_info,
    }
);

## Then sub should return true
ok( $is_ok, q{Parsed qc recipe data table metrics with no header} );

## Given data when header key is "header"
$regexp_key = q{header};
%regexp     = ( $recipe_name => { $regexp_key => [$value], } );

$is_ok = parse_qc_recipe_data(
    {
        infile              => $infile,
        qc_data_href        => \%qc_data,
        qc_header_href      => \%qc_header,
        qc_recipe_data_href => \%qc_recipe_data,
        recipe              => $recipe_name,
        regexp_href         => \%regexp,
        sample_id           => $sample_id,
        sample_info_href    => \%sample_info,
    }
);

## Then sub should return true
ok( $is_ok, q{Parsed qc recipe data table metrics with header} );

done_testing();
