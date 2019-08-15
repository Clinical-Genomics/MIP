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
        q{MIP::Qc_data}        => [qw{ parse_regexp_hash_and_collect }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qc_data qw{ parse_regexp_hash_and_collect };

diag(   q{Test parse_regexp_hash_and_collect from Qc_data.pm v}
      . $MIP::Qc_data::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a valid regexp when file exists and data is without headers
my $outdirectory = dirname($Bin);
my $outfile      = q{mip};
my %qc_recipe_data;
my %qc_header;
my $recipe_name = q{mip};
my $regexp = q{perl -nae 'if ($_ =~ /\bour\s\$VERSION\b/xms) {print q{Got version};}'};
my $regexp_key = q{test_add_data_from_regexp};
my %regexp     = ( $recipe_name => { $regexp_key => $regexp, }, );

parse_regexp_hash_and_collect(
    {
        outdirectory        => $outdirectory,
        outfile             => $outfile,
        qc_header_href      => \%qc_header,
        qc_recipe_data_href => \%qc_recipe_data,
        recipe              => $recipe_name,
        regexp_href         => \%regexp,
    }
);

my %expected_qc_data = ( $recipe_name => { $regexp_key => [ qw{Got version}, ], } );

## Then reg exp returned data should have been added to qc_data
is_deeply( \%qc_recipe_data, \%expected_qc_data,
    q{Parsed regexp hash and collected data without header} );

## Given a regexp when file exists and data comes with header
$regexp{$recipe_name}{header} = q{test};    # Should not return anything

my $is_ok = parse_regexp_hash_and_collect(
    {
        outdirectory        => $outdirectory,
        outfile             => $outfile,
        qc_header_href      => \%qc_header,
        qc_recipe_data_href => \%qc_recipe_data,
        recipe              => $recipe_name,
        regexp_href         => \%regexp,
    }
);
ok( $is_ok, q{Parsed regexp hash and collected data with header} );

done_testing();
