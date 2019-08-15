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
use MIP::Constants qw{ $COMMA $SPACE $TAB };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

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
        q{MIP::File::Format::Feature_file} => [qw{ parse_feature_file_header }],
        q{MIP::Vcfparser}                  => [qw{ define_select_data_headers }],
        q{MIP::Test::Fixtures}             => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Feature_file qw{ parse_feature_file_header };
use MIP::Vcfparser qw{ define_select_data_headers };

diag(   q{Test parse_feature_file_header from Feature_file.pm v}
      . $MIP::File::Format::Feature_file::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $HGNC_ID_NR     => 3;
Readonly my $HGNC_SYMBOL_NR => 4;

## Given
my @feature_columns   = ( $HGNC_ID_NR, $HGNC_SYMBOL_NR );
my %feature_data      = define_select_data_headers();
my $feature_file_type = q{select_file};
my $feature_file_path = q{a_select_file_path};
my @headers     = ( q{#chromosome}, qw{ gene_start gene_stop hgnc_id hgnc_symbol } );
my $header_line = join $TAB, @headers;

my $is_ok = parse_feature_file_header(
    {
        feature_columns_ref => \@feature_columns,
        feature_data_href   => \%feature_data,
        feature_file_type   => $feature_file_type,
        feature_file_path   => $feature_file_path,
        header_line         => $header_line,
    }
);

## Then
ok( $is_ok, q{Parsed feature file header } );

done_testing();
