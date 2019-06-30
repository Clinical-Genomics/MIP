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
use Modern::Perl qw{ 2017 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE $TAB };
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
        q{MIP::Vcfparser} => [qw{ define_select_data_headers parse_feature_file_header }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfparser qw{ define_select_data_headers parse_feature_file_header };

diag(   q{Test parse_feature_file_header from Vcfparser.pm v}
      . $MIP::Vcfparser::VERSION
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
my $feature_file_key  = q{select_file};
my $feature_file_path = q{a_select_file_path};
my $header_key        = q{Not present in feature data};
my @headers     = ( q{#chromosome}, qw{ gene_start gene_stop hgnc_id hgnc_symbol } );
my $header_line = join $TAB, @headers;

#say STDERR $header_line;
my $is_ok = parse_feature_file_header(
    {
        feature_columns_ref => \@feature_columns,
        feature_data_href   => \%feature_data,
        feature_file_key    => $feature_file_key,
        feature_file_path   => $feature_file_path,
        header_line         => $header_line,
    }
);

## Then
ok( $is_ok, q{Parsed feature file header } );

done_testing();
