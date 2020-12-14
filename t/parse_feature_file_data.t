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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::File::Format::Feature_file} => [qw{ parse_feature_file_data }],
        q{MIP::Vcfparser}                  => [qw{ define_select_data_headers }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Feature_file qw{ parse_feature_file_data };
use MIP::Vcfparser qw{ define_select_data_headers };

diag(   q{Test parse_feature_file_data from Feature_file.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $GENE_START                 => 1234;
Readonly my $GENE_STOP                  => 1235;
Readonly my $HGNC_ID_NR                 => 3;
Readonly my $HGNC_SYMBOL_NR             => 4;
Readonly my $MATCHING_COLUMN            => 3;
Readonly my $MATCHING_COLUMN_WITH_SPACE => 4;

## Given feature data line
my @data_features           = ( 1, $GENE_START, $GENE_STOP, $HGNC_ID_NR, q{a gene} );
my $data_line               = join $TAB, @data_features;
my @feature_columns         = ( $HGNC_ID_NR, $HGNC_SYMBOL_NR );
my %feature_data            = define_select_data_headers();
my $feature_file_type       = q{select_file};
my $padding                 = 1;
my $feature_matching_column = $MATCHING_COLUMN;
my %tree;

my $is_ok = parse_feature_file_data(
    {
        data_line               => $data_line,
        feature_columns_ref     => \@feature_columns,
        feature_data_href       => \%feature_data,
        feature_file_type       => $feature_file_type,
        feature_matching_column => $feature_matching_column,
        padding                 => $padding,
        tree_href               => \%tree,
    }
);
my %expected_feature_data = ( $HGNC_ID_NR => $HGNC_ID_NR );

## Then return true if parsed
ok( $is_ok, q{Parsed feature file data } );
is( $feature_data{$HGNC_ID_NR}, $expected_feature_data{$HGNC_ID_NR},
    q{Set feature data} );

## Given a data with whitespace
$feature_matching_column = $MATCHING_COLUMN_WITH_SPACE;
parse_feature_file_data(
    {
        data_line               => $data_line,
        feature_columns_ref     => \@feature_columns,
        feature_data_href       => \%feature_data,
        feature_file_type       => $feature_file_type,
        feature_matching_column => $feature_matching_column,
        padding                 => $padding,
        tree_href               => \%tree,
    }
);
$expected_feature_data{a_gene} = q{a_gene};

## Then underscore should replace whitespace in feature data
is( $feature_data{a_gene}, $expected_feature_data{a_gene}, q{Set feature data} );

done_testing();
