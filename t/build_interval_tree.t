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
use Set::IntervalTree;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SEMICOLON $SPACE };
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
        q{MIP::Vcfparser}      => [qw{ build_interval_tree }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfparser qw{ build_interval_tree };

diag(   q{Test build_interval_tree from Vcfparser.pm v}
      . $MIP::Vcfparser::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $FEATURE_COL_START      => 3;
Readonly my $FEATURE_COL_STOP       => 4;
Readonly my $GENE_START             => 1234;
Readonly my $GENE_STOP              => 1236;
Readonly my $HGNC_SYMBOL_COL_NR     => 3;
Readonly my $HGNC_ID_COL_NR         => 4;
Readonly my $SNV_START              => 1234;
Readonly my $SNV_STOP               => 1235;
Readonly my $SNV_OUTSIDE_TREE_START => 123_40;
Readonly my $SNV_OUTSIDE_TREE_STOP  => 123_41;
Readonly my $WITHESPACE_COL_NR      => 5;

## Given a line with features when overlapping
my $contig            = 1;
my @feature_columns   = ( $HGNC_SYMBOL_COL_NR, $HGNC_ID_COL_NR );
my @line_elements     = ( $contig, $GENE_START, $GENE_STOP, q{hgnc_symbol}, q{hgnc_id} );
my $padding           = 2;
my $feature_file_type = q{select_feature};
my %tree;

build_interval_tree(
    {
        feature_columns_ref => \@feature_columns,
        line_elements_ref   => \@line_elements,
        padding             => $padding,
        feature_file_type   => $feature_file_type,
        tree_href           => \%tree,
    }
);

my @expected_features = join $SEMICOLON,
  @line_elements[ $FEATURE_COL_START .. $FEATURE_COL_STOP ];
my $features_ref = $tree{$feature_file_type}{$contig}->fetch( $SNV_START, $SNV_STOP );

## Then keys for where to place the interval tree should be set
ok( exists $tree{$feature_file_type}{$contig}, q{Set keys for interval tree} );

## Then interval tree should return features array ref
is_deeply( \@{$features_ref}, \@expected_features, q{Set interval tree} );

## Given a feature when not overlapping
$features_ref = $tree{$feature_file_type}{$contig}
  ->fetch( $SNV_OUTSIDE_TREE_START, $SNV_OUTSIDE_TREE_STOP );

## Then no feature should be returned
is( @{$features_ref}, 0, q{No feature from interval tree} );

## Given feature with whitespace
push @line_elements,   q{w ithe space};
push @feature_columns, $WITHESPACE_COL_NR;

## New contig
$contig++;
$line_elements[0] = $contig;

build_interval_tree(
    {
        feature_columns_ref => \@feature_columns,
        line_elements_ref   => \@line_elements,
        padding             => $padding,
        feature_file_type   => $feature_file_type,
        tree_href           => \%tree,
    }
);

@expected_features = join $SEMICOLON,
  ( @line_elements[ $FEATURE_COL_START .. $FEATURE_COL_STOP ], q{w_ithe_space} );
$features_ref = $tree{$feature_file_type}{$contig}->fetch( $SNV_START, $SNV_STOP );

## Then interval tree should return features array ref
is_deeply( \@{$features_ref}, \@expected_features,
    q{Translated whitespace in feature string} );

done_testing();
