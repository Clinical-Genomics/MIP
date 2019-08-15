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
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::File::Format::Feature_file} => [qw{ read_feature_file }],
        q{MIP::Test::Fixtures}             => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Feature_file qw{ read_feature_file };

diag(   q{Test read_feature_file from Feature_file.pm v}
      . $MIP::File::Format::Feature_file::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $MATCHING_COL_NR => 3;
my $log = test_log( { no_screen => 1, } );

## Given
my $feature_file_path =
  catfile( $Bin, qw{ data 643594-miptest aggregated_gene_panel_test.txt } );
my $padding = 1;
my %select_data;
my @select_feature_annotation_columns;
my %tree;

my $is_ok = read_feature_file(
    {
        feature_columns_ref     => \@select_feature_annotation_columns,
        feature_data_href       => \%select_data,
        feature_file_path       => $feature_file_path,
        log                     => $log,
        padding                 => $padding,
        feature_file_type       => q{select_feature},
        feature_matching_column => $MATCHING_COL_NR,
        tree_href               => \%tree,
    }
);

## Then
ok( $is_ok, q{Read feature file} );

done_testing();
