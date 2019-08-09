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
use MIP::Constants qw{ $COMMA $DOT $EMPTY_STR $SPACE $TAB $UNDERSCORE };
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
        q{MIP::File::Format::Feature_file} =>
          [qw{ parse_feature_file_data parse_feature_file_header tree_annotations }],
        q{MIP::Vcfparser}      => [qw{ define_select_data_headers }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Feature_file
  qw{ parse_feature_file_data parse_feature_file_header tree_annotations };
use MIP::Vcfparser qw{ define_select_data_headers };

diag(   q{Test tree_annotations from Feature_file.pm v}
      . $MIP::File::Format::Feature_file::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $ADK_CHROM_NR    => 10;
Readonly my $ADK_GENE_START  => 1230;
Readonly my $ADK_GENE_STOP   => 1240;
Readonly my $HGNC_ID_COL_NR  => 3;
Readonly my $ADK_HGNC_ID_NR  => 257;
Readonly my $HGNC_SYMBOL_NR  => 4;
Readonly my $MATCHING_COLUMN => 3;
Readonly my $PAX_CHROM_NR    => 21;
Readonly my $PAX_GENE_START  => 21_705_659;
Readonly my $PAX_GENE_STOP   => 21_718_486;
Readonly my $PAX_HGNC_ID_NR  => 8615;
Readonly my $ADK_VAR_START   => 1235;

## Feature data and headers
my @feature_columns         = ( $HGNC_ID_COL_NR, $HGNC_SYMBOL_NR );
my $feature_file_path       = q{a_select_file_path};
my $feature_file_type       = q{select_file};
my $feature_matching_column = $MATCHING_COLUMN;
my @gene_adk_features =
  ( $ADK_CHROM_NR, $ADK_GENE_START, $ADK_GENE_STOP, $ADK_HGNC_ID_NR, q{ADK} );
my $gene_adk_line = join $TAB, @gene_adk_features;
my @gene_pax_features =
  ( $PAX_CHROM_NR, $PAX_GENE_START, $PAX_GENE_STOP, $PAX_HGNC_ID_NR, q{PAX1} );
my $gene_pax_line = join $TAB, @gene_pax_features;
my @gene_no_hgnc_id =
  ( $ADK_CHROM_NR, $ADK_GENE_START, $ADK_GENE_STOP, q{no_hgnc_id}, q{non_coding_region} );
my $gene_no_hgnc_id_line = join $TAB, @gene_no_hgnc_id;
my @headers     = ( q{#chromosome}, qw{ gene_start gene_stop hgnc_id hgnc_symbol } );
my $header_line = join $TAB, @headers;
my $padding     = 1;
my %select_data = define_select_data_headers();
my %vcf_record;
my %tree;

parse_feature_file_header(
    {
        feature_columns_ref => \@feature_columns,
        feature_data_href   => \%select_data,
        feature_file_type   => $feature_file_type,
        feature_file_path   => $feature_file_path,
        header_line         => $header_line,
    }
);
GENE:

foreach my $gene_line ( $gene_adk_line, $gene_pax_line, $gene_no_hgnc_id_line ) {

## Parse feature data and build interval tree annotations
    parse_feature_file_data(
        {
            data_line               => $gene_line,
            feature_columns_ref     => \@feature_columns,
            feature_data_href       => \%select_data,
            feature_file_type       => $feature_file_type,
            feature_matching_column => $feature_matching_column,
            padding                 => $padding,
            tree_href               => \%tree,
        }
    );

}
## Given a snv when no matching tree annotation
my @no_tree_variant_ann = ( 1, 1, $DOT, qw{ A G } );
tree_annotations(
    {
        data_href         => \%select_data,
        line_elements_ref => \@no_tree_variant_ann,
        feature_file_type => $feature_file_type,
        record_href       => \%vcf_record,
        tree_href         => \%tree,
    }
);
is( keys %vcf_record, 0, q{Did not set feature annotation to vcf record} );

## Given a snv variant line which matches feature line
my @snv_single_variant_elements = ( $ADK_CHROM_NR, $ADK_VAR_START, $DOT, qw{ A G } );
my %noid_region                 = tree_annotations(
    {
        data_href         => \%select_data,
        line_elements_ref => \@snv_single_variant_elements,
        feature_file_type => $feature_file_type,
        record_href       => \%vcf_record,
        tree_href         => \%tree,
    }
);

my %expected_vcf_record = (
    q{INFO_addition_}
      . $feature_file_type => {
        hgnc_id     => qq{$ADK_HGNC_ID_NR,no_hgnc_id},
        hgnc_symbol => q{ADK,non_coding_region},
      }
);

my $expected_variant_id = join $UNDERSCORE, ( $ADK_CHROM_NR, $ADK_VAR_START, qw{ A G } );
## Then set ADK annotation features to vcf hash
is_deeply( \%vcf_record, \%expected_vcf_record,
    q{Set ADK annotation features to vcf record hash} );

## Then no_id_rgion should have a variant_key_id
ok( $noid_region{$expected_variant_id}, q{Set variant id for no id region} );

done_testing();
