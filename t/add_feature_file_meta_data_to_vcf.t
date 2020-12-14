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
        q{MIP::Vcfparser}      => [qw{ add_feature_file_meta_data_to_vcf }],
        q{MIP::Vcfparser}      => [qw{ define_select_data_headers  }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfparser qw{ add_feature_file_meta_data_to_vcf };
use MIP::Vcfparser qw{ define_select_data_headers };

diag(   q{Test add_feature_file_meta_data_to_vcf from Feature_file.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $DISEASE_AS_TRNST_ID_COL_NR => 3;

## Given a annotation when not present in vcf header meta data
my $annotation                 = q{Disease_associated_transcript};
my @feature_annotation_columns = $DISEASE_AS_TRNST_ID_COL_NR;
my %feature_data               = define_select_data_headers();
my $feature_file_type          = q{select_file};
my $feature_file_path          = q{a_select_file_path};
my $header =
    q?##INFO=<ID=?
  . q?Disease_associated_transcript?
  . q?,Number=.,Type=String,Description="String taken from ?
  . $feature_file_path . q?">?;
my %meta_data;
$feature_data{present}{$annotation}{INFO} = $header;

add_feature_file_meta_data_to_vcf(
    {
        data_href                      => \%feature_data,
        feature_annotation_columns_ref => \@feature_annotation_columns,
        file_key                       => $feature_file_type,
        meta_data_href                 => \%meta_data,
    }
);

my %expected_meta_data =
  ( $feature_file_type => { INFO => { $annotation => $header, }, }, );

## Then add INFO field using feature data header
is_deeply(
    \%{ $meta_data{$feature_file_type} },
    \%{ $expected_meta_data{$feature_file_type} },
    q{Add arbitrary INFO field to meta data for vcf header}
);

## Given a annotation when present in vcf header meta data
my $existing_annotation = q{Disease_associated_transcript_2};
$feature_data{present}{$existing_annotation}{INFO} = $header;
$meta_data{INFO}{$existing_annotation} = $header;

add_feature_file_meta_data_to_vcf(
    {
        data_href                      => \%feature_data,
        feature_annotation_columns_ref => \@feature_annotation_columns,
        file_key                       => $feature_file_type,
        meta_data_href                 => \%meta_data,
    }
);

is( $meta_data{$feature_file_type}{INFO}{$existing_annotation},
    undef, q{Did not add already present INFO field to meta data for vcf header} );

done_testing();
