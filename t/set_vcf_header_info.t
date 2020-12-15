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
        q{MIP::File::Format::Feature_file} => [qw{ set_vcf_header_info }],
        q{MIP::Vcfparser}                  => [qw{ define_select_data_headers  }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Feature_file qw{ set_vcf_header_info };
use MIP::Vcfparser qw{ define_select_data_headers };

diag(   q{Test set_vcf_header_info from Feature_file.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a header key when not present in feature data
my $header = q?##INFO=<ID=Disease_associated_transcript,Number=.,?;
$header .= q?Type=String,Description="Known pathogenic transcript(s) for gene">?;

my $extract_columns_counter = 0;
my %feature_data            = define_select_data_headers();
my $feature_file_type       = q{select_file};
my $feature_file_path       = q{a_select_file_path};
my $header_key              = q{Not present in feature data};

set_vcf_header_info(
    {
        feature_file_type => $feature_file_type,
        feature_file_path => $feature_file_path,
        header_key        => $header_key,
        meta_data_href    => \%feature_data,
        position          => $extract_columns_counter,
    }
);

my %expected_feature_data = (
    present => {
        $header_key => {
            INFO => q?##INFO=<ID=?
              . $header_key
              . q?,Number=.,Type=String,Description="String taken from ?
              . $feature_file_path . q?">?,
            column_order => $extract_columns_counter,
        },
    },
);

## Then add INFO field using feature data header
is_deeply(
    \%{ $feature_data{present} },
    \%{ $expected_feature_data{present} },
    q{Add arbitrary INFO field using input header}
);

## Given a existing header key
my $existing_header_key = q{Disease_associated_transcript};

## Increment position
$extract_columns_counter++;

set_vcf_header_info(
    {
        feature_file_type => $feature_file_type,
        feature_file_path => $feature_file_path,
        header_key        => $existing_header_key,
        meta_data_href    => \%feature_data,
        position          => $extract_columns_counter,
    }
);

$expected_feature_data{present}{$existing_header_key}{INFO} = $header;
$expected_feature_data{present}{$existing_header_key}{column_order} =
  $extract_columns_counter;

is_deeply(
    \%{ $feature_data{present} },
    \%{ $expected_feature_data{present} },
    q{Add INFO field from predefined header}
);

done_testing();
