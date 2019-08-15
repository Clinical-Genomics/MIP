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
        q{MIP::Vcfparser}      => [qw{ add_transcript_to_feature_file }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfparser qw{ add_transcript_to_feature_file };

diag(   q{Test add_transcript_to_feature_file from Vcfparser.pm v}
      . $MIP::Vcfparser::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a gene, and transcript when not part of feature file
my $transcript = q{a transcript id};
my %vcf_record;

add_transcript_to_feature_file(
    {
        vcf_record_href => \%vcf_record,
        transcript      => $transcript,
    }
);

my %expected_vcf_record = ( range_transcripts => [ $transcript, ], );

## Then transcript should be added to range features only
is_deeply( \%vcf_record, \%expected_vcf_record,
    q{Added transcript to rang feature in hash when no select data hash} );

## Given a select data hash
my %select_data;

## Clean-up from previous test
delete $vcf_record{range_transcripts}[0];

add_transcript_to_feature_file(
    {
        vcf_record_href  => \%vcf_record,
        select_data_href => \%select_data,
        transcript       => $transcript,
    }
);

%expected_vcf_record = ( range_transcripts => [ $transcript, ], );

## Then transcript should be added to range features only
is_deeply( \%vcf_record, \%expected_vcf_record,
    q{Added transcript to rang feature in hash when supplied hgnc_id in select data} );

## Given a select data hash when undef $hgnc_id
my $hgnc_id;

## Clean-up from previous test
delete $vcf_record{range_transcripts}[0];

add_transcript_to_feature_file(
    {
        hgnc_id          => $hgnc_id,
        vcf_record_href  => \%vcf_record,
        select_data_href => \%select_data,
        transcript       => $transcript,
    }
);

## Then transcript should be added to range features only
is_deeply( \%vcf_record, \%expected_vcf_record,
    q{Added transcript to rang feature in hash when supplied undef hgnc_id in select data}
);

## Given a transcript with gene in select feature
$hgnc_id = 2;
$select_data{$hgnc_id} = 1;

## Clean-up from previous test
delete $vcf_record{range_transcripts}[0];

add_transcript_to_feature_file(
    {
        hgnc_id          => $hgnc_id,
        vcf_record_href  => \%vcf_record,
        select_data_href => \%select_data,
        transcript       => $transcript,
    }
);

push @{ $expected_vcf_record{select_transcripts} }, $transcript;

## Then transcript should be added to range and select feature
is_deeply( \%vcf_record, \%expected_vcf_record,
    q{Added transcript to rang and select feature in hash} );

done_testing();
