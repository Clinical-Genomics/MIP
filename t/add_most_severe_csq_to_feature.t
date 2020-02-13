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
        q{MIP::Vcfparser}      => [qw{ add_most_severe_csq_to_feature }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfparser qw{ add_most_severe_csq_to_feature };

diag(   q{Test add_most_severe_csq_to_feature from Vcfparser.pm v}
      . $MIP::Vcfparser::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a hgnc_id not part of select data when per gene is true
my $hgnc_id                 = 1;
my $most_severe_consequence = q{missense_variant};
my %most_severe_feature;
my $most_severe_transcript = q{transcript_id};
my %vcf_record;
my %select_data;

add_most_severe_csq_to_feature(
    {
        hgnc_id                  => $hgnc_id,
        most_severe_consequence  => $most_severe_consequence,
        most_severe_feature_href => \%most_severe_feature,
        most_severe_transcript   => $most_severe_transcript,
        per_gene                 => 1,
        vcf_record_href          => \%vcf_record,
        select_data_href         => \%select_data,
    }
);

my %expected_most_severe_feature = ( range => [ $most_severe_consequence, ], );
my %expected_vcf_record = ( range_transcripts => [ $most_severe_transcript, ], );

## Then add most severe consequence to range feature
is_deeply(
    \%most_severe_feature,
    \%expected_most_severe_feature,
    q{Added most severe consequence to range feature}
);

## Then add most severe transcript to vcf record
is_deeply( \%vcf_record, \%expected_vcf_record,
    q{Added most severe transcript to vcf record} );

## Given a select data with a hgnc_id when perl gene is true
$select_data{$hgnc_id} = 1;

## Clean-up from previous test
delete $vcf_record{range_transcripts};
delete $most_severe_feature{range};

add_most_severe_csq_to_feature(
    {
        hgnc_id                  => $hgnc_id,
        most_severe_consequence  => $most_severe_consequence,
        most_severe_feature_href => \%most_severe_feature,
        most_severe_transcript   => $most_severe_transcript,
        per_gene                 => 1,
        vcf_record_href          => \%vcf_record,
        select_data_href         => \%select_data,
    }
);

$expected_most_severe_feature{select}    = [$most_severe_consequence];
$expected_vcf_record{select_transcripts} = [$most_severe_transcript];

## Then add most severe consequence to features
is_deeply(
    \%most_severe_feature,
    \%expected_most_severe_feature,
    q{Added most severe consequence to features}
);

## Then add most severe transcript to vcf record for select feature
is_deeply( \%vcf_record, \%expected_vcf_record,
    q{Added most severe transcript to vcf record for select feature} );
done_testing();
