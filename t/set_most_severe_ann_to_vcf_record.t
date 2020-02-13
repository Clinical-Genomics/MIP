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
use MIP::Constants qw{ $COMMA $SPACE $UNDERSCORE };
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
        q{MIP::Vcfparser}      => [qw{ set_most_severe_ann_to_vcf_record }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfparser qw{ set_most_severe_ann_to_vcf_record };

diag(   q{Test set_most_severe_ann_to_vcf_record from Vcfparser.pm v}
      . $MIP::Vcfparser::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given most severe annotations when only range feature
my @feature_type_keys            = qw{ range select};
my $most_severe_consequence_term = q{most_severe_consequence=A|missense_variant|1};
my $most_severe_pli_score        = q{most_severe_pli=1};
my %most_severe_feature          = ( range => [$most_severe_consequence_term], );
my %most_severe_pli              = ( range => $most_severe_pli_score, );
my %vcf_record;
my $vcf_key_range = join $UNDERSCORE,
  ( qw{INFO addition}, $feature_type_keys[0], qw{ feature } );

set_most_severe_ann_to_vcf_record(
    {
        feature_type_keys_ref    => \@feature_type_keys,
        most_severe_feature_href => \%most_severe_feature,
        most_severe_pli_href     => \%most_severe_pli,
        vcf_record_href          => \%vcf_record,
    }
);

my %expected_vcf_record = (
    $vcf_key_range => {
        most_severe_pli         => $most_severe_pli_score,
        most_severe_consequence => $most_severe_consequence_term,
    },
);

## Then set annotations to only range feature in vcf record
is_deeply( \%vcf_record, \%expected_vcf_record,
    q{Set most severe annotations to vcf_record for range feature} );

## Given most severe annotations when range and select feature
# Clean-up from previous test
%vcf_record = ();

$most_severe_feature{select} = [$most_severe_consequence_term];
$most_severe_pli{select}     = $most_severe_pli_score;

my $vcf_key_select = join $UNDERSCORE,
  ( qw{INFO addition}, $feature_type_keys[1], qw{ feature } );

set_most_severe_ann_to_vcf_record(
    {
        feature_type_keys_ref    => \@feature_type_keys,
        most_severe_feature_href => \%most_severe_feature,
        most_severe_pli_href     => \%most_severe_pli,
        vcf_record_href          => \%vcf_record,
    }
);

$expected_vcf_record{$vcf_key_select}{most_severe_pli} = $most_severe_pli_score;
$expected_vcf_record{$vcf_key_select}{most_severe_consequence} =
  $most_severe_consequence_term;

## Then set annotations to only range feature in vcf record
is_deeply( \%vcf_record, \%expected_vcf_record,
    q{Set most severe annotations to vcf_record for range and select feature} );

done_testing();
