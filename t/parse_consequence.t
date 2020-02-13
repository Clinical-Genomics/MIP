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
        q{MIP::Vcfparser}      => [qw{ parse_consequence }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfparser qw{ parse_consequence };

diag(   q{Test parse_consequence from Vcfparser.pm v}
      . $MIP::Vcfparser::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given VEP CSQ info and select feaure gene
my $allele                  = q{A};
my @feature_type_keys       = qw{ range select};
my $hgnc_id                 = 1;
my $hgnc_symbol             = q{a_gene_symbol};
my %hgnc_map                = ( $hgnc_id => $hgnc_symbol, );
my $most_severe_consequence = q{missense_variant};
my %most_severe_feature;
my %most_severe_pli;

## Initilize pli score for feature keys
@most_severe_pli{@feature_type_keys} = 0;

my $most_severe_transcript = q{transcript_id};
my %consequence            = (
    $hgnc_id => {
        $allele => {
            most_severe_consequence => $most_severe_consequence,
            most_severe_transcript  => $most_severe_transcript
        },
    },
);
my $per_gene    = 1;
my %pli_score   = ( $hgnc_symbol => 1, );
my %select_data = ( $hgnc_id => 1, );
my %vcf_record;

my $is_ok = parse_consequence(
    {
        consequence_href         => \%consequence,
        hgnc_map_href            => \%hgnc_map,
        most_severe_feature_href => \%most_severe_feature,
        most_severe_pli_href     => \%most_severe_pli,
        per_gene                 => $per_gene,
        pli_score_href           => \%pli_score,
        select_data_href         => \%select_data,
        vcf_record_href          => \%vcf_record,
    }
);

## Then parse VEP consequence
ok( $is_ok, q{Parsed consequence and most severe annotations for feature file(s) } );

done_testing();
