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
use MIP::Constants qw{ $COMMA $PIPE $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Vcfparser}      => [qw{ parse_vep_csq }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfparser qw{ parse_vep_csq };

diag(   q{Test parse_vep_csq from Vcfparser.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $ALLELE_INDEX      => 0;
Readonly my $CONSEQUENCE_INDEX => 1;
Readonly my $TRANSCRIPT_INDEX  => 2;
Readonly my $HGNC_ID           => 257;
Readonly my $HGNC_ID_INDEX     => 3;
Readonly my $HGNC_SYMBOL_INDEX => 4;

## Given VEP CSQ info and select feaure gene
my $allele           = q{A};
my $hgnc_symbol      = q{ADK};
my %hgnc_map         = ( $HGNC_ID => $hgnc_symbol, );
my $consequence_term = q{missense_variant};
my %consequence;
my $per_gene      = 1;
my %pli_score     = ( $hgnc_symbol => 1, );
my %select_data   = ( $HGNC_ID => 1, );
my $transcript_id = q{enst_id};
my $transcript    = join $PIPE,
  ( $allele, $consequence_term, $transcript_id, $HGNC_ID, $hgnc_symbol );
my @transcripts             = ($transcript);
my %vcf_record              = ( INFO_key_value => { CSQ => [$transcript], } );
my %vep_format_field_column = (
    Allele      => $ALLELE_INDEX,
    Consequence => $CONSEQUENCE_INDEX,
    Feature     => $TRANSCRIPT_INDEX,
    HGNC_ID     => $HGNC_ID_INDEX,
    SYMBOL      => $HGNC_SYMBOL_INDEX,
);

my $is_ok = parse_vep_csq(
    {
        consequence_href             => \%consequence,
        per_gene                     => $per_gene,
        pli_score_href               => \%pli_score,
        record_href                  => \%vcf_record,
        select_data_href             => \%select_data,
        vep_format_field_column_href => \%vep_format_field_column,
    }
);

## Then parse VEP consequence
ok( $is_ok, q{Parsed vep consequence } );

done_testing();
