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
        q{MIP::Vcfparser}      => [qw{ parse_vep_csq_schema }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfparser qw{ parse_vep_csq_schema };

diag(   q{Test parse_vep_csq_schema from Vcfparser.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $HGNC_ID_CSQ_FIELD_INDEX => 22;

## Given a VEP CSQ line in meta data when parse_vep is true
my $vep_csq_schema_header =
q{##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|REFSEQ_MATCH|SOURCE|GIVEN_REF|USED_REF|BAM_EDIT|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|ExACpLI|LoFtool|MES-NCSS_downstream_acceptor|MES-NCSS_downstream_donor|MES-NCSS_upstream_acceptor|MES-NCSS_upstream_donor|MES-SWA_acceptor_alt|MES-SWA_acceptor_diff|MES-SWA_acceptor_ref|MES-SWA_acceptor_ref_comp|MES-SWA_donor_alt|MES-SWA_donor_diff|MES-SWA_donor_ref|MES-SWA_donor_ref_comp|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|genomic_superdups_frac_match">};
my %meta_data = ( INFO => { CSQ => $vep_csq_schema_header, }, );
my $most_severe_consequence_header =
q{##INFO=<ID=most_severe_consequence,Number=.,Type=String,Description="Most severe genomic consequence.">};
my %vep_format_field_column;

my $is_ok = parse_vep_csq_schema(
    {
        meta_data_href               => \%meta_data,
        parse_vep                    => 1,
        vep_format_field_column_href => \%vep_format_field_column,
    }
);

## Then return true
ok( $is_ok, q{Parsed VEP CSQ header line} );

## Then HGNC ID key and corresponding index should have been added to hash
is( $vep_format_field_column{HGNC_ID},
    $HGNC_ID_CSQ_FIELD_INDEX, q{Set HGNC ID field in hash} );

## Then "most severe consequence" VCF header line should been added to meta data
is(
    $meta_data{INFO}{most_severe_consequence},
    $most_severe_consequence_header,
    q{Added moste severe consequence header}
);

done_testing();
