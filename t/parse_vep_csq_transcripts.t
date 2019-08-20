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
        q{MIP::Vcfparser}      => [qw{ parse_vep_csq_transcripts }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfparser qw{ parse_vep_csq_transcripts };

diag(   q{Test parse_vep_csq_transcripts from Vcfparser.pm v}
      . $MIP::Vcfparser::VERSION
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

## Given a gene, and transcript when not part of feature file
my $allele = q{A};
my %consequence;
my $consequence_term = q{missense_variant};
my $transcript_id    = q{enst_id};
my $hgnc_symbol      = q{ADK};
my %hgnc_map         = ( $HGNC_ID => $hgnc_symbol, );
my $per_gene         = 1;
my $transcript       = join $PIPE,
  ( $allele, $consequence_term, $transcript_id, $HGNC_ID, $hgnc_symbol );
my @transcripts = ($transcript);
my %select_data = ( $HGNC_ID => 1, );
my %vcf_record;
my %vep_format_field_column = (
    Allele      => $ALLELE_INDEX,
    Consequence => $CONSEQUENCE_INDEX,
    Feature     => $TRANSCRIPT_INDEX,
    HGNC_ID     => $HGNC_ID_INDEX,
    SYMBOL      => $HGNC_SYMBOL_INDEX,
);

my $is_ok = parse_vep_csq_transcripts(
    {
        consequence_href             => \%consequence,
        hgnc_map_href                => \%hgnc_map,
        per_gene                     => $per_gene,
        select_data_href             => \%select_data,
        transcripts_ref              => \@transcripts,
        vcf_record_href              => \%vcf_record,
        vep_format_field_column_href => \%vep_format_field_column,
    }
);

## Then return true
ok( $is_ok, q{Parsed vep csq transcrips} );

done_testing();
