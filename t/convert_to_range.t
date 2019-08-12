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
use MIP::Constants qw{ $COMMA $DOT $SPACE };
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
        q{MIP::File::Format::Vcf} => [qw{ convert_to_range }],
        q{MIP::Test::Fixtures}    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Vcf qw{ convert_to_range };

diag(   q{Test convert_to_range from Vcf.pm v}
      . $MIP::File::Format::Vcf::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $ALT_ALLELE_FIELD_INDEX => 4;
Readonly my $REF_ALLELE_FIELD_INDEX => 3;
Readonly my $CHROM_NR               => 1;
Readonly my $DELETION_STOP          => 3_897_459;
Readonly my $INSERTION_STOP         => 3_897_460;
Readonly my $SNV_INS_STOP           => 3_897_459;
Readonly my $SNV_STOP               => 3_897_457;
Readonly my $START                  => 3_897_456;

## Given no alternative alleles
my @no_variant_elements = ( $CHROM_NR, $START, $DOT, q{A}, $DOT );

my $stop = convert_to_range(
    {
        alt_allele_field => $no_variant_elements[$ALT_ALLELE_FIELD_INDEX],
        reference_allele => $no_variant_elements[$REF_ALLELE_FIELD_INDEX],
        start_position   => $no_variant_elements[1],
    }
);

## Then stop should be undef
is( $stop, undef, q{Stop position when no alt allele} );

## Given a SNV when single alternative allele
my @snv_single_variant_elements = ( $CHROM_NR, $START, $DOT, qw{ A G } );

$stop = convert_to_range(
    {
        alt_allele_field => $snv_single_variant_elements[$ALT_ALLELE_FIELD_INDEX],
        reference_allele => $snv_single_variant_elements[$REF_ALLELE_FIELD_INDEX],
        start_position   => $snv_single_variant_elements[1],
    }
);

## Then stop should equal
is( $stop, $SNV_STOP, q{SNV stop position when single allele} );

## Given a SNV when multiple alternative SNV alleles
my @snv_multiple_variant_elements = ( $CHROM_NR, $START, $DOT, q{A}, q{G,T} );

$stop = convert_to_range(
    {
        alt_allele_field => $snv_multiple_variant_elements[$ALT_ALLELE_FIELD_INDEX],
        reference_allele => $snv_multiple_variant_elements[$REF_ALLELE_FIELD_INDEX],
        start_position   => $snv_multiple_variant_elements[1],
    }
);

## Then stop should equal snv stop
is( $stop, $SNV_STOP, q{SNV stop position when multiple SNV alleles} );

## Given a SNV/insertion when multiple alternative SNV/insertion (GTC) alleles
my @snv_insert_variant_elements = ( $CHROM_NR, $START, $DOT, q{A}, q{G,AGTC,T} );

$stop = convert_to_range(
    {
        alt_allele_field => $snv_insert_variant_elements[$ALT_ALLELE_FIELD_INDEX],
        reference_allele => $snv_insert_variant_elements[$REF_ALLELE_FIELD_INDEX],
        start_position   => $snv_insert_variant_elements[1],
    }
);

## Then stop should shifted 3 bases downstream
is( $stop, $SNV_INS_STOP, q{SNV/insertion stop position when multiple SNV alleles} );

## Given a microsatelite when a deletion of 2 bases (TC) and an insertion of one base (T)
my @microsatelite_variant_elements = ( $CHROM_NR, $START, $DOT, q{GTC}, q{G,GTCT} );

$stop = convert_to_range(
    {
        alt_allele_field => $microsatelite_variant_elements[$ALT_ALLELE_FIELD_INDEX],
        reference_allele => $microsatelite_variant_elements[$REF_ALLELE_FIELD_INDEX],
        start_position   => $microsatelite_variant_elements[1],
    }
);

## Then stop should shifted 3 base downstream
is( $stop, $DELETION_STOP,
    q{Microsatelite stop position when multiple SNV alleles with longest deletion} );

## Given a microsatelite when a deletion of 1 bases (T) and an insertion of 4 base (TCTA)
@microsatelite_variant_elements = ( $CHROM_NR, $START, $DOT, q{GT}, q{A,G,GTCTA} );

$stop = convert_to_range(
    {
        alt_allele_field => $microsatelite_variant_elements[$ALT_ALLELE_FIELD_INDEX],
        reference_allele => $microsatelite_variant_elements[$REF_ALLELE_FIELD_INDEX],
        start_position   => $microsatelite_variant_elements[1],
    }
);

## Then stop should shifted 3 base downstream
is( $stop, $INSERTION_STOP,
    q{Microsatelite stop position when multiple SNV alleles with long insertion} );

done_testing();
