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
        q{MIP::File::Format::Vcf} => [qw{ get_transcript_effects }],
        q{MIP::Test::Fixtures}    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Vcf qw{ get_transcript_effects };

diag(   q{Test get_transcript_effects from Vcf.pm v}
      . $MIP::File::Format::Vcf::VERSION
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

## Given transcripts effects and vep format field map
my $allele           = q{A};
my $consequence_term = q{missense_variant};
my $transcript_id    = q{enst_id};
my $hgnc_symbol      = q{ADK};
my @transcript_effects =
  ( $allele, $consequence_term, $transcript_id, $HGNC_ID, $hgnc_symbol );
my %vep_format_field_column = (
    Allele      => $ALLELE_INDEX,
    Consequence => $CONSEQUENCE_INDEX,
    Feature     => $TRANSCRIPT_INDEX,
    HGNC_ID     => $HGNC_ID_INDEX,
    SYMBOL      => $HGNC_SYMBOL_INDEX,
);
my %csq = get_transcript_effects(
    {
        transcript_effects_ref       => \@transcript_effects,
        vep_format_field_column_href => \%vep_format_field_column
    }
);

my %expected_csq = (
    allele            => $allele,
    consequence_field => $consequence_term,
    transcript_id     => $transcript_id,
    hgnc_id           => $HGNC_ID,
    hgnc_symbol       => $hgnc_symbol,
);
## Then return transcript csq fields hash
is_deeply( \%csq, \%expected_csq, q{Got transcript csq fields} );

done_testing();
