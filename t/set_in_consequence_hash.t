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
use MIP::Constants qw{ $COMMA $SPACE %SO_CONSEQUENCE_SEVERITY };
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
        q{MIP::File::Format::Vcf} => [qw{ set_in_consequence_hash }],
        q{MIP::Test::Fixtures}    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Vcf qw{ set_in_consequence_hash };

diag(   q{Test set_in_consequence_hash from Vcf.pm v}
      . $MIP::File::Format::Vcf::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given SO terms, rank and keys to set
my $hgnc_id = 1;
my $allele  = q{A};
my %consequence;
my $consequence_term = q{missense_variant};
my $consequence_rank = $SO_CONSEQUENCE_SEVERITY{$consequence_term}{rank};
my $transcript_id    = q{a_transcript_id};
my %set_key          = (
    most_severe_consequence => $consequence_term,
    most_severe_transcript  => $transcript_id,
    rank                    => $consequence_rank,
);

set_in_consequence_hash(
    {
        allele           => $allele,
        consequence_href => \%consequence,
        hgnc_id          => $hgnc_id,
        set_key_href     => \%set_key,
    }
);

## Then keys should be have values in cosequence hash
my %expected_consequence = (
    $hgnc_id => {
        $allele => {
            most_severe_consequence => $consequence_term,
            most_severe_transcript  => $transcript_id,
            rank                    => $consequence_rank,
        },
    },
);

is_deeply( \%consequence, \%expected_consequence, q{Set score in consequence hash} );

done_testing();
