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
        q{MIP::File::Format::Vcf} => [qw{ set_info_key_pairs_in_vcf_record }],
        q{MIP::Test::Fixtures}    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Vcf qw{ set_info_key_pairs_in_vcf_record };

diag(   q{Test set_info_key_pairs_in_vcf_record from Vcf.pm v}
      . $MIP::File::Format::Vcf::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given an VCF INFO record element
my %vcf_record = ( INFO => q{AF=1;HGNC_ID=ADK}, );

set_info_key_pairs_in_vcf_record( { vcf_record_href => \%vcf_record, } );

my %expected_vcf_record = (
    INFO           => q{AF=1;HGNC_ID=ADK},
    INFO_key_value => {
        AF      => 1,
        HGNC_ID => q{ADK},
    },
);

## Then key-value pairs should be added under INFO_key_value
is_deeply( \%vcf_record, \%expected_vcf_record, q{Set INFO key-value pair to hash} );

done_testing();
