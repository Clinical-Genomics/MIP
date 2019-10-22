#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
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
use MIP::Constants qw{ $COMMA $COLON $NEWLINE $SPACE };
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
        q{MIP::File::Format::Vcf} => [qw{ parse_vcf_header }],
        q{MIP::Test::Fixtures}    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Vcf qw{ parse_vcf_header };

diag(   q{Test parse_vcf_header from Vcf.pm v}
      . $MIP::File::Format::Vcf::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given
my $vcf_path = catfile( $Bin, qw{ data test_data 643594-miptest_gent_vrecal.vcf } );
my %header_data;
## Set data under these key value pairs (i.e. $vcf_schema_key => vcf_id  or
## $vcf_schema_key => $vcf_schema_key)
my %vcf_schema_key = (
    contig => q{contig},
    other  => q{other},
);

my %vcf_schema_and_id_key = (
    ALT        => q{NON_REF},
    fileformat => q{fileformat},
    FILTER     => q{PASS},         # $vcf_schema_key => vcf_id
    FORMAT     => q{AD},
    INFO       => q{AC},
);

open my $filehandle, q{<}, $vcf_path
  or croak( q{Cannot open } . $vcf_path . $COLON . $OS_ERROR, $NEWLINE );
my $is_ok;

LINE:
while (<$filehandle>) {

    chomp;

    my $line = $_;

    if ( $line =~ /\A [#]{2}/sxm ) {

        $is_ok = parse_vcf_header(
            {
                meta_data_href   => \%header_data,
                meta_data_string => $line,
            }
        );
    }
}
close $filehandle;
my %expected_header_data = (
    ALT => {
        NON_REF =>
q{##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">},
    },
    contig     => { contig     => [q{##contig=<ID=1,length=249250621>}], },
    fileformat => { q{VCFv4.2} => q{##fileformat=VCFv4.2}, },
    FILTER => { q{PASS} => q{##FILTER=<ID=PASS,Description="All filters passed">}, },
    FORMAT => {
        q{AD} =>
q{##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">},
    },
    INFO => {
        q{AC} =>
q{##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">},
    },
    other => { other => [q{##source=ApplyVQSR}], },
);

## Then sub should return true
ok( $is_ok, q{Parsed vcf header} );

while ( my ($vcf_schema) = each %vcf_schema_key ) {

    is(
        $header_data{$vcf_schema}{$vcf_schema}[0],
        $expected_header_data{$vcf_schema}{$vcf_schema}[0],
        qq{Added header key and data $vcf_schema to hash}
    );
}

while ( my ( $vcf_schema, $vcf_id ) = each %vcf_schema_and_id_key ) {

    is(
        $header_data{$vcf_schema}{$vcf_id},
        $expected_header_data{$vcf_schema}{$vcf_id},
        qq{Added header key and data $vcf_schema to hash}
    );
}
done_testing();
