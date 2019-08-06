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
use Modern::Perl qw{ 2017 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $NEWLINE $SPACE };
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
        q{MIP::File::Format::Vcf} => [qw{ parse_vcf_header }],
        q{MIP::Vcfparser}         => [qw{ write_meta_data }],
        q{MIP::Test::Fixtures}    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Vcf qw{ parse_vcf_header };
use MIP::Vcfparser qw{ write_meta_data };

diag(   q{Test write_meta_data from Vcfparser.pm v}
      . $MIP::Vcfparser::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a unsorted vcf file header
my $vcf_path = catfile( $Bin, qw{ data test_data test.vcf } );

open my $VCF_FH, q{<}, $vcf_path
  or croak( q{Cannot open } . $vcf_path . $COLON . $OS_ERROR, $NEWLINE );

my %meta_data;

## Read file and add header to meta data
LINE:
while (<$VCF_FH>) {

    chomp;

    ## Unpack line
    my $line = $_;

    ## Skip blank lines
    next LINE if ( $line =~ /^\s+$/sxm );

    ## Header meta data
    if ( $line =~ /\A [#]{2}/sxm ) {

        parse_vcf_header(
            {
                meta_data_href   => \%meta_data,
                meta_data_string => $line,
            }
        );
        next;
    }
}

close $VCF_FH;

## Add specific feature file annotations
$meta_data{range}{INFO}{Reduced_penetrance} =
q{##INFO=<ID=Reduced_penetrance,Number=.,Type=String,Description="Pathogenic gene which can exhibit reduced penetrance">};
$meta_data{select}{INFO}{HGNC_ID} =
  q{##INFO=<ID=HGNC_symbol,Number=.,Type=String,Description="The HGNC gene symbol">};

my $file_content;
my $select_file_content;

## Open files for writing
open my $SELECT_VCF_FH, q{>}, \$select_file_content
  or croak( q{Cannot open } . $select_file_content . $COLON . $OS_ERROR, $NEWLINE );

open my $VCF_OUT_FH, q{>}, \$file_content
  or croak( q{Cannot open } . $file_content . $COLON . $OS_ERROR, $NEWLINE );

write_meta_data(
    {
        FILEHANDLE       => $VCF_OUT_FH,
        meta_data_href   => \%meta_data,
        SELECTFILEHANDLE => $SELECT_VCF_FH,
    }
);

close $SELECT_VCF_FH;
close $VCF_OUT_FH;

my $expected_header = <<'EOF';
##fileformat=VCFv4.2
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=Reduced_penetrance,Number=.,Type=String,Description="Pathogenic gene which can exhibit reduced penetrance">
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##source=HaplotypeCaller
##bcftools_normVersion=1.9+htslib-1.9
##Software=<ID=vcfparser.pl,Version=1.2.16,Date=2019-08-05
EOF

my $expected_select_header = <<'EOF';
##fileformat=VCFv4.2
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=HGNC_symbol,Number=.,Type=String,Description="The HGNC gene symbol">
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##source=HaplotypeCaller
##bcftools_normVersion=1.9+htslib-1.9
##Software=<ID=vcfparser.pl,Version=1.2.16,Date=2019-08-05
EOF

## Then header should be written to vcf out file
is( $file_content, $expected_header, q{Wrote meta-data for vcf outfile} );

is( $select_file_content, $expected_select_header,
    q{Wrote meta-data for select vcf outfile} );

done_testing();
