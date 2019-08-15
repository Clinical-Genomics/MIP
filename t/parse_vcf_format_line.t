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
use MIP::Constants qw{ $COLON $COMMA $NEWLINE $SPACE $TAB };
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
        q{MIP::Vcfparser}      => [qw{ parse_vcf_format_line }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfparser qw{ parse_vcf_format_line };

diag(   q{Test parse_vcf_format_line from Vcfparser.pm v}
      . $MIP::Vcfparser::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given format line and filhandles
my @vcf_format_columns = (
    q{#CHROM}, qw{ POS ID REF ALT QUAL FILTER INFO FORMAT ACC5339A2 ACC5339A3 ACC5339A4 }
);
my $vcf_format_line = join $TAB, @vcf_format_columns;
my $file_content;
my $select_file_content;

## Open files for writing
open my $SELECT_VCF_FH, q{>}, \$select_file_content
  or croak( q{Cannot open } . $select_file_content . $COLON . $OS_ERROR, $NEWLINE );

open my $VCF_OUT_FH, q{>}, \$file_content
  or croak( q{Cannot open } . $file_content . $COLON . $OS_ERROR, $NEWLINE );

my @ret_format_columns = parse_vcf_format_line(
    {
        FILEHANDLE       => $VCF_OUT_FH,
        format_line      => $vcf_format_line,
        SELECTFILEHANDLE => $SELECT_VCF_FH,
    }
);

close $SELECT_VCF_FH;
close $VCF_OUT_FH;

my $expected_vcf_format_line = $vcf_format_line . $NEWLINE;

## Then vcf format line columns should be returned
is_deeply( \@ret_format_columns, \@vcf_format_columns, q{Got vcf format columns} );

## Then vcf format line should be written to vcf outfile and select file
is( $file_content, $expected_vcf_format_line, q{Wrote vcf format line to vcf outfile} );

is( $select_file_content, $expected_vcf_format_line,
    q{Wrote vcf format line to select vcf outfile} );

done_testing();
