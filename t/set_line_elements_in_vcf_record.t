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
        q{MIP::File::Format::Vcf} => [qw{ set_line_elements_in_vcf_record }],
        q{MIP::Test::Fixtures}    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Vcf qw{ set_line_elements_in_vcf_record };

diag(   q{Test set_line_elements_in_vcf_record from Vcf.pm v}
      . $MIP::File::Format::Vcf::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given
my @line_elements = qw{ 1 1 rsid1 A G . PASS AF=1 GT  0/1 0/1 1/1};
my %vcf_record;
my @vcf_format_columns = (
    q{#CHROM}, qw{ POS ID REF ALT QUAL FILTER INFO FORMAT ACC5339A2 ACC5339A3 ACC5339A4 }
);

set_line_elements_in_vcf_record(
    {
        line_elements_ref      => \@line_elements,
        vcf_record_href        => \%vcf_record,
        vcf_format_columns_ref => \@vcf_format_columns,
    }
);

my %expected_vcf_record;
@expected_vcf_record{@vcf_format_columns} = @line_elements;

## Then
is_deeply( \%vcf_record, \%expected_vcf_record, q{Set line elements to vcf record} );

done_testing();
