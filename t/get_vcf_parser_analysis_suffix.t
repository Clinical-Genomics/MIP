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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA     => q{,};
Readonly my $EMPTY_STR => q{};
Readonly my $SPACE     => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Get::Analysis}  => [qw{ get_vcf_parser_analysis_suffix }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::Analysis qw{ get_vcf_parser_analysis_suffix };

diag(   q{Test get_vcf_parser_analysis_suffix from Analysis.pm v}
      . $MIP::Get::Analysis::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my @expected_suffix    = ($EMPTY_STR);
my @expected_sv_suffix = ( $EMPTY_STR, q{selected} );

## Given research outfile
my $vcfparser_outfile_count = 1;
my @analysis_suffix         = get_vcf_parser_analysis_suffix(
    { vcfparser_outfile_count => $vcfparser_outfile_count, } );

## Then return only $EMPTY_STR suffix
is_deeply( \@analysis_suffix, \@expected_suffix, q{Get single analysis suffix} );

## Given research and select outfiles
my $sv_vcfparser_outfile_count = 2;
my @analysis_sv_suffix         = get_vcf_parser_analysis_suffix(
    { vcfparser_outfile_count => $sv_vcfparser_outfile_count, } );

## Then return both suffixes
is_deeply( \@analysis_sv_suffix, \@expected_sv_suffix, q{Get both analysis suffixes} );

done_testing();
