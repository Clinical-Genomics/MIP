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
use MIP::Constants qw{ $COMMA $EMPTY_STR $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Analysis}       => [qw{ get_vcf_parser_analysis_suffix }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Analysis qw{ get_vcf_parser_analysis_suffix };

diag(   q{Test get_vcf_parser_analysis_suffix from Analysis.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

Readonly my $TWO => 2;

my @expected_suffixes       = ($EMPTY_STR);
my @expected_sv_suffixes    = ( $EMPTY_STR, q{selected} );
my @expected_panel_suffixes = qw{ all selected };

## Given research outfile
my @analysis_suffixes = get_vcf_parser_analysis_suffix(
    {
        vcfparser_outfile_count => 1,
    }
);

## Then return only $EMPTY_STR suffix
is_deeply( \@analysis_suffixes, \@expected_suffixes, q{Get single analysis suffix} );

## Given research and select outfiles
my @analysis_sv_suffixes = get_vcf_parser_analysis_suffix(
    {
        vcfparser_outfile_count => $TWO,
    }
);

## Then return both suffixes
is_deeply( \@analysis_sv_suffixes, \@expected_sv_suffixes,
    q{Get both analysis suffixes} );

my @analysis_panel_suffixes = get_vcf_parser_analysis_suffix(
    {
        analysis_type           => q{panel},
        vcfparser_outfile_count => $TWO,
    }
);

## Then return both suffixes
is_deeply( \@analysis_panel_suffixes, \@expected_panel_suffixes,
    q{Get both panel suffixes} );

done_testing();
