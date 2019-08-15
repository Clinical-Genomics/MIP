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
        q{MIP::Vcfparser}      => [qw{ set_most_severe_pli }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfparser qw{ set_most_severe_pli };

diag(   q{Test set_most_severe_pli from Vcfparser.pm v}
      . $MIP::Vcfparser::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $LOW_PLI_SCORE => 0.001;
Readonly my $MAX_PLI_SCORE => 1.0;

## Given a gene and undef pli score
my $hgnc_symbol     = q{ADK};
my %most_severe_pli = ( range_pli => 0 );
my $pli_score;
my %select_data;

my $is_ok = set_most_severe_pli(
    {
        hgnc_symbol          => $hgnc_symbol,
        most_severe_pli_href => \%most_severe_pli,
        pli_score            => $pli_score,
        select_data_href     => \%select_data,
    }
);

## Then return false
is( $is_ok, undef, q{Return if no pli score} );

## Given a defined pli score when a gene is not in select_data
set_most_severe_pli(
    {
        hgnc_symbol          => $hgnc_symbol,
        most_severe_pli_href => \%most_severe_pli,
        pli_score            => $LOW_PLI_SCORE,
        select_data_href     => \%select_data,
    }
);

my %expected_most_severe_pli = ( range_pli => $LOW_PLI_SCORE );

## Then set most severe pli score for range feature
is_deeply( \%most_severe_pli, \%expected_most_severe_pli,
    q{Set most severe pli score for range feature} );

## Given a defined pli score when a gene is in select_data
$select_data{$hgnc_symbol} = 1;
set_most_severe_pli(
    {
        hgnc_symbol          => $hgnc_symbol,
        most_severe_pli_href => \%most_severe_pli,
        pli_score            => $MAX_PLI_SCORE,
        select_data_href     => \%select_data,
    }
);

$expected_most_severe_pli{range_pli}  = $MAX_PLI_SCORE;
$expected_most_severe_pli{select_pli} = $MAX_PLI_SCORE;

## Then set most severe pli score for range feature
is_deeply( \%most_severe_pli, \%expected_most_severe_pli,
    q{Set most severe pli score for features} );

done_testing();
