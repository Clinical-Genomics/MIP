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
        q{MIP::List}           => [qw{ get_splitted_lists }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::List qw{ get_splitted_lists };

diag(   q{Test get_splitted_lists from List.pm v}
      . $MIP::List::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given bam contigs
my @grch38_contigs = qw{
  chr1 chr2 chr3 chr4 chr5 chr6
  chr7 chr8 chr9 chr10 chr11 chr12
  chr13 chr14 chr15 chr16 chr17 chr18
  chr19 chr20 chr21 chr22 chrX chrY
  chrM };

my ( $match_contigs_ref, $no_match_contigs_ref ) = get_splitted_lists(
    {
        regexp   => qr/M/,
        list_ref => \@grch38_contigs,
    }
);

my @expected_match    = qw{ chrM };
my @expected_no_match = qw{
  chr1 chr2 chr3 chr4 chr5 chr6
  chr7 chr8 chr9 chr10 chr11 chr12
  chr13 chr14 chr15 chr16 chr17 chr18
  chr19 chr20 chr21 chr22 chrX chrY };

## Then
is_deeply( $match_contigs_ref,    \@expected_match,    q{Return match list} );
is_deeply( $no_match_contigs_ref, \@expected_no_match, q{Return no match list} );

done_testing();
