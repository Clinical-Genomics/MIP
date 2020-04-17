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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Contigs}        => [qw{ delete_male_contig }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Contigs qw{ delete_male_contig };

diag(   q{Test delete_male_contig from Contigs.pm v}
      . $MIP::Contigs::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given contigs, when male not present
my @contigs = qw{ 1 2 3 4 Y};

my $is_male = 0;

my @no_male_contigs = delete_male_contig(
    {
        contigs_ref      => \@contigs,
        contig_names_ref => [qw{ Y }],
        found_male       => $is_male,
    }
);

## Define expected outcome
my @expected_contigs = qw{ 1 2 3 4 };

## Then removed the male contig
is_deeply( \@no_male_contigs, \@expected_contigs, q{Removed male contig} );

## Given contigs, when male is presnt
$is_male = 1;

my @has_male_contigs = delete_male_contig(
    {
        contigs_ref      => \@contigs,
        contig_names_ref => [qw{ Y }],
        found_male       => $is_male,
    }
);

## Define expected outcome
my @expected_male_contigs = qw{ 1 2 3 4 Y };

## Then did not remove male contig
is_deeply( \@has_male_contigs, \@expected_male_contigs, q{Did not remove male contig} );

done_testing();
