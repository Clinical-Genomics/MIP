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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Contigs}        => [qw{ delete_contig_elements }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Contigs qw{ delete_contig_elements };

diag(   q{Test delete_contig_elements from Contigs.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given contigs, when no prefix
my @contigs        = qw{ 1 2 3 4 Y MT };
my @remove_contigs = qw{ 2 4 Y MT };

my @cleansed_contigs = delete_contig_elements(
    {
        contigs_ref        => \@contigs,
        remove_contigs_ref => \@remove_contigs,
    }
);

## Define expected outcome
my @expected_contigs = qw{ 1 3 };

## Then remove the contigs
is_deeply( \@cleansed_contigs, \@expected_contigs, q{Removed contigs} );

## Given contigs, when prefix
my @chr_contigs = qw{ chr1 chr2 chr3 chr4 chrY chrM };

@cleansed_contigs = delete_contig_elements(
    {
        contigs_ref        => \@chr_contigs,
        remove_contigs_ref => \@remove_contigs,
    }
);

## Define expected outcome
@expected_contigs = qw{ chr1 chr3 };

## Then remove the contigs irrespective of prefix
is_deeply( \@cleansed_contigs, \@expected_contigs, q{Removed contigs with chr prefix} );

## Given contigs, when prefix in remove
@chr_contigs = qw{ chr1 chr2 chr3 chr4 chrY chrM };
my @chr_remove_contigs = qw{ chr1 chrY chrM };

@cleansed_contigs = delete_contig_elements(
    {
        contigs_ref        => \@chr_contigs,
        remove_contigs_ref => \@chr_remove_contigs,
    }
);

## Define expected outcome
@expected_contigs = qw{ chr2 chr3 chr4 };

## Then remove the contigs irrespective of prefix
is_deeply( \@cleansed_contigs, \@expected_contigs,
    q{Removed contigs independent of prefix} );

done_testing();
