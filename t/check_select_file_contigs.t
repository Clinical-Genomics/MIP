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
use Test::Trap;

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
        q{MIP::Contigs}        => [qw{ check_select_file_contigs }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Contigs qw{ check_select_file_contigs };

diag(   q{Test check_select_file_contigs from Contigs.pm v}
      . $MIP::Contigs::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given select and primary contigs
my @primary_contigs     = qw{ chr1 chr2 chr3 };
my @select_file_contigs = qw{ chr1 chr2 };

my $is_ok = check_select_file_contigs(
    {
        contigs_ref             => \@primary_contigs,
        select_file_contigs_ref => \@select_file_contigs,
    }
);

## Then return true
ok( $is_ok, q{Select contigs is subset of priamry contigs} );

## Given select and primary contigs when not matching
push @select_file_contigs, q{chr_unknown};
trap {
    check_select_file_contigs(
        {
            contigs_ref             => \@primary_contigs,
            select_file_contigs_ref => \@select_file_contigs,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if not a primary contig subset } );
like( $trap->stderr, qr/Is\s+not\s+a\s+subset/xms, q{Throw fatal log message} );

done_testing();
