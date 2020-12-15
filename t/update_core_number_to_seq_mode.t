#! /usr/bin/env perl

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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Cluster} => [qw{ update_core_number_to_seq_mode }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Cluster qw{ update_core_number_to_seq_mode };

diag(   q{Test update_core_number_to_seq_mode from Cluster.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a core number and sequence run types
Readonly my $CORE_NUMBER => 1;

my @sequence_run_types = qw{ paired-end single-end };

my @returned_core_numbers;

foreach my $sequence_run_type (@sequence_run_types) {

    push @returned_core_numbers,
      update_core_number_to_seq_mode(
        {
            core_number       => $CORE_NUMBER,
            sequence_run_type => $sequence_run_type,
        }
      );
}

## Then update core number or sequence mode
is( $returned_core_numbers[0], 3, q{Updated core number to paired-end mode} );

is( $returned_core_numbers[1], 2, q{Updated core number to single-end mode} );

done_testing();
