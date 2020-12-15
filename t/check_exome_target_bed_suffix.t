#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Reference} => [qw{ check_exome_target_bed_suffix }],
        q{MIP::Test::Fixtures}   => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Reference qw{ check_exome_target_bed_suffix };

diag(   q{Test check_exome_target_bed_suffix from Reference.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log({});

## Given a proper file suffix
my $exome_target_file = catfile( q{path_to_bed_file}, q{test_bed_file.bed} );

my $is_ok = check_exome_target_bed_suffix(
    {
        path           => $exome_target_file,
    }
);

## Then return true
ok( $is_ok, q{Provided file has a .bed extension} );

## Given a not valid file suffix
my $wrong_file_suffix = catfile( q{path_to_bed_file}, q{test_bed_file.gz} );

trap {
    check_exome_target_bed_suffix(
        {
            path           => $wrong_file_suffix,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if wrong file suffix} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message if wrong file suffix} );

done_testing();
