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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::Language::Perl} => [qw{ check_modules_existance }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Language::Perl qw{ check_modules_existance };

diag(   q{Test check_modules_existance from Perl.pm v}
      . $MIP::Language::Perl::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a perl module that exists
my @modules = qw{ warnings };

## When checking if perl module exists
my $was_found = check_modules_existance(
    {
        modules_ref  => \@modules,
        program_name => $PROGRAM_NAME,
    }
);

## Then return
ok( $was_found, q{Required perl module} );

## Given faulty module
@modules = qw{ not_a_perl_module };

## When checking if perl module exists
trap {
    check_modules_existance(
        {
            modules_ref  => \@modules,
            program_name => $PROGRAM_NAME,
        }
    )
};

## Then exit and throw FATAL log message
like( $trap->stderr, qr/FATAL/xms, q{Could not find module - croaked and exited} );

done_testing();
