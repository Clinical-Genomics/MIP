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
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Environment::User} => [qw{ check_email_address }],
        q{MIP::Test::Fixtures}    => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::User qw{ check_email_address };

diag(   q{Test check_email_address from User.pm v}
      . $MIP::Environment::User::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

## Given valid email
my $valid_email = q{mip@scilifelab.se};

my $is_ok = check_email_address(
    {
        email => $valid_email,
    }
);
## Then return true
ok( $is_ok, q{Valid email} );

## Given not valid email
my $not_valid_email = q{mip@clinical_genomics.se};

trap {
    check_email_address(
        {
            email => $not_valid_email,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if not valid email} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message if not valid email} );

## Given undefined email
my $not_defined_email = undef;

my $is_undef = check_email_address(
    {
        email => $not_defined_email,
    }
);

## Then return false
is( $is_undef, undef, q{Returned if email is undefined} );

done_testing();
