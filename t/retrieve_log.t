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
use MIP::Test::Fixtures qw{ test_log };

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Log::MIP_log4perl} => [qw{ retrieve_log }],
        q{MIP::Test::Fixtures}    => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Log::MIP_log4perl qw{ retrieve_log };

diag(   q{Test retrieve_log from MIP_log4perl.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a log name and debug log level
my $log = test_log( {} );

$log = retrieve_log(
    {
        log_name => q{TEST},
        quiet    => 0,
        verbose  => 1,
    }
);

## Then level should be debug
ok( $log->is_debug(), q{Got log for debug level} );

## Given warn log level
$log = retrieve_log(
    {
        log_name => q{TEST},
        quiet    => 1,
        verbose  => 0,
    }
);

## Then level should be warn
ok( $log->is_warn(), q{Got log for warn level} );

## Given warn log level
$log = retrieve_log(
    {
        log_name => q{TEST},
        quiet    => 0,
        verbose  => 0,
    }
);

## Then level should be info
ok( $log->is_info(), q{Got log for info level} );

done_testing();
