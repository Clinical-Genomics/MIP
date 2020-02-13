#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname basename fileparse };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use Time::Piece;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $DOT $COMMA $SPACE };
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
        q{MIP::Log::MIP_log4perl} => [qw{ get_log }],
        q{MIP::Test::Fixtures}    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Log::MIP_log4perl qw{ get_log };

diag(   q{Test get_log from MIP_log4perl.pm v}
      . $MIP::Log::MIP_log4perl::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Add date_time_stamp
my $date_time       = localtime;
my $date_time_stamp = $date_time->datetime;
my $date            = $date_time->ymd;

## Create temp logger
my $test_dir = File::Temp->newdir();

# Catches script name and removes ending
my $script = fileparse( basename( $PROGRAM_NAME, $DOT . q{t} ) );

my %active_parameter = (
    log_file    => undef,
    outdata_dir => catdir($test_dir),
);

## Given no log file from user and a log name
my $log = get_log(
    {
        active_parameter_href => \%active_parameter,
        date                  => $date,
        date_time_stamp       => $date_time_stamp,
        log_name              => uc q{mip_analyse},
        script                => $script,
    }
);

## Then get log object using default log file
ok( $log, q{Got log} );

done_testing();
