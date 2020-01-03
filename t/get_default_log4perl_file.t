#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname basename fileparse };
use File::Spec::Functions qw{ catdir catfile };
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
        q{MIP::Log::MIP_log4perl} => [qw{ get_default_log4perl_file }],
        q{MIP::Test::Fixtures}    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Log::MIP_log4perl qw{ get_default_log4perl_file };

diag(   q{Test get_default_log4perl_file from MIP_log4perl.pm v}
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
my $test_dir      = File::Temp->newdir();
my $test_log_path = catfile( $test_dir, q{test.log} );

# Catches script name and removes ending
my $script = fileparse( basename( $PROGRAM_NAME, $DOT . q{t} ) );

my %active_parameter = ( log_file => undef, );

## Given no log file from user
$active_parameter{log_file} = get_default_log4perl_file(
    {
        cmd_input       => $active_parameter{log_file},
        date            => $date,
        date_time_stamp => $date_time_stamp,
        outdata_dir     => catfile($test_dir),
        script          => $script,
    }
);

## Then get the default Log4perl file using supplied dynamic parameters
ok( $active_parameter{log_file}, q{Got default log file} );

## Given a log file from user
$active_parameter{log_file} = $test_log_path;

$active_parameter{log_file} = get_default_log4perl_file(
    {
        cmd_input       => $active_parameter{log_file},
        script          => $script,
        date            => $date,
        date_time_stamp => $date_time_stamp,
    }
);

## Then get the Log4perl file using cmd input
is( $active_parameter{log_file}, $test_log_path, q{Got user log file} );

done_testing();
