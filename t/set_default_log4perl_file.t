#!/usr/bin/env perl

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ dirname basename fileparse };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
use FindBin qw{ $Bin };
use Getopt::Long;
use Params::Check qw{ check allow last_error };
use Test::More;
use Time::Piece;
use utf8;
use warnings qw{ FATAL utf8 };
use 5.018;

## CPANM
use autodie;
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $COMMA   => q{,};
Readonly my $DOT     => q{.};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

### User Options
GetOptions(

    # Display help text
    q{h|help} => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },

    # Display version number
    q{v|version} => sub {
        done_testing();
        say {*STDOUT} $NEWLINE
          . basename($PROGRAM_NAME)
          . $SPACE
          . $VERSION
          . $NEWLINE;
        exit;
    },
    q{vb|verbose} => $VERBOSE,
  )
  or (
    done_testing(),
    help(
        {
            USAGE     => $USAGE,
            exit_code => 1,
        }
    )
  );

BEGIN {

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Script::Utils} => [qw{ help }], );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Log::MIP_log4perl});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Log::MIP_log4perl qw{ set_default_log4perl_file };

diag(   q{Test set_default_log4perl_file from MIP_log4perl.pm v}
      . $MIP::Log::MIP_log4perl::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Add date_time_stamp for later use in log and qc_metrics yaml file
my $date_time       = localtime;
my $date_time_stamp = $date_time->datetime;
my $date            = $date_time->ymd;

## Create temp logger
my $test_dir = File::Temp->newdir();
my $test_log_path = catfile( $test_dir, q{test.log} );

# Catches script name and removes ending
my $script = fileparse( basename( $PROGRAM_NAME, $DOT . q{t} ) );

my %active_parameter = (
    log_file    => undef,
    outdata_dir => $test_dir
);

## Set the default Log4perl file using supplied dynamic parameters.
$active_parameter{log_file} = set_default_log4perl_file(
    {
        active_parameter_href => \%active_parameter,
        cmd_input             => $active_parameter{log_file},
        script                => $script,
        date                  => $date,
        date_time_stamp       => $date_time_stamp,
    }
);

## Test

ok( $active_parameter{log_file}, q{Set default log file} );

## Reset for new test
$active_parameter{log_file} = $test_log_path;

## Set the default Log4perl file using supplied dynamic parameters.
$active_parameter{log_file} = set_default_log4perl_file(
    {
        active_parameter_href => \%active_parameter,
        cmd_input             => $active_parameter{log_file},
        script                => $script,
        date                  => $date,
        date_time_stamp       => $date_time_stamp,
    }
);

is( $active_parameter{log_file}, $test_log_path, q{Did not set default} );

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

## Function  : Build the USAGE instructions
## Returns   :
## Arguments : $program_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $program_name;

    my $tmpl = {
        program_name => {
            default     => basename($PROGRAM_NAME),
            strict_type => 1,
            store       => \$program_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
