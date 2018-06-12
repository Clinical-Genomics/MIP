#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
use FindBin qw{ $Bin };
use Getopt::Long;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use Test::Trap;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2014 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $COMMA   => q{,};
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
    ## Modules with import
    my %perl_module = (
        q{MIP::Log::MIP_log4perl} => [qw{ initiate_logger }],
        q{MIP::Script::Utils}     => [qw{ help }],
    );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

    ## Modules
    my @modules = (q{MIP::Get::Parameter});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Get::Parameter qw{ get_programs_for_shell_installation };

diag(   q{Test get_programs_for_shell_installation from Parameter.pm v}
      . $MIP::Get::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create temp logger
my $test_dir = File::Temp->newdir();
my $test_log_path = catfile( $test_dir, q{test.log} );

## Creates log object
my $log = initiate_logger(
    {
        file_path => $test_log_path,
        log_name  => q{TEST},
    }
);

my @returned_shell_programs;
my $expected_shell_programs_ref;

## Given an install set where two programs overlapp betwen bioconda and shell.
## No general shell preference is given,
## but one of the programs have been selected for shell installation.
my %parameter = (
    env_1 => {
        bioconda => {
            test_1 => 123,
            test_2 => 123,
            test_3 => 123,
        },
        shell => {
            test_2 => 234,
            test_3 => 234,
            test_4 => 234,

        },
    },
    prefer_shell  => 0,
    shell_install => [qw{ test_3 }],
);

## When the subroutine is executed
@returned_shell_programs = get_programs_for_shell_installation(
    {
        conda_programs_href        => $parameter{env_1}{bioconda},
        log                        => $log,
        prefer_shell               => $parameter{prefer_shell},
        shell_install_programs_ref => $parameter{shell_install},
        shell_programs_href        => $parameter{env_1}{shell},
    }
);
@returned_shell_programs = sort @returned_shell_programs;

## Then
$expected_shell_programs_ref = [qw{ test_3 test_4 }];
is_deeply( \@returned_shell_programs, $expected_shell_programs_ref,
    q{Testing shell_install} );

## Given a request to install a program via shell,
## for wich there are no shell installation recipe.
@parameter{shell_install} = [qw{ test_1 }];

## When the subroutine is executed
trap {
    @returned_shell_programs = get_programs_for_shell_installation(
        {
            conda_programs_href        => $parameter{env_1}{bioconda},
            log                        => $log,
            prefer_shell               => $parameter{prefer_shell},
            shell_install_programs_ref => $parameter{shell_install},
            shell_programs_href        => $parameter{env_1}{shell},
        }
      )
};

## Then, exit and print fatal message
ok( $trap->exit, q{Return exit signal} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message} );

## Given a general preference for shell installation
$parameter{prefer_shell}  = 1;
@parameter{shell_install} = [];

## When the subroutine is executed
@returned_shell_programs = get_programs_for_shell_installation(
    {
        conda_programs_href        => $parameter{env_1}{bioconda},
        log                        => $log,
        prefer_shell               => $parameter{prefer_shell},
        shell_install_programs_ref => $parameter{shell_install},
        shell_programs_href        => $parameter{env_1}{shell},
    }
);
@returned_shell_programs = sort @returned_shell_programs;

## Then
$expected_shell_programs_ref = [qw{ test_2 test_3 test_4 }];
is_deeply( \@returned_shell_programs, $expected_shell_programs_ref,
    q{Prefer shell} );

## Given no special requests
$parameter{prefer_shell}  = 0;
@parameter{shell_install} = [];

## When the subroutine is executed
@returned_shell_programs = get_programs_for_shell_installation(
    {
        conda_programs_href        => $parameter{env_1}{bioconda},
        log                        => $log,
        prefer_shell               => $parameter{prefer_shell},
        shell_install_programs_ref => $parameter{shell_install},
        shell_programs_href        => $parameter{env_1}{shell},
    }
);
@returned_shell_programs = sort @returned_shell_programs;

## Then
$expected_shell_programs_ref = [qw{ test_4 }];
is_deeply( \@returned_shell_programs, $expected_shell_programs_ref,
    q{No special requests} );

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
            store       => \$program_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help     Display this help message
    -v/--version  Display version
END_USAGE
}
