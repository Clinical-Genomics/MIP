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
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2014 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
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
    my @modules = (q{MIP::Check::Installation});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Check::Installation qw{ check_and_add_dependencies };
use MIP::Test::Fixtures qw{ test_log };

diag(   q{Test check_and_add_dependencies from Check::Installation.pm v}
      . $MIP::Check::Installation::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create temp logger
my $log = test_log();

## Given a parameter hash specifying an installation with no conflicting dependencies
my %installation = (
    conda => {
        bio_prog_1   => 1,
        dependency_1 => 1,
    },
    shell => {
        bio_prog_2 => {
            conda_dependency => {
                dependency_1 => 1,
                dependency_2 => 2,
            },
            version => 1,
        },
    },
);

## When subroutine is executed
check_and_add_dependencies(
    {
        conda_program_href => $installation{conda},
        dependency_href => $installation{shell}{bio_prog_2}{conda_dependency},
        log             => $log,
        shell_program   => q{bio_prog_2},
    }
);

## Then the shell dependency is added to the conda installation
my %dependency_added = (
    conda => {
        bio_prog_1   => 1,
        dependency_1 => 1,
        dependency_2 => 2,
    },
    shell => {
        bio_prog_2 => {
            conda_dependency => {
                dependency_1 => 1,
                dependency_2 => 2,
            },
            version => 1,
        },
    },
);
is_deeply( \%installation, \%dependency_added, q{Add dependencies} );

## Given a parameter hash specifying an installation with conflicting dependencies
$installation{conda}{dependency_1} = 2;

## When subroutine is executed
trap {
    check_and_add_dependencies(
        {
            conda_program_href => $installation{conda},
            dependency_href =>
              $installation{shell}{bio_prog_2}{conda_dependency},
            log           => $log,
            shell_program => q{bio_prog_2},
        }
      )
};

## Then print FATAL log message and exit
like( $trap->stderr, qr/FATAL/xms, q{Print fatal log message} );
ok( $trap->exit, q{Exit signal} );

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
