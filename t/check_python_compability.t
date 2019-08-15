#!/usr/bin/env perl

use 5.026;
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
use Modern::Perl qw{ 2018 };
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
Readonly my $TWO     => 2;

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
        say {*STDOUT} $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $VERSION . $NEWLINE;
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

use MIP::Check::Installation qw{ check_python_compability };

diag(   q{Test check_pyton3_compability from Installation.pm v}
      . $MIP::Check::Installation::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create temp logger
my $test_dir      = File::Temp->newdir();
my $test_log_path = catfile( $test_dir, q{test.log} );

## Creates log object
my $log = initiate_logger(
    {
        file_path => $test_log_path,
        log_name  => q{TEST},
    }
);

## Given nothing special
my %parameter = (
    python3_programs => [qw{ py3_prog_1 py3_prog_2 }],
    select_programs  => [],
    python_env       => {
        bioconda => {
            bio_prog_1 => 1,
            bio_prog_2 => 1,
        },
        conda => {
            conda_prog_1 => 1,
            conda_prog_2 => 1,
        },
        pip => {
            py_prog_1 => 1,
            py_prog_2 => 1,
        },
        shell => {
            bio_prog_1 => 1,
        },
    },
);

## When subroutine is executed
trap {
    check_python_compability(
        {
            installation_set_href => $parameter{python_env},
            python3_programs_ref  => $parameter{python3_programs},
            python_version        => $parameter{python_env}{conda}{python},
            select_programs_ref   => $parameter{select_programs},
            log                   => $log,
        }
    )
};

## Then nothing special happens
ok( $trap->return, q{Returns on normal} );

## Given python not part of installation
delete $parameter{python_env}{conda}{python};

## When subroutine is executed
trap {
    check_python_compability(
        {
            installation_set_href => $parameter{python_env},
            python3_programs_ref  => $parameter{python3_programs},
            python_version        => $parameter{python_env}{conda}{python},
            select_programs_ref   => $parameter{select_programs},
            log                   => $log,
        }
    )
};

## Then print a warning to the log
like( $trap->stderr, qr/WARN/xms, q{Warn when pyhon is missing from installation} );

## Given a python version other than 2 or 3
foreach my $python_version (qw{ 2,3 4 hmm }) {
    ## When subroutine is executed
    trap {
        check_python_compability(
            {
                installation_set_href => $parameter{python_env},
                python3_programs_ref  => $parameter{python3_programs},
                python_version        => $python_version,
                select_programs_ref   => $parameter{select_programs},
                log                   => $log,
            }
        )
    };

    ## Then print a fatal message to the log and exit
    like( $trap->stderr, qr/FATAL/xms,
        q{Throw fatal message for python version: } . $python_version );
    ok( $trap->exit, q{Exit on python version: } . $python_version );
}

## Given a list of selected programs which are not compatible with python 2
@parameter{select_programs} = [qw{ py3_prog_1 py3_prog_2 }];
$parameter{python_env}{conda}{python} = $TWO;

## When subroutine is executed
trap {
    check_python_compability(
        {
            installation_set_href => $parameter{python_env},
            python3_programs_ref  => $parameter{python3_programs},
            python_version        => $parameter{python_env}{conda}{python},
            select_programs_ref   => $parameter{select_programs},
            log                   => $log,
        }
    )
};

## Then print a fatal message to the log and exit
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal message for python incompability} );
ok( $trap->exit, q{Exit when python version is incompatible with programs} );

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
