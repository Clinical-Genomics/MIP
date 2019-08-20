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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Check::Installation} => [qw{ check_python_compability }],
        q{MIP::Test::Fixtures}      => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Installation qw{ check_python_compability };

diag(   q{Test check_pyton_compability from Installation.pm v}
      . $MIP::Check::Installation::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

## Given nothing special
my %active_parameter = (
    python3_programs => [qw{ py3_prog_1 py3_prog_2 }],
    select_programs  => [],
    python_env       => {
        conda => {
            python       => 2.7,
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

trap {
    check_python_compability(
        {
            installation_set_href => $active_parameter{python_env},
            log                   => $log,
            python3_programs_ref  => $active_parameter{python3_programs},
            python_version        => $active_parameter{python_env}{conda}{python},
            select_programs_ref   => $active_parameter{select_programs},
        }
    )
};

## Then nothing special happens
ok( $trap->return, q{Returns on normal} );

## Given python not part of installation
delete $active_parameter{python_env}{conda}{python};

trap {
    check_python_compability(
        {
            installation_set_href => $active_parameter{python_env},
            log                   => $log,
            python3_programs_ref  => $active_parameter{python3_programs},
            python_version        => $active_parameter{python_env}{conda}{python},
            select_programs_ref   => $active_parameter{select_programs},
        }
    )
};

## Then print a warning to the log
like( $trap->stderr, qr/WARN/xms, q{Warn when pyhon is missing from installation} );

## Given a python version other than 2 or 3
foreach my $python_version (qw{ 4 hmm }) {
    trap {
        check_python_compability(
            {
                installation_set_href => $active_parameter{python_env},
                log                   => $log,
                python3_programs_ref  => $active_parameter{python3_programs},
                python_version        => $python_version,
                select_programs_ref   => $active_parameter{select_programs},
            }
        )
    };

    ## Then print a fatal message to the log and exit
    like( $trap->stderr, qr/FATAL/xms,
        q{Throw fatal message for python version: } . $python_version );
    ok( $trap->exit, q{Exit on python version: } . $python_version );
}

## Given a list of selected programs which are not compatible with python 2
@active_parameter{select_programs} = [qw{ py3_prog_1 py3_prog_2 }];
$active_parameter{python_env}{conda}{python} = q{2.7};

trap {
    check_python_compability(
        {
            installation_set_href => $active_parameter{python_env},
            log                   => $log,
            python3_programs_ref  => $active_parameter{python3_programs},
            python_version        => $active_parameter{python_env}{conda}{python},
            select_programs_ref   => $active_parameter{select_programs},
        }
    );
};

## Then print a fatal message to the log and exit
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal message for python incompability} );
ok( $trap->exit, q{Exit when python version is incompatible with programs} );

done_testing();
