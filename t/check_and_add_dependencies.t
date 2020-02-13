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
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = '1.0.0';

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Check::Installation} => [qw{ check_and_add_dependencies }],
        q{MIP::Test::Fixtures}      => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Installation qw{ check_and_add_dependencies };

diag(   q{Test check_and_add_dependencies from Check::Installation.pm v}
      . $MIP::Check::Installation::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create temp logger
my $log = test_log( {} );

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
        dependency_href    => $installation{shell}{bio_prog_2}{conda_dependency},
        log                => $log,
        shell_program      => q{bio_prog_2},
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
            dependency_href    => $installation{shell}{bio_prog_2}{conda_dependency},
            log                => $log,
            shell_program      => q{bio_prog_2},
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
