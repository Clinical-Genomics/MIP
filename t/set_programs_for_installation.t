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
use Clone qw{ clone };
use Modern::Perl qw{ 2014 };
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
        q{MIP::Set::Parameter} => [qw{ set_programs_for_installation }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::Parameter qw{ set_programs_for_installation };

diag(   q{Test set_programs_for_installation from Set::Parameter.pm v}
      . $MIP::Set::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create temp logger
my $log = test_log();

## Create starting hash
my $parameter_href = {
    installations    => [qw{ test_env }],
    python3_programs => [],
    select_programs  => [],
    shell_install    => [],
    skip_programs    => [],
    test_env         => {
        conda => {
            bio_prog_1   => 1,
            dependency_1 => 1,
            python       => 2.7,
            snpeff       => 1,
            snpsift      => 1,
        },
        pip => {
            py_prog_1 => 1,
            py_prog_2 => 1,
        },
        shell => {
            bio_prog_1 => {
                conda_dependency => {
                    dependency_1 => 1,
                },
                version => 1,
            },
            bio_prog_2 => {
                conda_dependency => {
                    dependency_1 => 1,
                },
                version => 2,
            },
            snpeff => {
                snpeff_genome_versions => [qw{ GRCh37 }],
                version                => 2,
            },
        },
    },
};

## Copy starting hash to working copy
my $parameter_copy_href = clone($parameter_href);

## Given a parameter hash with conflicting options
$parameter_copy_href->{select_programs} = [qw{ bio_prog_1 }];
$parameter_copy_href->{skip_programs}   = [qw{ bio_prog_1 bio_prog_2 }];

## When subroutine is executed
trap {
    set_programs_for_installation(
        {
            installation   => q{test_env},
            parameter_href => $parameter_copy_href,
            log            => $log,
        }
      )
};

## Then print FATAL log message and exit
like( $trap->stderr, qr/FATAL/xms, q{Fatal log message} );
ok( $trap->exit, q{Exit signal} );

## Given a parameter hash with the select_program option and installation of several environements
$parameter_copy_href                    = clone($parameter_href);
$parameter_copy_href->{select_programs} = [qw{ bio_prog_1 }];
$parameter_copy_href->{installations}   = [qw{ test_env test_env2 }];

## When subroutine is executed
trap {
    set_programs_for_installation(
        {
            installation   => q{test_env},
            parameter_href => $parameter_copy_href,
            log            => $log,
        }
      )
};

## Then print FATAL log message for trying to use --select_program together with and installation of more than one environment
like( $trap->stderr, qr/FATAL/xms, q{Fatal log message} );
ok( $trap->exit, q{Exit signal} );

## Given a parameter hash with a request to skip programs
$parameter_copy_href = clone($parameter_href);
$parameter_copy_href->{skip_programs} = [qw{ py_prog_1 bio_prog_1 python }];

## When subroutine is executed
trap {
    set_programs_for_installation(
        {
            installation   => q{test_env},
            parameter_href => $parameter_copy_href,
            log            => $log,
        }
      )
};

##Then warn for no python and solve the installation as such
my $installation_href = clone($parameter_copy_href);
$installation_href->{test_env}{snpeff_genome_versions} = [qw{ GRCh37 }];
delete $installation_href->{test_env}{conda}{qw{ bio_prog_1 python }};
delete $installation_href->{test_env}{pip}{py_prog_1};
delete $installation_href->{test_env}{shell}{qw{ bio_prog_1 snpeff }};

like( $trap->stderr, qr/WARN/xms, q{Warn for no python} );
is_deeply( $parameter_copy_href, $installation_href, q{Solve installation} );

## Given a selective installation
$parameter_copy_href                    = clone($parameter_href);
$parameter_copy_href->{select_programs} = [qw{ python snpeff bio_prog_2 }];
$parameter_copy_href->{shell_install}   = [qw{ snpeff }];

## When subroutine is executed
set_programs_for_installation(
    {
        installation   => q{test_env},
        parameter_href => $parameter_copy_href,
        log            => $log,
    }
);

## Then solve the installation as such
$installation_href = clone($parameter_copy_href);
$installation_href->{test_env}{snpeff_genome_versions} = [qw{ GRCh37 }];
delete $installation_href->{test_env}{conda}{qw{ bio_prog_1 snpeff snpsift }};
delete $installation_href->{test_env}{pip}{qw{ py_prog_1 py_prog_2 }};
delete $installation_href->{test_env}{shell}{bio_prog_1};

is_deeply( $parameter_copy_href, $installation_href, q{Solve installation} );

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
