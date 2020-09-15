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
use Modern::Perl qw{ 2018 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.04;

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
my $log = test_log( {} );

## Create starting hash
my $active_parameter_href = {
    select_programs => [],
    shell_install   => [],
    skip_programs   => [],
    conda           => {
        bio_prog_1   => 1,
        dependency_1 => 1,
        picard       => 1,
        python       => 2.7,
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
    },
    container => {
        multiqc => {
            executable => undef,
            uri        => q{docker://ewels/multiqc:v1.7},
        },
    },
};

## Copy starting hash to working copy
my $active_parameter_copy_href = clone($active_parameter_href);

## Given a parameter hash with conflicting options
$active_parameter_copy_href->{select_programs} = [qw{ bio_prog_1 }];
$active_parameter_copy_href->{skip_programs}   = [qw{ bio_prog_1 bio_prog_2 }];

## When subroutine is executed
trap {
    set_programs_for_installation(
        {
            active_parameter_href => $active_parameter_copy_href,
        }
    )
};

## Then print FATAL log message and exit
like( $trap->stderr, qr/mutually\sexclusive/xms, q{Fatal log message} );
ok( $trap->exit, q{Exit signal} );

## Given a parameter hash with a request to skip programs
$active_parameter_copy_href                    = clone($active_parameter_href);
$active_parameter_copy_href->{select_programs} = [];
$active_parameter_copy_href->{skip_programs}   = [qw{ py_prog_1 bio_prog_1 python }];

## When subroutine is executed
set_programs_for_installation(
    {
        active_parameter_href => $active_parameter_copy_href,
    }
);

## Then warn for no python and solve the installation as such
my $installation_href = clone($active_parameter_copy_href);
delete $installation_href->{conda}{qw{ bio_prog_1 python }};
delete $installation_href->{pip}{py_prog_1};
delete $installation_href->{shell}{qw{ bio_prog_1 }};

is_deeply( $active_parameter_copy_href, $installation_href, q{Solve installation} );

## Given a selective installation
$active_parameter_copy_href                    = clone($active_parameter_href);
$active_parameter_copy_href->{select_programs} = [qw{ python picard bio_prog_2 }];
$active_parameter_copy_href->{shell_install}   = [qw{ picard }];

## When subroutine is executed
set_programs_for_installation(
    {
        active_parameter_href => $active_parameter_copy_href,
    }
);

## Then solve the installation as such
$installation_href = clone($active_parameter_copy_href);
delete $installation_href->{conda}{qw{ bio_prog_1 picard }};
delete $installation_href->{pip}{qw{ py_prog_1 py_prog_2 }};
delete $installation_href->{shell}{bio_prog_1};
delete $installation_href->{container}{multiqc};

is_deeply( $active_parameter_copy_href, $installation_href, q{Solve installation} );

done_testing();

