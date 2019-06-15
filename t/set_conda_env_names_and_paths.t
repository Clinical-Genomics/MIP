#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
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
use MIP::Constants qw{ $COMMA $NEWLINE $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli test_log };

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
        q{MIP::Set::Parameter} => [qw{ set_conda_env_names_and_paths }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::Parameter qw{ set_conda_env_names_and_paths };

diag(   q{Test set_conda_env_names_and_paths from Parameter.pm v}
      . $MIP::Set::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create log object
my $log = test_log( {} );

## Given a test hash with a defined name for MIP main env (emip) and two sub
## environments whereof one is named.
my %parameter = (
    conda_dir_path   => catdir( $Bin, qw{ data modules miniconda } ),
    environment_name => {
        emip        => q{Test},
        etest_env_1 => undef,
        etest_env_2 => q{test_2},
    },
    installations         => [qw{ emip etest_env_1 etest_env_2 }],
    environment_base_name => q{Test_base},
    no_environment_date   => 1,
);

## When the subroutine is executed
set_conda_env_names_and_paths(
    {
        log            => $log,
        parameter_href => \%parameter,
    }
);
## Then:
## The unnamed sub environment shall be named with the main environment name
## followed by the sub environment
is( $parameter{environment_name}{etest_env_1},
    q{Test_test_env_1}, q{Set undefined environment name} );
## The named sub environment shall retain it's name
is( $parameter{environment_name}{etest_env_2},
    q{test_2}, q{Leave defined environment name} );
## The path to MIP's conda env shall be set
is(
    $parameter{emip}{conda_prefix_path},
    catdir( $parameter{conda_dir_path}, qw{ envs Test } ),
    q{Set main environment specific conda paths}
);
## The path to the MIP's sub environment shall be set
is(
    $parameter{etest_env_2}{conda_prefix_path},
    catdir( $parameter{conda_dir_path}, qw{ envs test_2 } ),
    q{Set sub environment specific conda paths}
);

## Given that MIP's main and sub environments are unnamed
$parameter{environment_name}{emip}        = undef;
$parameter{environment_name}{etest_env_1} = undef;

## When the subroutine is executed
trap {
    set_conda_env_names_and_paths(
        {
            log            => $log,
            parameter_href => \%parameter,
        }
      )
};

## Then:
## The name of base env plus sub environemnt shall be the environemnt
is( $parameter{environment_name}{etest_env_1},
    q{Test_base_test_env_1},
    q{Set undefined environment name when emip env name is undefined} );
## The conda path for MIP's main environment shall be MIP's base environment
is(
    $parameter{emip}{conda_prefix_path},
    catdir( $parameter{conda_dir_path}, qw{ envs Test_base } ),
    q{Set main environment to base }
);
## The Log shall warn
like( $trap->stderr, qr/WARN/xms, q{Warn when using default base name} );

## Given to prefix and suffix to environemnt name
$parameter{environment_prefix} = q{D};
$parameter{environment_suffix} = q{JD};

set_conda_env_names_and_paths(
    {
        log            => $log,
        parameter_href => \%parameter,
    }
);

## Then return formatted environemnt name
is( $parameter{environment_name}{etest_env_2},
    q{D_test_2_JD}, q{Preppend prefix and append suffix} );

## Given request to add environemnt creation date
%parameter = (
    conda_dir_path   => catdir( $Bin, qw{ data modules miniconda } ),
    environment_name => {
        emip => q{Test},
    },
    installations         => [qw{ emip }],
    environment_base_name => q{Test_base},
    no_environment_date   => 0,
);

set_conda_env_names_and_paths(
    {
        log            => $log,
        parameter_href => \%parameter,
    }
);

## Then add six figure date
like( $parameter{environment_name}{emip}, qr/_\d{6}/xms, q{Add date} );

done_testing();
