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
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Recipes::Install::Conda} => [qw{ install_conda_packages }],
        q{MIP::Test::Fixtures}          => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Install::Conda qw{ install_conda_packages };

diag(   q{Test install_conda_packages from Conda.pm v}
      . $MIP::Recipes::Install::Conda::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open my $FILEHANDLE, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## Given packages when no conda env
my %conda_packages = (
    picard => undef,
    pip    => q{19.1.0},
    python => q{2.7},
);
my $conda_env;
my $conda_env_path = catdir(qw{ a non existing conda env path });

trap {
    install_conda_packages(
        {
            conda_env           => $conda_env,
            conda_env_path      => $conda_env_path,
            conda_packages_href => \%conda_packages,
            FILEHANDLE          => $FILEHANDLE,
        }
    )
};

## Then broadcast installing in conda root message
like(
    $trap->stderr,
    qr/updating\s+packages\s+in\s+conda\s+root/xms,
    q{Install in conda root}
);

## Given packages when conda env is defined into a new env
$conda_env = q{test_env};

trap {
    install_conda_packages(
        {
            conda_env           => $conda_env,
            conda_env_path      => $conda_env_path,
            conda_packages_href => \%conda_packages,
            FILEHANDLE          => $FILEHANDLE,
        }
    )
};

## Then broadcast installing in conda root message
like(
    $trap->stderr,
    qr/Writing\s+installation\s+instructions\s+for\s+environmen/xms,
    q{Install in new conda env }
);

## Given packages when conda env is defined and env path exits
$conda_env      = q{test_env};
$conda_env_path = catdir( $Bin, qw{ data modules miniconda envs test_env } );

trap {
    install_conda_packages(
        {
            conda_env           => $conda_env,
            conda_env_path      => $conda_env_path,
            conda_packages_href => \%conda_packages,
            FILEHANDLE          => $FILEHANDLE,
        }
    )
};

## Then broadcast installing in existing conda environment
like(
    $trap->stderr,
    qr/install\s+packages\s+into\s+existing\s+environment/xms,
    q{Install in existing conda env }
);

## Given a package when not following semantic versioning
$conda_packages{picard} = q{not_semantic_version};

trap {
    install_conda_packages(
        {
            conda_env           => $conda_env,
            conda_env_path      => $conda_env_path,
            conda_packages_href => \%conda_packages,
            FILEHANDLE          => $FILEHANDLE,
        }
    )
};

## Then exit and throw FATAL log message
is( $trap->leaveby, q{die}, q{Exit if version is defined and not semantic version} );
like(
    $trap->die,
    qr/The\s+version\s+number\s+does\s+not\s+match/xms,
    q{Throw error if version is defined and not semantic version}
);

## Close the filehandle
close $FILEHANDLE;

done_testing();
