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
        q{MIP::Recipes::Install::Pip} => [qw{ install_pip_packages }],
        q{MIP::Test::Fixtures}        => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Install::Pip qw{ install_pip_packages };

diag(   q{Test install_pip_packages from Pip.pm v}
      . $MIP::Recipes::Install::Pip::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open my $FILEHANDLE, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

my $log = test_log( { no_screen => 1, } );

## Given no existing conda env
my $conda_env;
my %pip_packages = (
    genmod         => q{1.2.3},
    random_package => undef,
);

install_pip_packages(
    {
        conda_env         => $conda_env,
        FILEHANDLE        => $FILEHANDLE,
        pip_packages_href => \%pip_packages,
    }
);

## Close the filehandle
close $FILEHANDLE;

## Then install in conda main
my ($returned_command) =
  $file_content =~ /^(##\s+Install\s+PIP\s+packages\s+in\s+conda\s+main)/ms;

ok( $returned_command, q{Installed pip packages in conda main} );

## Given existing conda env

$conda_env = catdir( $Bin, qw{ data modules miniconda envs mip_travis } );

## Store file content in memory by using referenced variable
open $FILEHANDLE, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

install_pip_packages(
    {
        conda_env         => $conda_env,
        FILEHANDLE        => $FILEHANDLE,
        pip_packages_href => \%pip_packages,
    }
);

## Close the filehandle
close $FILEHANDLE;

# Then install in conda environment
($returned_command) =
  $file_content =~ /^(##\s+Install\s+PIP\s+packages\s+in\s+conda\s+environment)/ms;

ok( $returned_command, q{Installed pip packages in conda environment} );

## Given existing conda env when not following semantic versioning
$pip_packages{stranger} = q{not_semantic_versioning};

## Store file content in memory by using referenced variable
open $FILEHANDLE, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

trap {
    install_pip_packages(
        {
            conda_env         => $conda_env,
            FILEHANDLE        => $FILEHANDLE,
            pip_packages_href => \%pip_packages,
        }
    )
};

is( $trap->leaveby, q{die}, q{Exit if not semantic versioning} );
like(
    $trap->die,
    qr/The\s+version\s+number\s+does\s+not\s+match/xms,
    q{Throw error if not semantic versioning}
);

done_testing();
