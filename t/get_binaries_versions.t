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
use Modern::Perl qw{ 2018 };
use Readonly;
use Test::Trap qw{ :stderr:output(systemsafe) };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Environment::Executable} => [qw{ get_binaries_versions }],
        q{MIP::Test::Fixtures}          => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Executable qw{ get_binaries_versions };

diag(   q{Test get_binaries_versions from Executable.pm v}
      . $MIP::Environment::Executable::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given a binary and version when exists
my $binary      = q{mip};
my $binary_path = catfile( $Bin, qw{ data modules miniconda envs mip_travis bin mip } );
my %binary_info = ( $binary => $binary_path, );

my %binary_version = get_binaries_versions( { binary_info_href => \%binary_info, } );

my %expected_binary_version = (
    mip => {
        path    => $binary_path,
        version => q{v7.1.4},
    },
);

## Then mip binary name and version should be returned
is_deeply( \%binary_version, \%expected_binary_version, q{Got binary version} );

## Given a binary and version when it does not exists
$binary_info{just_to_enable_testing} = q{does_not_exist};

trap {
    get_binaries_versions( { binary_info_href => \%binary_info, } )
};

like( $trap->stderr, qr/Could\s+not\s+find\s+version/xms, q{Throw warning} );

done_testing();
