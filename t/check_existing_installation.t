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
use MIP::Constants qw{$COLON $COMMA $SPACE };
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
        q{MIP::Check::Installation} => [qw{ check_existing_installation }],
        q{MIP::Test::Fixtures}      => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Installation qw{ check_existing_installation };

diag(   q{Test check_existing_installation from Installation.pm v}
      . $MIP::Check::Installation::VERSION
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

my $conda_prefix_path = catdir( $Bin, qw{ data modules miniconda } );
my $log               = test_log( {} );
my $program_directory_path =
  catdir( $conda_prefix_path, qw{ envs test_env share picard } );

## Given existing installation
trap {
    check_existing_installation(
        {
            conda_environment      => q{test_env},
            conda_prefix_path      => $conda_prefix_path,
            FILEHANDLE             => $FILEHANDLE,
            log                    => $log,
            program_directory_path => $program_directory_path,
            program_name           => q{Picard},
        }
    );
};
close $FILEHANDLE;
## Then overwrite
like( $trap->stderr, qr/WARN/xms, q{Warn when installation exists} );
like( $file_content, qr/Removing\sold\sPicard/xms, q{Write instructions to file} );

## Store file content in memory by using referenced variable
open $FILEHANDLE, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

$program_directory_path =
  catdir( $conda_prefix_path, qw{ envs test_env share dummy_program } );

## Given non-existing installation
my $is_installed = check_existing_installation(
    {
        conda_environment      => q{test_env},
        conda_prefix_path      => $conda_prefix_path,
        FILEHANDLE             => $FILEHANDLE,
        log                    => $log,
        program_directory_path => $program_directory_path,
        program_name           => q{dummy_program},
    }
);
close $FILEHANDLE;

## Then return 0
is( $is_installed, 0, q{Return 0 when not installed} );

done_testing();
