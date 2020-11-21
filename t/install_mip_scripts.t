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
use Path::Tiny;
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_constants test_log test_mip_hashes test_standard_cli };

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
        q{MIP::Recipes::Install::Mip_scripts} => [qw{ install_mip_scripts }],
        q{MIP::Test::Fixtures} =>
          [qw{ test_constants test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Install::Mip_scripts qw{ install_mip_scripts };

diag(   q{Test install_mip_scripts from Mip_scripts.pm v}
      . $MIP::Recipes::Install::Mip_scripts::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( { no_screen => 1, } );
test_constants( {} );

## Given install parameters
my %active_parameter =
  test_mip_hashes( { mip_hash_name => q{install_active_parameter}, } );
my $conda_prefix_path = catdir( $Bin, qw{data modules miniconda envs install_test} );

$active_parameter{conda_prefix_path} = $conda_prefix_path;

my $is_ok = install_mip_scripts(
    {
        active_parameter_href => \%active_parameter,
    }
);

## Then return TRUE
ok( $is_ok, q{Executed install mip_scripts recipe} );

## Given error in caching
my %process_return = (
    buffers_ref        => [],
    error_message => q{Error message},
    stderrs_ref        => [],
    stdouts_ref        => [],
    success            => 0,
);
test_constants( { test_process_return_href => \%process_return }, );
trap {
    install_mip_scripts(
        {
            active_parameter_href => \%active_parameter,
        }
    )
};

## Then exit and print error message
is( $trap->leaveby, q{die}, q{Error in case of install failure} );
like( $trap->die, qr/Error \s+ message /xms, q{Print error} );

path($conda_prefix_path)->remove_tree;

done_testing();
