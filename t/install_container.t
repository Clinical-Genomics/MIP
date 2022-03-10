#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp qw{ tempdir };
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
use MIP::Test::Fixtures qw{ test_constants test_log test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Install::Container} => [qw{ install_containers }],
        q{MIP::Test::Fixtures}              => [qw{ test_constants test_log test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Install::Container qw{ install_containers };

diag(   q{Test install_containers from Container.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $test_dir = File::Temp->newdir();

test_log( { no_screen => 1, } );
test_constants( {} );

## Given install parameters
my %active_parameter = test_mip_hashes( { mip_hash_name => q{install_active_parameter}, } );
$active_parameter{reference_dir}          = catdir( $test_dir, qw{ a dir } );
$active_parameter{container_manager}      = q{docker};
$active_parameter{conda_environment_path} = tempdir( CLEANUP => 1 );

my $is_ok = install_containers(
    {
        active_parameter_href => \%active_parameter,
        container_href        => $active_parameter{container},
    }
);

## Then return TRUE
ok( $is_ok, q{Executed install container recipe} );

## Given error in caching
my %process_return = (
    buffers_ref   => [],
    error_message => q{Error message},
    stderrs_ref   => [],
    stdouts_ref   => [],
    success       => 0,
);
test_constants( { test_process_return_href => \%process_return }, );
trap {
    install_containers(
        {
            active_parameter_href => \%active_parameter,
            container_href        => $active_parameter{container},
        }
    )
};

## Then exit and print error message
is( $trap->leaveby, q{die}, q{Error in case of caching failure} );
like( $trap->die, qr/Error \s+ message /xms, q{Print error} );

done_testing();
