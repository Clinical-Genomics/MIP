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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import test_mip_hashes };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Active_parameter} => [qw{ set_programs_for_installation }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ set_programs_for_installation };

diag(   q{Test set_programs_for_installation from Active_parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( {} );

my %active_parameter = test_mip_hashes( { mip_hash_name => q{install_active_parameter}, } );

## Given neither select_programs or skip_programs input
## Copy to local copy
my %active_parameter_copy  = %{ clone( \%active_parameter ) };
my %expected_all_container = %{ $active_parameter{container} };

## When setting programs to install
set_programs_for_installation(
    {
        active_parameter_href => \%active_parameter_copy,
    }
);

## Then all programs in container should be installed
is_deeply( $active_parameter_copy{container},
    \%expected_all_container, q{Solve installation for pipeline} );

## Regenerate local copy
%active_parameter_copy = %{ clone( \%active_parameter ) };

## Given a parameter hash with conflicting options
$active_parameter_copy{select_programs} = [qw{ bwa }];
$active_parameter_copy{skip_programs}   = [qw{ htslib }];

## When setting programs to install
trap {
    set_programs_for_installation(
        {
            active_parameter_href => \%active_parameter_copy,
        }
    )
};

## Then print FATAL log message and exit
like(
    $trap->stderr,
    qr/mutually \s+ exclusive/xms,
    q{Thorw fatal log message if conflicting options}
);
ok( $trap->exit, q{Exit if conflicting options} );

## Given an active_parameter hash with a request to skip programs
%active_parameter_copy                  = %{ clone( \%active_parameter ) };
$active_parameter_copy{select_programs} = [];
$active_parameter_copy{skip_programs}   = [qw{ htslib }];

## When setting programs to install
set_programs_for_installation(
    {
        active_parameter_href => \%active_parameter_copy,
    }
);

## Then htslib should be deleted from active_parameters
my %expected_container = %{ clone( $active_parameter{container} ) };
delete $expected_container{htslib};

is_deeply( $active_parameter_copy{container},
    \%expected_container, q{Solve installation when skipping programs} );

## Given a request to install a specific program
%active_parameter_copy = %{ clone( \%active_parameter ) };
$active_parameter_copy{select_programs} = [qw{ bwa }];

## When setting programs to install
set_programs_for_installation(
    {
        active_parameter_href => \%active_parameter_copy,
    }
);

## Then only the selected program should be present in the active_parameters hash
%expected_container = ();
$expected_container{bwa} = $active_parameter{container}{bwa};

is_deeply( $active_parameter_copy{container},
    \%expected_container, q{Solve installation when selecting programs} );

done_testing();

