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
our $VERSION = 1.05;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

BEGIN {

    use MIP::Test::Fixtures qw{ test_import test_mip_hashes };

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

my $log = test_log( {} );

my %active_parameter =
  test_mip_hashes( { mip_hash_name => q{install_active_parameter}, } );

## Copy starting hash to working copy
my %active_parameter_copy = %{ clone( \%active_parameter ) };

## Given a parameter hash with conflicting options
$active_parameter_copy{select_programs} = [qw{ bwa }];
$active_parameter_copy{skip_programs}   = [qw{ htslib }];

## When subroutine is executed
trap {
    set_programs_for_installation(
        {
            active_parameter_href => \%active_parameter_copy,
        }
    )
};

## Then print FATAL log message and exit
like( $trap->stderr, qr/mutually\sexclusive/xms, q{Fatal log message} );
ok( $trap->exit, q{Exit signal} );

## Given a parameter hash with a request to skip programs
%active_parameter_copy                  = %{ clone( \%active_parameter ) };
$active_parameter_copy{select_programs} = [];
$active_parameter_copy{skip_programs}   = [qw{ htslib }];

## When subroutine is executed
set_programs_for_installation(
    {
        active_parameter_href => \%active_parameter_copy,
    }
);

## Then solve the installation as such
my %expected_container = %{ clone( $active_parameter{container} ) };
delete $expected_container{htslib};

is_deeply( $active_parameter_copy{container},
    \%expected_container, q{Solve installation} );

## Given a selective installation
%active_parameter_copy = %{ clone( \%active_parameter ) };
$active_parameter_copy{select_programs} = [qw{ bwa }];

## When subroutine is executed
set_programs_for_installation(
    {
        active_parameter_href => \%active_parameter_copy,
    }
);

## Then solve the installation as such
%expected_container = ();
$expected_container{bwa} = $active_parameter{container}{bwa};

is_deeply( $active_parameter_copy{container},
    \%expected_container, q{Solve installation} );

done_testing();

