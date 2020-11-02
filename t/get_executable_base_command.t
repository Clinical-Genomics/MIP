#! /usr/bin/env perl

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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::Constants}               => [qw{ set_container_cmd }],
        q{MIP::Environment::Executable} => [qw{ get_executable_base_command }],
        q{MIP::Test::Fixtures}          => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Constants qw{ set_container_cmd };
use MIP::Environment::Executable qw{ get_executable_base_command };

diag(   q{Test get_executable_base_command from Executable.pm v}
      . $MIP::Environment::Executable::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a base command
my $base_command = q{chanjo};

## When the "base_command" do not exist in CONTAINER_CMD
my $returned_base_command =
  get_executable_base_command( { base_command => $base_command, } );

## Then the base command is returned unchanged
is( $base_command, $returned_base_command, q{Return original base command} );

## Given an existing base command in CONTAINER_CMD
my $container_base_command =
  q{singularity exec docker.io/clinicalgenomics/chanjo:4.2.0 chanjo};
my %container_cmd = ( $base_command => $container_base_command, );

set_container_cmd( { container_cmd_href => \%container_cmd, } );

## When the "base_command" do not exist in CONTAINER_CMD
$returned_base_command =
  get_executable_base_command( { base_command => $base_command, } );

## Then the base command is returned with container cmds
is( $returned_base_command, $container_base_command,
    q{Return base command with container cmds} );

done_testing();
