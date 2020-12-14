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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Constants}      => [qw{ set_container_cmd }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Constants qw{ set_container_cmd };

diag(   q{Test set_container_cmd from Constants.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given an existing base command in CONTAINER_CMD
my $base_command = q{chanjo};
my $container_base_command =
  q{singularity exec docker.io/clinicalgenomics/chanjo:4.2.0 chanjo};
my %container_cmd = ( $base_command => $container_base_command, );

## When setting the global constant CONTAINER_CMD
set_container_cmd( { container_cmd_href => \%container_cmd, } );

## Then the chanjo base command should have a container command
is(
    $container_base_command,
    $MIP::Constants::CONTAINER_CMD{$base_command},
    q{CONTAINER_CMD is set}
);

done_testing();
