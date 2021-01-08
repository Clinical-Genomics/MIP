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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Active_parameter} => [qw{ get_package_env_cmds }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ get_package_env_cmds };

diag(   q{Test get_package_env_cmds from Active_parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a program
my $program_name = q{genmod};

## Given a load_env hash and a method
my %active_parameter = (
    load_env => {
        q{mip_pyv3.6} => {
            cnvnator      => undef,
            method        => q{conda},
            $program_name => undef,
        },
    },
);

## When no prior to load cmd
my @program_environment_cmds = get_package_env_cmds(
    {
        active_parameter_href => \%active_parameter,
        package_name          => $program_name,
    }
);

## Then return command to activate environment for package
is_deeply(
    \@program_environment_cmds,
    [qw{ conda activate mip_pyv3.6 }],
    q{Got package environment command}
);

done_testing();
