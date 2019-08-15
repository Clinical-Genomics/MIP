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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Get::Parameter} => [qw{ get_package_source_env_cmds }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::Parameter qw{ get_package_source_env_cmds };

diag(   q{Test get_package_source_env_cmds from Parameter.pm v}
      . $MIP::Get::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a program, when no prior to load cmd
my $program_name = q{genmod};

my %active_parameter = (
    load_env => {
        q{mip_pyv3.6} => {
            cnvnator      => q{prior_to_load_cmd;},
            $program_name => undef,
            method        => q{conda},
        },
    },
);

my @program_source_environment_cmds = get_package_source_env_cmds(
    {
        active_parameter_href => \%active_parameter,
        package_name          => $program_name,
    }
);

## Then return only load command
is_deeply(
    \@program_source_environment_cmds,
    [qw{ conda activate mip_pyv3.6 }],
    q{Got package source environment command}
);

## Given program, when using priors
my @program_source_environment_cmd_and_priors = get_package_source_env_cmds(
    {
        active_parameter_href => \%active_parameter,
        package_name          => q{cnvnator},
    }
);

## Then return load command with priors
is_deeply(
    \@program_source_environment_cmd_and_priors,
    [qw{ prior_to_load_cmd; conda activate mip_pyv3.6 }],
    q{Got package source environment command with priors}
);

done_testing();
