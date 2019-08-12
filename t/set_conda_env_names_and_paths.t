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
use Modern::Perl qw{ 2014 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $NEWLINE $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli test_log };

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
        q{MIP::Set::Parameter} => [qw{ set_conda_env_names_and_paths }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::Parameter qw{ set_conda_env_names_and_paths };

diag(   q{Test set_conda_env_names_and_paths from Parameter.pm v}
      . $MIP::Set::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a test hash with an undefined environment name and date addition request.
my %active_parameter = (
    conda_path       => catdir( $Bin, qw{ data modules miniconda } ),
    environment_name => {
        etest_env => undef,
    },
    installations         => [qw{ etest_env }],
    environment_base_name => q{base},
    add_environment_date  => 1,
);

set_conda_env_names_and_paths(
    {
        active_parameter_href => \%active_parameter,
    }
);
## Then construct and set conda environment name
like( $active_parameter{environment_name}{etest_env},
    qr/base_test_env_\d{6}/xms, q{Set environment name} );

is(
    $active_parameter{etest_env}{conda_prefix_path},
    catdir(
        $active_parameter{conda_path}, q{envs},
        $active_parameter{environment_name}{etest_env}
    ),
    q{Set environment conda path}
);

done_testing();
