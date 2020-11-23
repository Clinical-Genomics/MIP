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
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.04;

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
        q{MIP::Get::Parameter} => [qw{ get_install_parameter_attribute }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::Parameter qw{ get_install_parameter_attribute };

diag(   q{Test get_install_parameter_attribute from Parameter.pm v}
      . $MIP::Get::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given parameters with scalar, array and hashes
my %parameter = test_mip_hashes( { mip_hash_name => q{install_active_parameter}, } );

## Hash attribute
my %singularity_container = get_install_parameter_attribute(
    {
        parameter_href => \%parameter,
        parameter_name => q{container},
    }
);

## Then all singularity containers should have been loaded
is_deeply( \%singularity_container, \%{ $parameter{container} }, q{Got hash attribute } );

## Array attribute
my @vep_plugins = get_install_parameter_attribute(
    {
        parameter_href => \%parameter,
        parameter_name => q{vep_plugins},
    }
);

## Then all VEP plugins should have been loaded
is_deeply( \@vep_plugins, \@{ $parameter{vep_plugins} }, q{Got array attibute} );

## Scalar attribute
my $bash_set_errexit = get_install_parameter_attribute(
    {
        parameter_href => \%parameter,
        parameter_name => q{bash_set_errexit},
    }
);

## Then scalar attribute value should be returned
is( 1, $bash_set_errexit, q{Got scalar attibute} );

## Given non-existing key
my $does_not_exist = q{this key does not exist};
my $is_ok          = trap {
    get_install_parameter_attribute(
        {
            parameter_href => \%parameter,
            parameter_name => $does_not_exist,
        }
    )
};

## Then exit and throw FATAL log message
is( $trap->leaveby, q{die}, q{Exit if the key cannot be found} );
like( $trap->die, qr/Could\s+not/xms, q{Throw error} );

done_testing();
