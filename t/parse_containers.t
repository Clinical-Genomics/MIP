#! /usr/bin/env perl

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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Environment::Container} => [qw{ parse_containers }],
        q{MIP::Test::Fixtures}         => [qw{ test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Container qw{ parse_containers };

diag(   q{Test parse_containers from Container.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
    }
);

my %parameter = test_mip_hashes(
    {
        mip_hash_name => q{recipe_parameter},
    }
);

## Given an installation config
my $container_config_path = catfile( $Bin, qw{ data test_data miptest_container_config.yaml } );

## Given a container manager and an install config
$active_parameter{container_config_file} = $container_config_path;
$active_parameter{container_manager}     = q{singularity};

## When parsing containers
my $is_ok = parse_containers(
    {
        active_parameter_href => \%active_parameter,
        parameter_href        => \%parameter,
    }
);

## Then return true
is( $is_ok, 1, q{Parsed containers} );

done_testing();
