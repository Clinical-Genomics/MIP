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
        q{MIP::Config}                 => [qw{ get_install_containers }],
        q{MIP::Environment::Container} => [qw{ set_executable_container_cmd }],
        q{MIP::Test::Fixtures}         => [qw{ test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Config qw{ get_install_containers };
use MIP::Environment::Container qw{ set_executable_container_cmd };

diag(   q{Test set_executable_container_cmd from Container.pm}
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
$active_parameter{container} =
  { get_install_containers( { container_config_file => $container_config_path, } ) };

## Given a container manager and some bind paths
my $container_manager = q{singularity};

## When parsing containers
my %container_cmd = set_executable_container_cmd(
    {
        active_parameter_href => \%active_parameter,
        container_href        => $active_parameter{container},
        container_manager     => $container_manager,
        parameter_href        => \%parameter,
    }
);

my $expected_arriba_cmd =
q{singularity exec --bind reference_dir!/a_dir:opt/conda/share/a_dir --cleanenv --no-home docker://docker.io/uhrigs/arriba:2.1.0 /arriba_v2.1.0/arriba};

## Then return command for how to execute arriba using singularity
is( $container_cmd{arriba}, $expected_arriba_cmd, q{Set singularity cmd for executable} );

## Given docker as container manager
$container_manager = q{docker};

## When parsing containers
%container_cmd = set_executable_container_cmd(
    {
        active_parameter_href => \%active_parameter,
        container_href        => $active_parameter{container},
        container_manager     => $container_manager,
        parameter_href        => \%parameter,
    }
);

my $expected_arriba_docker_cmd =
q{docker run --volume reference_dir!/a_dir:opt/conda/share/a_dir --rm --entrypoint "" docker://docker.io/uhrigs/arriba:2.1.0 /arriba_v2.1.0/arriba};

## Then return command for how to execute arriba using docker
is( $container_cmd{arriba}, $expected_arriba_docker_cmd, q{Set docker cmd for executable} );

## Given "no_executable_in_image" as executable value
my $no_executable_in_image = q{bwakit};

## When parsing containers
%container_cmd = set_executable_container_cmd(
    {
        active_parameter_href => \%active_parameter,
        container_href        => $active_parameter{container},
        container_manager     => $container_manager,
        parameter_href        => \%parameter,
    }
);

my $expected_no_executable_in_image =
  q{docker run --rm --entrypoint "" docker://docker.io/jemten/bwakit:0.7.17};
## Then set no command for executable
is(
    $container_cmd{$no_executable_in_image},
    $expected_no_executable_in_image,
    q{Set no cmd for "no_executable_in_image"}
);

## Given an executable with no executable_value
my $only_executable = q{run-bwamem};

my $expected_only_executable_cmd =
  q{docker run --rm --entrypoint "" docker://docker.io/jemten/bwakit:0.7.17 run-bwamem};

## Then return command for how to execute the executable
is( $container_cmd{$only_executable},
    $expected_only_executable_cmd, q{Set cmd for only executable} );

done_testing();
