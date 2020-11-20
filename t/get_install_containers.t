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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

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
        q{MIP::Config}         => [qw{ get_install_containers }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Config qw{ get_install_containers };

diag(   q{Test get_install_containers from Config.pm v}
      . $MIP::Config::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given an installation config
my $install_config_path =
  catfile( $Bin, qw{ data test_data install_active_parameters.yaml } );

## When getting containers from install config file
my %container =
  get_install_containers( { install_config_file => $install_config_path, } );

## Then return container info from installation config
my %expected_container = (
    arriba => {
        bind_path => {
            arriba => q{reference_dir!/a_dir:opt/conda/share/a_dir}
        },
        executable => {
            arriba            => q{/arriba_v1.2.0/arriba},
            q{draw_fusions.R} => q{/arriba_v1.2.0/draw_fusions.R},
        },
        uri => q{docker.io/uhrigs/arriba:1.2.0},
    },
);

is_deeply( $container{arriba}, $expected_container{arriba}, q{Got install containers} );

done_testing();
