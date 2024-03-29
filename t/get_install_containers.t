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

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Config} => [qw{ get_install_containers }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Config qw{ get_install_containers };

diag(   q{Test get_install_containers from Config.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given an installation config
my $container_config_path = catfile( $Bin, qw{ data test_data miptest_container_config.yaml } );

## When getting containers from install config file
my %container = get_install_containers( { container_config_file => $container_config_path, } );

## Then return container info from installation config
my %expected_container = (
    arriba => {
        bind_path => {
            arriba => q{reference_dir!/a_dir:opt/conda/share/a_dir}
        },
        executable => {
            arriba            => q{/arriba_v2.1.0/arriba},
            q{draw_fusions.R} => q{/arriba_v2.1.0/draw_fusions.R},
        },
        uri => q{docker.io/uhrigs/arriba:2.1.0},
    },
);

is_deeply( $container{arriba}, $expected_container{arriba}, q{Got install containers} );

done_testing();
