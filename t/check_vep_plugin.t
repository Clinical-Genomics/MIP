#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
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
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Io::Read qw{ read_from_file };
use MIP::Test::Fixtures qw{ test_log };
use MIP::Update::Parameters qw{ update_with_dynamic_config_parameters };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Vep}            => [qw{ check_vep_plugin }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vep qw{ check_vep_plugin };

diag(   q{Test check_vep_plugin from Vep.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 1, } );

## Read config for most up to date format
my %rd_dna_config = read_from_file(
    {
        format => q{yaml},
        path   => catfile( dirname($Bin), qw{ templates mip_rd_dna_config.yaml } ),
    }
);
my %dynamic_parameter =
  ( cluster_constant_path => catfile( dirname($Bin), qw{ t data } ), );
update_with_dynamic_config_parameters(
    {
        active_parameter_href  => \%rd_dna_config,
        dynamic_parameter_href => \%dynamic_parameter,
    }
);
## Given an undefined vep plugin hash
my %vep_plugin;

my $is_not_ok = check_vep_plugin(
    {
        parameter_name       => q{vep_plugin},
        vep_plugin_href      => \%vep_plugin,
        vep_plugins_dir_path => $rd_dna_config{vep_plugins_dir_path},
    }
);
## Then there is nothing to check - return false
is( $is_not_ok, 0, q{Nothing to check} );

## Given a plugin hash
$vep_plugin{dbNSFP} = $rd_dna_config{vep_plugin}{dbNSFP};

my $is_ok = check_vep_plugin(
    {
        parameter_name       => q{vep_plugin},
        vep_plugin_href      => \%vep_plugin,
        vep_plugins_dir_path => $rd_dna_config{vep_plugins_dir_path},
    }
);

## Then return true
ok( $is_ok, q{Checked vep plugin hash} );

## Given a not valid hash ref
$vep_plugin{not_valid_annotation} = [q{not_a_valid_ref}];

trap {
    check_vep_plugin(
        {
            parameter_name       => q{vep_plugin},
            vep_plugin_href      => \%vep_plugin,
            vep_plugins_dir_path => $rd_dna_config{vep_plugins_dir_path},
        }
    )
};

## Then exit and throw FATAL log message
is( $trap->leaveby, q{die}, q{Exit if not a hash ref } );
like( $trap->die, qr/Is\s+not\s+a/xms, q{Not a hash ref} );

delete $vep_plugin{not_valid_annotation};

## Given a plugin that is missing
$vep_plugin{MockPlugin} = $rd_dna_config{vep_plugin}{ExACpLI};
trap {
    check_vep_plugin(
        {
            parameter_name       => q{vep_plugin},
            vep_plugin_href      => \%vep_plugin,
            vep_plugins_dir_path => $rd_dna_config{vep_plugins_dir_path},
        }
    )
};

## Then exit
is( $trap->leaveby, q{exit}, q{Exit if the plugin doesn't exist } );

done_testing();
