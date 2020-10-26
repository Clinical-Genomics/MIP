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
use MIP::Test::Fixtures qw{ test_log test_standard_cli };
use MIP::Update::Parameters qw{ update_with_dynamic_config_parameters };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Active_parameter} => [qw{ parse_vep_plugin }],
        q{MIP::Test::Fixtures}   => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ parse_vep_plugin };

diag(   q{Test parse_vep_plugin from Active_parameter.pm v}
      . $MIP::Active_parameter::VERSION
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

## Given a plugin hash
my %active_parameter = (
    vep_plugin => {
        dbNSFP => $rd_dna_config{vep_plugin}{dbNSFP},
    },
    vep_plugins_dir_path => $rd_dna_config{vep_plugins_dir_path},
);

#my %active_parameter_href = (
#    vep_plugin => {
#        dbNSFP => {
#            exists_check =>
#                q{directory},
#            path         => cwd(),
#            parameters   => [qw{ param_1 param_2 }],
#        },
#    },
#);
my @mip_vep_plugins = qw{ vep_plugin };

my $is_ok = parse_vep_plugin(
    {
        active_parameter_href => \%active_parameter,
        mip_vep_plugins_ref   => \@mip_vep_plugins,
    }
);

## Then return true
ok( $is_ok, q{Parsed vep plugin} );

done_testing();
