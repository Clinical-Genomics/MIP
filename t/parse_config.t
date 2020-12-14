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
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Config}         => [qw{ parse_config }],
        q{MIP::Io::Read}       => [qw{ read_from_file }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Config qw{ parse_config };
use MIP::Io::Read qw{ read_from_file };

diag(   q{Test parse_config from Config.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given
my @definition_file_paths = (
    catfile( dirname($Bin), qw{ definitions mip_parameters.yaml } ),
    catfile( dirname($Bin), qw{ definitions analyse_parameters.yaml } ),
    catfile( dirname($Bin), qw{ definitions rd_dna_parameters.yaml } ),
);
## Loads a YAML file into an arbitrary hash and returns it.
my %parameter;
DEFINITION_FILE:
foreach my $definition_file_path (@definition_file_paths) {

    %parameter = (
        %parameter,
        read_from_file(
            {
                format => q{yaml},
                path   => $definition_file_path,
            }
        ),
    );
}
my %active_parameter =
  ( config_file => catfile( dirname($Bin), qw{ templates mip_rd_dna_config.yaml } ), );
$active_parameter{case_id}                = q{case_1};
$active_parameter{cluster_constant_path}  = q{cluster_constant_path};
$active_parameter{analysis_constant_path} = q{analysis_constant_path};

my $is_ok = parse_config(
    {
        active_parameter_href => \%active_parameter,
        parameter_href        => \%parameter,
    }
);

## Then
ok( $is_ok, q{Parsed config} );

done_testing();
