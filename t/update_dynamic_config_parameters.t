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
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::Update::Parameters} => [qw{ update_dynamic_config_parameters }],
        q{MIP::Test::Fixtures}     => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Update::Parameters qw{ update_dynamic_config_parameters };

diag(   q{Test update_dynamic_config_parameters from Update::Parameters.pm v}
      . $MIP::Update::Parameters::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %active_parameter = (
    cluster_constant_path  => catfile(qw{ root dir_1 dir_2 case_id! }),
    analysis_constant_path => q{analysis},
    case_id                => q{case_1},
    pedigree_file    => catfile(qw{ cluster_constant_path! case_id!_pedigree.yaml }),
    sample_info_file => catfile(
        qw{ cluster_constant_path! analysis_constant_path! case_id!_qc_sample_info.yaml }
    ),
    sv_vep_plugin => {
        ExACpLI => {
            exists_check => q{file},
            path => catfile(qw{ cluster_constant_path! analysis_constant_path! pli.txt }),
        },
        LofTool => {
            exists_check => q{file},
            path =>
              catfile(qw{ cluster_constant_path! analysis_constant_path! loftool.txt }),
        },
    },
);

## Given a cluster_constant_path when containing case_id!
my @dynamic_parameters = qw{ cluster_constant_path analysis_constant_path };

## Loop through all dynamic parameters and update info
PARAMETER:
foreach my $parameter_name (@dynamic_parameters) {

    ## Updates the active parameters to particular user/cluster for dynamic config parameters following specifications. Leaves other entries untouched.
    update_dynamic_config_parameters(
        {
            active_parameter_href  => \%active_parameter,
            dynamic_parameter_href => { case_id => $active_parameter{case_id}, },
            parameter_name         => $parameter_name,
        }
    );
}

## Then cluster constant path should be updated with supplied case id
my $updated_cluster_constant_path = catdir(qw{ root dir_1 dir_2 case_1 });
is(
    $active_parameter{cluster_constant_path},
    $updated_cluster_constant_path,
    q{Updated cluster constant path with case_id}
);

## Given parameters with dynamic parameters
my @order_parameters = qw{ pedigree_file sample_info_file sv_vep_plugin };

my %dynamic_parameter = (
    cluster_constant_path  => $active_parameter{cluster_constant_path},
    analysis_constant_path => $active_parameter{analysis_constant_path},
    case_id                => $active_parameter{case_id},
);

## Loop through all parameters and update info
PARAMETER:
foreach my $parameter_name (@order_parameters) {

    ## Updates the active parameters to particular user/cluster for dynamic config parameters following specifications. Leaves other entries untouched.
    update_dynamic_config_parameters(
        {
            active_parameter_href  => \%active_parameter,
            dynamic_parameter_href => \%dynamic_parameter,
            parameter_name         => $parameter_name,
        }
    );
}

## Then all parameters containing dynamic parameters should have been updated
my $updated_pedigree_file = catfile(qw{ root dir_1 dir_2 case_1 case_1_pedigree.yaml });
is( $active_parameter{pedigree_file},
    $updated_pedigree_file, q{Updated pedigree file path} );

my $updated_sample_info_file =
  catfile(qw{ root dir_1 dir_2 case_1 analysis case_1_qc_sample_info.yaml });
is( $active_parameter{sample_info_file},
    $updated_sample_info_file, q{Updated sample_info_file path} );

my $updated_sv_vep_pli = catfile(qw{ root dir_1 dir_2 case_1 analysis pli.txt });
is( $active_parameter{sv_vep_plugin}{ExACpLI}{path},
    $updated_sv_vep_pli, q{Updated sv_vep_plugin pli path} );

my $updated_sv_vep_loftool = catfile(qw{ root dir_1 dir_2 case_1 analysis loftool.txt });
is( $active_parameter{sv_vep_plugin}{LofTool}{path},
    $updated_sv_vep_loftool, q{Updated sv_vep_plugin loftool path} );

done_testing();
