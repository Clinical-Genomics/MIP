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
        q{MIP::Update::Parameters} => [qw{ update_with_dynamic_config_parameters }],
        q{MIP::Test::Fixtures}     => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Update::Parameters qw{ update_with_dynamic_config_parameters };

diag(   q{Test update_with_dynamic_config_parameters from Update::Parameters.pm v}
      . $MIP::Update::Parameters::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a data structure with dynamic parameters (!) and a hash map over those
my %active_parameter = (
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
    some_array => [
        catfile(qw{ cluster_constant_path! analysis_constant_path! loftool.txt }),
        undef,
        { a_hash       => catfile(q{cluster_constant_path!}) },
        { another_hash => [ catfile(q{cluster_constant_path!}), ] },
    ],
);

my %dynamic_parameter = (
    cluster_constant_path  => catdir(qw{ root dir_1 dir_2 case_1 }),
    analysis_constant_path => q{analysis},
    case_id                => q{case_1},
);

## Then updates the active parameters with dynamic config parameters following specifications.
update_with_dynamic_config_parameters(
    {
        active_parameter_href  => \%active_parameter,
        dynamic_parameter_href => \%dynamic_parameter,
    }
);

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

my $some_array_ref = [
    catfile(qw{ root dir_1 dir_2 case_1 analysis loftool.txt }),
    undef,
    { a_hash       => catfile(qw{ root dir_1 dir_2 case_1 }) },
    { another_hash => [ catfile(qw{ root dir_1 dir_2 case_1 }) ] }
];
is_deeply( $active_parameter{some_array}, $some_array_ref, q{Updated array reference} );
done_testing();
