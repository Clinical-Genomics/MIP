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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_add_io_for_recipe test_log test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.03;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $INDEL_TS_FILTER_LEVEL => 99.9;
Readonly my $RECIPE_CORE_NUMBER    => 16;
Readonly my $SNV_TS_FILTER_LEVEL   => 99.9;

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Analysis::Gatk_variantrecalibration} =>
          [qw{ analysis_gatk_variantrecalibration_wgs }],
        q{MIP::Test::Fixtures} => [qw{ test_add_io_for_recipe test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Gatk_variantrecalibration
  qw{ analysis_gatk_variantrecalibration_wgs };

diag(   q{Test analysis_gatk_variantrecalibration_wgs from Gatk_variantrecalibration.pm v}
      . $MIP::Recipes::Analysis::Gatk_variantrecalibration::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

## Given analysis parameters
my $recipe_name    = q{gatk_variantrecalibration_wgs};
my $slurm_mock_cmd = catfile( $Bin, qw{ data modules slurm-mock.pl } );

my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
        recipe_name   => $recipe_name,
    }
);
$active_parameter{$recipe_name}                     = 1;
$active_parameter{recipe_core_number}{$recipe_name} = $RECIPE_CORE_NUMBER;
$active_parameter{recipe_time}{$recipe_name}        = 1;
my $case_id   = $active_parameter{case_id};
my $sample_id = $active_parameter{sample_ids}[0];
@{ $active_parameter{gatk_variantrecalibration_annotations} } = (qw{ DP });
$active_parameter{gatk_variantrecalibration_resource_snv} = { resource => q{a_resource} };
$active_parameter{gatk_variantrecalibration_resource_indel} =
  { resource => q{a_resource} };
$active_parameter{gatk_variantrecalibration_snv_tsfilter_level} = $SNV_TS_FILTER_LEVEL;
$active_parameter{gatk_variantrecalibration_indel_tsfilter_level} =
  $INDEL_TS_FILTER_LEVEL;
$active_parameter{gatk_calculategenotypeposteriors}              = 1;
$active_parameter{gatk_variantrecalibration_snv_max_gaussians}   = 1;
$active_parameter{gatk_variantrecalibration_indel_max_gaussians} = 1;
$active_parameter{gatk_variantrecalibration_ts_tranches}         = [qw{ 100 99.9 }];

my %file_info = test_mip_hashes(
    {
        mip_hash_name => q{file_info},
        recipe_name   => $recipe_name,
    }
);

my %job_id;
my %parameter = test_mip_hashes(
    {
        mip_hash_name => q{recipe_parameter},
        recipe_name   => $recipe_name,
    }
);

test_add_io_for_recipe(
    {
        file_info_href    => \%file_info,
        id                => $case_id,
        parameter_href    => \%parameter,
        recipe_name       => $recipe_name,
        step              => q{vcf},
    }
);

$parameter{cache}{consensus_analysis_type} = q{wgs};
$parameter{cache}{trio}                    = 1;

my %sample_info = test_mip_hashes(
    {
        mip_hash_name => q{qc_sample_info},
        recipe_name   => $recipe_name,
    }
);

my $is_ok = analysis_gatk_variantrecalibration_wgs(
    {
        active_parameter_href => \%active_parameter,
        case_id               => $case_id,
        file_info_href        => \%file_info,
        job_id_href           => \%job_id,
        parameter_href        => \%parameter,
        profile_base_command  => $slurm_mock_cmd,
        recipe_name           => $recipe_name,
        sample_info_href      => \%sample_info,
    }
);

## Then return TRUE
ok( $is_ok, q{ Executed analysis recipe } . $recipe_name );

## Given a mixed consensus and single sample
$parameter{cache}{consensus_analysis_type} = q{mixed};
@{ $active_parameter{sample_ids} } = $sample_id;
$parameter{cache}{trio}                  = undef;
$sample_info{sample}{$sample_id}{mother} = 0;
$sample_info{sample}{$sample_id}{father} = 0;

$is_ok = analysis_gatk_variantrecalibration_wgs(
    {
        active_parameter_href => \%active_parameter,
        case_id               => $case_id,
        file_info_href        => \%file_info,
        job_id_href           => \%job_id,
        parameter_href        => \%parameter,
        profile_base_command  => $slurm_mock_cmd,
        recipe_name           => $recipe_name,
        sample_info_href      => \%sample_info,
    }
);

## Then return TRUE
ok( $is_ok, q{ Executed analysis recipe } . $recipe_name . q{ mixed and single sample} );

done_testing();
