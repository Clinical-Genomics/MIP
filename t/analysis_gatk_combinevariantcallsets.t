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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures
  qw{ test_add_io_for_recipe test_log test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.06;

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
        q{MIP::Recipes::Analysis::Gatk_combinevariantcallsets} =>
          [qw{ analysis_gatk_combinevariantcallsets }],
        q{MIP::Test::Fixtures} =>
          [qw{ test_add_io_for_recipe test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Gatk_combinevariantcallsets
  qw{ analysis_gatk_combinevariantcallsets };

diag(   q{Test analysis_gatk_combinevariantcallsets from Gatk_combinevariantcallsets.pm v}
      . $MIP::Recipes::Analysis::Gatk_combinevariantcallsets::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( { log_name => q{MIP}, no_screen => 1, } );

## Given analysis parameters and multiple callers
my $recipe_name    = q{gatk_combinevariantcallsets};
my $slurm_mock_cmd = catfile( $Bin, qw{ data modules slurm-mock.pl } );

my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
        recipe_name   => $recipe_name,
    }
);
$active_parameter{$recipe_name}                     = 1;
$active_parameter{gatk_variantrecalibration}        = 1;
$active_parameter{glnexus_merge}                    = 1;
$active_parameter{recipe_core_number}{$recipe_name} = 1;
$active_parameter{recipe_time}{$recipe_name}        = 1;
my $case_id = $active_parameter{case_id};
$active_parameter{gatk_path}                            = q{gatk.jar};
$active_parameter{gatk_combinevariantcallsets_bcf_file} = 1;

my @variant_callers = qw{ gatk_variantrecalibration deepvariant };

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

CALLER:
foreach my $caller (@variant_callers) {

    test_add_io_for_recipe(
        {
            file_info_href => \%file_info,
            id             => $case_id,
            parameter_href => \%parameter,
            recipe_name    => $caller,
            step           => q{vcf},
        }
    );
    push @{ $parameter{cache}{variant_callers} }, $caller;
}

my @order_recipes = ( qw{ gatk_variantrecalibration glnexus_merge }, $recipe_name );

test_add_io_for_recipe(
    {
        file_info_href    => \%file_info,
        id                => $case_id,
        parameter_href    => \%parameter,
        order_recipes_ref => \@order_recipes,
        recipe_name       => $recipe_name,
        step              => q{vcf},
    }
);

my %sample_info;

my $is_ok = analysis_gatk_combinevariantcallsets(
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
ok( $is_ok, q{ Executed analysis recipe } . $recipe_name . q{ multiple callers} );

## Given analysis parameters and single callers
@{ $parameter{cache}{variant_callers} } = $parameter{cache}{variant_callers}[0];
$is_ok = analysis_gatk_combinevariantcallsets(
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
ok( $is_ok, q{ Executed analysis recipe } . $recipe_name . q{ single caller} );

done_testing();
