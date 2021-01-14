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
use MIP::Test::Fixtures qw{ test_add_io_for_recipe test_log test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Analysis::Sv_combinevariantcallsets} =>
          [qw{ analysis_sv_combinevariantcallsets }],
        q{MIP::Test::Fixtures} => [qw{ test_add_io_for_recipe test_log test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Sv_combinevariantcallsets
  qw{ analysis_sv_combinevariantcallsets };

diag(   q{Test analysis_sv_combinevariantcallsets from Sv_combinevariantcallsets.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

## Given analysis parameters
my $recipe_name    = q{sv_combinevariantcallsets};
my $slurm_mock_cmd = catfile( $Bin, qw{ data modules slurm-mock.pl } );

my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
        recipe_name   => $recipe_name,
    }
);
$active_parameter{$recipe_name}                       = 1;
$active_parameter{tiddit}                             = 1;
$active_parameter{manta}                              = 1;
$active_parameter{delly_reformat}                     = 1;
$active_parameter{cnvnator_ar}                        = 1;
$active_parameter{recipe_core_number}{$recipe_name}   = 1;
$active_parameter{recipe_time}{$recipe_name}          = 1;
$active_parameter{sv_combinevariantcallsets_bcf_file} = 1;

my $case_id                    = $active_parameter{case_id};
my @structural_variant_callers = qw{ tiddit manta delly_reformat cnvnator_ar };

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
$parameter{cache}{structural_variant_callers} = [@structural_variant_callers];

my @order_recipes = ( @structural_variant_callers, $recipe_name );

SV_CALLER:
foreach my $sv_caller (@structural_variant_callers) {
    test_add_io_for_recipe(
    {
        file_info_href    => \%file_info,
        id                => $case_id,
        parameter_href    => \%parameter,
        recipe_name       => $sv_caller,
        step              => q{vcf},
    }
);
}

    test_add_io_for_recipe(
    {
        file_info_href    => \%file_info,
        id                => $case_id,
        order_recipes_ref => \@order_recipes,
        parameter_href    => \%parameter,
        recipe_name       => $recipe_name,
        step              => q{vcf},
    }
);

my %sample_info;

my $is_ok = analysis_sv_combinevariantcallsets(
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

## Given a single sample
@{ $active_parameter{sample_ids} } = $active_parameter{sample_ids}[0];

$is_ok = analysis_sv_combinevariantcallsets(
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
ok( $is_ok, q{ Executed analysis recipe } . $recipe_name . q{ single sample} );

done_testing();
