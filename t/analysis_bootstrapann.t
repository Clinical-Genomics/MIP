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
  qw{ test_add_io_for_recipe test_log test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Analysis::BootstrapAnn} => [qw{ analysis_bootstrapann }],
        q{MIP::Test::Fixtures} =>
          [qw{ test_add_io_for_recipe test_log test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::BootstrapAnn qw{ analysis_bootstrapann };

diag(   q{Test analysis_bootstrapann from BootstrapAnn.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

## Given build parameters
my $recipe_name    = q{bootstrapann};
my $slurm_mock_cmd = catfile( $Bin, qw{ data modules slurm-mock.pl } );

my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
        recipe_name   => $recipe_name,
    }
);
$active_parameter{$recipe_name}                     = 1;
$active_parameter{recipe_core_number}{$recipe_name} = 1;
$active_parameter{recipe_time}{$recipe_name}        = 1;
my $sample_id = $active_parameter{sample_ids}[0];

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

## Special case since Bootstrap needs to collect from recipe not immediate upstream
test_add_io_for_recipe(
    {
        file_info_href => \%file_info,
        id             => $sample_id,
        parameter_href => \%parameter,
        recipe_name    => q{dna_vcf_reformat},
        step           => q{vcf},
    }
);

test_add_io_for_recipe(
    {
        file_info_href => \%file_info,
        id             => $sample_id,
        parameter_href => \%parameter,
        recipe_name    => q{gatk_asereadcounter},
        step           => q{vcf},
    }
);
test_add_io_for_recipe(
    {
        file_info_href    => \%file_info,
        id                => $sample_id,
        parameter_href    => \%parameter,
        order_recipes_ref => [
            qw{ gatk_variantfiltration dna_vcf_reformat gatk_asereadcounter },
            $recipe_name
        ],
        recipe_name => $recipe_name,
        step        => q{vcf},
    }
);
my %sample_info;

my $is_ok = analysis_bootstrapann(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        job_id_href           => \%job_id,
        parameter_href        => \%parameter,
        profile_base_command  => $slurm_mock_cmd,
        recipe_name           => $recipe_name,
        sample_id             => $sample_id,
        sample_info_href      => \%sample_info,
    }
);

## Then return TRUE
ok( $is_ok, q{Executed analysis recipe } . $recipe_name );

## Given a dna vcf file
$active_parameter{dna_vcf_file} = catfile(qw{ a dna file.vcf });

$is_ok = analysis_bootstrapann(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        job_id_href           => \%job_id,
        parameter_href        => \%parameter,
        profile_base_command  => $slurm_mock_cmd,
        recipe_name           => $recipe_name,
        sample_id             => $sample_id,
        sample_info_href      => \%sample_info,
    }
);

## Then return TRUE
ok( $is_ok, qq{Executed analysis recipe $recipe_name with dna vcf} );

## Given a sample id where ASE has been turned off
delete $active_parameter{dna_vcf_file};
$active_parameter{no_ase_samples} = [qw{ ADM1059A1 }];
$is_ok = analysis_bootstrapann(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        job_id_href           => \%job_id,
        parameter_href        => \%parameter,
        profile_base_command  => $slurm_mock_cmd,
        recipe_name           => $recipe_name,
        sample_id             => $sample_id,
        sample_info_href      => \%sample_info,
    }
);

## Then return FALSE
is( $is_ok, 0, q{When no ASE, turn off } . $recipe_name );

done_testing();
