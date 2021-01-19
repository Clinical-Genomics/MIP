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
        q{MIP::Recipes::Analysis::Vcf2cytosure} => [qw{ analysis_vcf2cytosure }],
        q{MIP::Test::Fixtures} => [qw{ test_add_io_for_recipe test_log test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Vcf2cytosure qw{ analysis_vcf2cytosure };

diag(   q{Test analysis_vcf2cytosure from Vcf2cytosure.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

## Given analysis parameters
my $recipe_name    = q{vcf2cytosure};
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
my $case_id = $active_parameter{case_id};
$active_parameter{vcf2cytosure_maxbnd}   = 1;
$active_parameter{vcf2cytosure_var_size} = 1;

my %file_info = test_mip_hashes(
    {
        mip_hash_name => q{file_info},
        recipe_name   => $recipe_name,
    }
);
my %parameter = test_mip_hashes(
    {
        mip_hash_name => q{recipe_parameter},
        recipe_name   => $recipe_name,
    }
);

## Special case since vcf2cytosure needs to collect from recipes not immediate upstream
my @order_recipes = ( qw{ markduplicates sv_annotate }, $recipe_name );

SAMPLE_ID:
foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

    test_add_io_for_recipe(
        {
            file_info_href => \%file_info,
            id             => $sample_id,
            parameter_href => \%parameter,
            recipe_name    => q{markduplicates},
            step           => q{bam},
        }
    );
}

test_add_io_for_recipe(
    {
        file_info_href => \%file_info,
        id             => $case_id,
        parameter_href => \%parameter,
        recipe_name    => q{sv_annotate},
        step           => q{vcf},
    }
);

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

my %job_id;

my %sample_info = (
    sample => {
        ADM1059A1 => {
            sample_display_name => q{boris},
            sex                 => q{male},
        },
        ADM1059A2 => {
            sample_display_name => q{stefan},
            sex                 => q{unknown},
        },
        ADM1059A3 => {
            sample_display_name => q{venus},
            sex                 => q{female},
        },
    }
);

my $is_ok = analysis_vcf2cytosure(
    {
        active_parameter_href => \%active_parameter,
        bin_size              => 1,
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

done_testing();
