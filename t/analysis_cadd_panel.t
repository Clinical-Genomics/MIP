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
use MIP::Test::Fixtures qw{ test_add_io_for_recipe test_log test_mip_hashes };

## Constants
Readonly my $GENOME_BUILD_VERSION_37 => 37;
Readonly my $GENOME_BUILD_VERSION_38 => 38;

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Analysis::Cadd} => [qw{ analysis_cadd_panel }],
        q{MIP::Test::Fixtures} => [qw{ test_add_io_for_recipe test_log test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Cadd qw{ analysis_cadd_panel };

diag(   q{Test analysis_cadd_panel from Cadd.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

## Given build parameters
my $recipe_name    = q{cadd};
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
$active_parameter{cadd_column_names} = [qw{ col_1 col_2 }];

my %file_info = test_mip_hashes(
    {
        mip_hash_name => q{file_info},
        recipe_name   => $recipe_name,
    }
);

$file_info{human_genome_reference_source} = q{grch};

my %job_id;
my %parameter = test_mip_hashes(
    {
        mip_hash_name => q{recipe_parameter},
        recipe_name   => $recipe_name,
    }
);

test_add_io_for_recipe(
    {
        file_info_href => \%file_info,
        id             => $case_id,
        parameter_href => \%parameter,
        recipe_name    => $recipe_name,
        step           => q{vcf},
    }
);

my %sample_info;

GENOME_VERSION:
foreach my $genome_version ( $GENOME_BUILD_VERSION_37, $GENOME_BUILD_VERSION_38 ) {

    $file_info{human_genome_reference_version} = $genome_version;

    my $is_ok = analysis_cadd_panel(
        {
            active_parameter_href   => \%active_parameter,
            case_id                 => $case_id,
            file_info_href          => \%file_info,
            job_id_href             => \%job_id,
            parameter_href          => \%parameter,
            profile_base_command    => $slurm_mock_cmd,
            recipe_name             => $recipe_name,
            sample_info_href        => \%sample_info,
        }
    );
    ## Then return TRUE
    ok( $is_ok, qq{Executed analysis recipe $recipe_name for grch$genome_version} );
}
done_testing();
