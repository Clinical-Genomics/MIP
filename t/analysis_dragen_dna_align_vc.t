#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ remove_tree };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
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
use MIP::Test::Fixtures qw{ test_add_io_for_recipe test_log test_mip_hashes test_standard_cli };

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
        q{MIP::Recipes::Analysis::Dragen_dna} => [qw{ analysis_dragen_dna_align_vc }],
        q{MIP::Test::Fixtures} => [qw{ test_add_io_for_recipe test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Dragen_dna qw{ analysis_dragen_dna_align_vc };

diag(   q{Test analysis_dragen_dna_align_vc from Dragen_dna.pm v}
      . $MIP::Recipes::Analysis::Dragen_dna::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $test_dir = File::Temp->newdir();
my $log      = test_log( { log_name => q{MIP}, no_screen => 1, } );

## Given analysis parameters
my $recipe_name    = q{analysis_dragen_dna_align_vc};
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
$active_parameter{dragen_hash_ref_dir_path}    = q{a_hash_dir};
$active_parameter{platform}                    = q{ILLUMINA};
$active_parameter{dragen_fastq_list_file_path} = q{an_dragen_fastq_file_path};
$active_parameter{dragen_analysis_dir}         = $test_dir;
$active_parameter{dragen_user_at_hostname}     = q{dragen@hostname};

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
        file_info_href => \%file_info,
        id             => $sample_id,
        parameter_href => \%parameter,
        recipe_name    => $recipe_name,
        step           => q{vcf},
    }
);

my %sample_info = (
    sample => {
        $sample_id => {
            file => {
                ADM1059A1_161011_HHJJCCCXY_NAATGCGC_lane7 => {
                    sequence_run_type   => q{single-end},
                    read_direction_file => {
                        ADM1059A1_161011_HHJJCCCXY_NAATGCGC_lane7_1 => {
                            flowcell       => q{HHJJCCCXY},
                            lane           => q{7},
                            sample_barcode => q{NAATGCGC},
                        },
                    },
                },
            },
        },
    },
);

my $is_ok = analysis_dragen_dna_align_vc(
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
ok( $is_ok, q{ Executed analysis recipe } . $recipe_name );

## Clean-up
remove_tree( $active_parameter{dragen_fastq_list_file_path} );

done_testing();
