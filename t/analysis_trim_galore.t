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
use MIP::Test::Fixtures qw{ test_add_io_for_recipe test_log test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

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
        q{MIP::Recipes::Analysis::Trim_galore} => [qw{ analysis_trim_galore }],
        q{MIP::Test::Fixtures} => [qw{ test_add_io_for_recipe test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Trim_galore qw{ analysis_trim_galore };

diag(   q{Test analysis_trim_galore from Trim_galore.pm v}
      . $MIP::Recipes::Analysis::Trim_galore::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

## Given analysis parameters
my $recipe_name    = q{trim_galore_ar};
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

test_add_io_for_recipe(
    {
        file_info_href    => \%file_info,
        id                => $sample_id,
        parameter_href    => \%parameter,
        recipe_name       => $recipe_name,
        step              => q{fastq},
    }
);

my %sample_info = (
    sample => {
        $sample_id => {
            file => {
                ADM1059A1_161011_TestFilev2_GAGATTCC_lane1 => {
                    sequence_run_type   => q{single-end},
                    read_direction_file => {
                        ADM1059A1_161011_TestFilev2_GAGATTCC_lane1_1 => {
                            flowcell       => q{TestFilev2},
                            lane           => q{1},
                            sample_barcode => q{GAGATTC},
                        },
                    },
                },
                ADM1059A1_161012_TestFilev3_GAGATTGG_lane1 => {
                    sequence_run_type   => q{paired-end},
                    read_direction_file => {
                        ADM1059A1_161012_TestFilev3_GAGATTGG_lane1_1 => {
                            flowcell       => q{TestFilev3},
                            lane           => q{1},
                            sample_barcode => q{GAGATTG},
                        },
                        ADM1059A1_161012_TestFilev3_GAGATTGG_lane1_2 => {
                            flowcell       => q{TestFilev3},
                            lane           => q{1},
                            sample_barcode => q{GAGATTG},
                        },
                    },
                },
            },
        },
    },
);

my $is_ok = analysis_trim_galore(
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

done_testing();
