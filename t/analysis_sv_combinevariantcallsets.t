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
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

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
        q{MIP::Recipes::Analysis::Sv_combinevariantcallsets} =>
          [qw{ analysis_sv_combinevariantcallsets }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Sv_combinevariantcallsets
  qw{ analysis_sv_combinevariantcallsets };

diag(   q{Test analysis_sv_combinevariantcallsets from Sv_combinevariantcallsets.pm v}
      . $MIP::Recipes::Analysis::Sv_combinevariantcallsets::VERSION
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

%{ $file_info{io}{TEST}{$case_id}{$recipe_name} } = test_mip_hashes(
    {
        mip_hash_name => q{io},
    }
);

SV_CALLER:
foreach my $sv_caller (@structural_variant_callers) {
    %{ $file_info{io}{TEST}{$case_id}{$sv_caller} } = test_mip_hashes(
        {
            mip_hash_name => q{io},
        }
    );
}

my %infile_lane_prefix;
my %job_id;
my %parameter = test_mip_hashes(
    {
        mip_hash_name => q{recipe_parameter},
        recipe_name   => $recipe_name,
    }
);
@{ $parameter{cache}{order_recipes_ref} } =
  ( qw{ tiddit manta delly_reformat cnvnator_ar }, $recipe_name );
$parameter{cache}{structural_variant_callers} =
  [qw{ tiddit manta delly_reformat cnvnator_ar }];

SV_CALLER:
foreach my $sv_caller (@structural_variant_callers) {
    $parameter{$sv_caller}{chain} = uc q{test};
}
$parameter{$recipe_name}{outfile_suffix} = q{.vcf};

SV_CALLER:
foreach my $sv_caller (@structural_variant_callers) {

    $parameter{$sv_caller}{outfile_suffix} = q{.vcf};
}

my %sample_info;

my $is_ok = analysis_sv_combinevariantcallsets(
    {
        active_parameter_href   => \%active_parameter,
        case_id                 => $case_id,
        file_info_href          => \%file_info,
        infile_lane_prefix_href => \%infile_lane_prefix,
        job_id_href             => \%job_id,
        parameter_href          => \%parameter,
        profile_base_command    => $slurm_mock_cmd,
        recipe_name             => $recipe_name,
        sample_info_href        => \%sample_info,
    }
);

## Then return TRUE
ok( $is_ok, q{ Executed analysis recipe } . $recipe_name );

## Given a single samples
@{ $active_parameter{sample_ids} } = $active_parameter{sample_ids}[0];

$is_ok = analysis_sv_combinevariantcallsets(
    {
        active_parameter_href   => \%active_parameter,
        case_id                 => $case_id,
        file_info_href          => \%file_info,
        infile_lane_prefix_href => \%infile_lane_prefix,
        job_id_href             => \%job_id,
        parameter_href          => \%parameter,
        profile_base_command    => $slurm_mock_cmd,
        recipe_name             => $recipe_name,
        sample_info_href        => \%sample_info,
    }
);

## Then return TRUE
ok( $is_ok, q{ Executed analysis recipe } . $recipe_name . q{ single sample} );

done_testing();
