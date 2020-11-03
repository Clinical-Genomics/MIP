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
        q{MIP::Recipes::Analysis::Glnexus} => [qw{ analysis_glnexus }],
        q{MIP::Test::Fixtures} =>
          [qw{ test_add_io_for_recipe test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Glnexus qw{ analysis_glnexus };

diag(   q{Test analysis_glnexus from Glnexus.pm v}
      . $MIP::Recipes::Analysis::Glnexus::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

## Given analysis parameters
my $recipe_name    = q{Glnexus};
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

SAMPLE_ID:
foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

    test_add_io_for_recipe(
        {
            file_info_href => \%file_info,
            id             => $sample_id,
            parameter_href => \%parameter,
            recipe_name    => $recipe_name,
            step           => q{vcf},
        }
    );
}

my %sample_info = test_mip_hashes(
    {
        mip_hash_name => q{qc_sample_info},
        recipe_name   => $recipe_name,
    }
);

CONSENSUS_ANALYSIS:
foreach my $analysis_type (qw { mixed panel wgs wes }) {
    $parameter{cache}{consensus_analysis_type} = $analysis_type;
    my $is_ok = analysis_glnexus(
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
    ok( $is_ok,
        q{ Executed analysis recipe } . $recipe_name . q{ with type } . $analysis_type );
}

# Test with one sample
my %single_sample_info = test_mip_hashes(
    {
        mip_hash_name => q{qc_single_sample_info},
        recipe_name   => $recipe_name,
    }
);

my %single_sample_active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter_single_sample},
        recipe_name   => $recipe_name,
    }
);

my %single_sample_file_info = test_mip_hashes(
    {
        mip_hash_name => q{file_info_single_sample},
        recipe_name   => $recipe_name,
    }
);

$single_sample_file_info{io}{TEST}{"643594-miptest"}{Glnexus}{out}{file_paths} = $file_info{io}{TEST}{"643594-miptest"}{Glnexus}{out}{file_paths};
$parameter{cache}{consensus_analysis_type} = q{wgs};
my $is_ok = analysis_glnexus(
    {
        active_parameter_href => \%single_sample_active_parameter,
        case_id               => $case_id,
        file_info_href        => \%single_sample_file_info,
        job_id_href           => \%job_id,
        parameter_href        => \%parameter,
        profile_base_command  => $slurm_mock_cmd,
        recipe_name           => $recipe_name,
        sample_info_href      => \%single_sample_info,
    }
);
## Then return TRUE
ok( $is_ok, q{ Executed analysis recipe } . $recipe_name . q{ with a single sample } );

done_testing();
