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
use Test::Trap qw{ :stderr :output(systemsafe) };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_add_io_for_recipe test_log test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Analysis::Sv_annotate} => [qw{ analysis_sv_annotate }],
        q{MIP::Test::Fixtures} => [qw{ test_add_io_for_recipe test_log test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Sv_annotate qw{ analysis_sv_annotate };

diag(   q{Test analysis_sv_annotate from Sv_annotate.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $FREQ_CUTOFF => 0.40;

my $log = test_log( { log_name => q{MIP} } );

## Given analysis parameters
my $recipe_name    = q{sv_annotate};
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
$active_parameter{sv_frequency_filter}           = 1;
$active_parameter{fqf_bcftools_filter_threshold} = $FREQ_CUTOFF;
$active_parameter{sv_vcfanno_config} =
  catfile( $Bin, qw{ data references grch37_frequency_vcfanno_filter_config_-v1.0-.toml } );
$active_parameter{sv_svdb_query}          = 1;
$active_parameter{sv_svdb_query_db_files} = { a_file => q{a_file|AF|AC|in_AF|in_AC|1}, };
@{ $active_parameter{sv_fqa_vcfanno_filters} } = (qw{ GNOMADAF_popmax GNOMADAF });

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
        id             => $case_id,
        parameter_href => \%parameter,
        recipe_name    => $recipe_name,
        step           => q{vcf},
    }
);

my %sample_info;

my @return = trap {
    analysis_sv_annotate(
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
    )
};
## Then return TRUE
ok( $return[0], q{Executed analysis recipe} . $recipe_name );

## Given request to filter with missing vcfanno annotation
push @{ $active_parameter{sv_fqa_vcfanno_filters} }, q{missingtag};

trap {
    analysis_sv_annotate(
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
    )
};

## Then print warning and which annotation that is missing to log
is( $trap->leaveby, q{return}, q{Don't fail missing annotation} );
like( $trap->stderr, qr/missingtag/xms, q{Print missing annotation tag} );

done_testing();
