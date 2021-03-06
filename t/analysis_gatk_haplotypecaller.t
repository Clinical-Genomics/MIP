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
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.04;

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
        q{MIP::Recipes::Analysis::Gatk_haplotypecaller} =>
          [qw{ analysis_gatk_haplotypecaller }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Gatk_haplotypecaller qw{ analysis_gatk_haplotypecaller };

diag(   q{Test analysis_gatk_haplotypecaller from Gatk_haplotypecaller.pm v}
      . $MIP::Recipes::Analysis::Gatk_haplotypecaller::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

## Given analysis parameters
my $recipe_name    = q{gatk_haplotypecaller};
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
$active_parameter{gatk_haplotypecaller_emit_ref_confidence}    = q{GVCF};
$active_parameter{gatk_haplotypecaller_linked_de_bruijn_graph} = 0;

my %file_info = test_mip_hashes(
    {
        mip_hash_name => q{file_info},
        recipe_name   => $recipe_name,
    }
);

SAMPLE_ID:
foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

    %{ $file_info{io}{TEST}{$sample_id}{$recipe_name} } = test_mip_hashes(
        {
            mip_hash_name => q{io},
        }
    );
}
## Set correct outfile path
CONTIG:
foreach my $contig ( @{ $file_info{contigs} } ) {

    $file_info{io}{TEST}{$sample_id}{$recipe_name}{in}{file_path_href}{$contig} =
      q{a_file.bam};
}

my %job_id;
my %parameter = test_mip_hashes(
    {
        mip_hash_name => q{recipe_parameter},
        recipe_name   => $recipe_name,
    }
);
@{ $parameter{cache}{order_recipes_ref} } = ($recipe_name);
$parameter{$recipe_name}{outfile_suffix} = q{.vcf};

my %sample_info = test_mip_hashes(
    {
        mip_hash_name => q{qc_sample_info},
        recipe_name   => $recipe_name,
    }
);

my $is_ok = analysis_gatk_haplotypecaller(
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
