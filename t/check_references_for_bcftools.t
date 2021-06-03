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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Check::Reference} => [qw{ check_references_for_bcftools }],
        q{MIP::Test::Fixtures}   => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Reference qw{ check_references_for_bcftools };

diag(   q{Test check_references_for_bcftools from Reference.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Create log object
my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

## Given reference vcf files associated with different recipes, check if the files need bcftools processing
my %active_parameter_test = (
    binary_path                        => { bcftools => q{bcftools}, },
    gatk_baserecalibration             => 1,
    gatk_baserecalibration_known_sites => [
        catfile( $Bin, qw{ data references grch37_dbsnp_-138-.vcf } ),
        catfile( $Bin, qw{ data references grch37_1000g_indels_-phase1-.vcf } ),
        catfile( $Bin, qw{ data references grch37_mills_and_1000g_indels_-gold_standard-.vcf } )
    ],
    gatk_variantevalall    => 1,
    gatk_varianteval_dbsnp => catfile( $Bin, qw{ data references grch37_dbsnp_-138_esa_129-.vcf } ),
    gatk_variantevalexome  => 1,
    gatk_variantrecalibration                => 1,
    gatk_variantrecalibration_resource_indel => {
        q{grch37_dbsnp_-138-.vcf} => q{dbsnp,known=true,training=false,truth=false,prior=2.0},
        q{grch37_mills_and_1000g_-gold_standard_indels-.vcf} =>
          q{mills,known=false,training=true,truth=true,prior=12.0},
    },
    variant_annotation => 1,
    vcfanno_config     =>
      catfile( $Bin, qw{ data references grch37_frequency_vcfanno_filter_config_-v1.0-.toml } ),
);

my %parameter_test = (
    gatk_baserecalibration_known_sites => {
        associated_recipe => [qw{ gatk_baserecalibration }],
        data_type         => q{ARRAY},
    },
    gatk_varianteval_dbsnp => {
        associated_recipe => [qw{ gatk_variantevalall gatk_variantevalexome }],
        data_type         => q{SCALAR},
    },
    gatk_variantrecalibration_resource_indel =>
      { associated_recipe => [qw{ gatk_variantrecalibration }], data_type => q{HASH}, },
    vcfanno_config => {
        associated_recipe => [qw{ variant_annotation }],
        data_type         => q{SCALAR},
    },
);

my @bcftools_references_test =
  qw{ gatk_baserecalibration_known_sites gatk_varianteval_dbsnp gatk_varianteval_dbsnp gatk_variantrecalibration_resource_indel vcfanno_config };

my @refs_to_process = check_references_for_bcftools(
    {
        active_parameter_href   => \%active_parameter_test,
        parameter_href          => \%parameter_test,
        bcftools_references_ref => \@bcftools_references_test,
    }
);

isnt( @refs_to_process, 0, q{Test passed, } . @refs_to_process . q{ references to process} );

done_testing();
