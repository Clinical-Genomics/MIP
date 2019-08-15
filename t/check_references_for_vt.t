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
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Check::Reference} => [qw{ check_references_for_vt }],
        q{MIP::Test::Fixtures}   => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Reference qw{ check_references_for_vt };

diag(   q{Test check_references_for_vt from Reference.pm v}
      . $MIP::Check::Reference::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Create log object
my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

my %active_parameter_test = (
    fqf_vcfanno_config => catfile(
        $Bin, qw{ data references grch37_frequency_vcfanno_filter_config_-v1.0-.toml }
    ),
    frequency_filter                   => 1,
    gatk_baserecalibration             => 1,
    gatk_baserecalibration_known_sites => [
        catfile( $Bin, qw{ data references grch37_dbsnp_-138-.vcf } ),
        catfile( $Bin, qw{ data references grch37_1000g_indels_-phase1-.vcf } ),
        catfile(
            $Bin, qw{ data references grch37_mills_and_1000g_indels_-gold_standard-.vcf }
        )
    ],
    gatk_variantevalall => 1,
    gatk_varianteval_dbsnp =>
      catfile( $Bin, qw{ data references grch37_dbsnp_-138_esa_129-.vcf } ),
    gatk_variantevalexome    => 1,
    snpeff                   => 1,
    snpsift_annotation_files => {
        catfile(
            $Bin, qw{ data references grch37_anon-swegen_snp_-1000samples-.vcf.gz }
        ) => q{SWEREF},
        catfile( $Bin, qw{ data references grch37_exac_reheader_-r0.3.1-.vcf.gz } ) =>
          q{EXAC},
        catfile( $Bin, qw{ data references grch37_gnomad.genomes_-r2.0.1-.vcf.gz } ) =>
          q{GNOMAD},
    },
);

my %parameter_test = (
    fqf_vcfanno_config => {
        associated_recipe => [qw{ frequency_filter }],
        data_type         => q{SCALAR},
    },
    gatk_baserecalibration_known_sites => {
        associated_recipe => [qw{ gatk_baserecalibration }],
        data_type         => q{ARRAY},
    },
    gatk_varianteval_dbsnp => {
        associated_recipe => [qw{ gatk_variantevalall gatk_variantevalexome }],
        data_type         => q{SCALAR},
    },
    snpsift_annotation_files =>
      { associated_recipe => [qw{ snpeff }], data_type => q{HASH}, },
);

my @vt_references_test =
  qw{ fqf_vcfanno_config gatk_baserecalibration_known_sites gatk_varianteval_dbsnp gatk_varianteval_dbsnp snpsift_annotation_files};

my @refs_to_process = check_references_for_vt(
    {
        active_parameter_href => \%active_parameter_test,
        log                   => $log,
        parameter_href        => \%parameter_test,
        vt_references_ref     => \@vt_references_test,
    }
);

isnt( @refs_to_process, 0,
    q{Test passed, } . @refs_to_process . q{ references to process} );

done_testing();
