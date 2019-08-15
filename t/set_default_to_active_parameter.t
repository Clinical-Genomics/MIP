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
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.6;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA   => q{,};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Set::Parameter} => [qw{ set_default_to_active_parameter }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::Parameter qw{ set_default_to_active_parameter };
use MIP::File::Format::Yaml qw{ load_yaml };

diag(   q{Test set_default_to_active_parameter from Set::Parameter.pm v}
      . $MIP::Set::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create log object
my $log = test_log( {} );

my @order_parameters =
  qw{ bcftools_mpileup_filter_variant bwa_mem bwa_mem_bamstats gatk_genotypegvcfs_ref_gvcf gatk_variantrecalibration_resource_indel markduplicates sv_vcfparser_range_feature_file };

my %active_parameter = (
    bwa_mem_bamstats                          => 0,
    expansionhunter                           => 1,
    gatk_genotypegvcfs_ref_gvcf               => q{test_file},
    markduplicates_picardtools_markduplicates => 1,
    mip                                       => 1,
    bwa_mem                                   => 0,
    gatk_baserecalibration                    => 1,
    gatk_genotypegvcfs                        => 1,
    gatk_variantrecalibration                 => 1,
    sv_vcfparser                              => 1,
);

my %parameter = load_yaml(
    {
        yaml_file => catfile( dirname($Bin), qw{ definitions rd_dna_parameters.yaml} ),
    }
);

$parameter{cache}{consensus_analysis_type} = q{wgs};

PARAMETER:
foreach my $parameter_name (@order_parameters) {

    set_default_to_active_parameter(
        {
            active_parameter_href => \%active_parameter,
            associated_recipes_ref =>
              \@{ $parameter{$parameter_name}{associated_recipe} },
            log            => $log,
            parameter_href => \%parameter,
            parameter_name => $parameter_name,
        }
    );
}

is( $active_parameter{gatk_genotypegvcfs_ref_gvcf},
    q{test_file}, q{Returned for not required exome mode parameter} );

is( $active_parameter{bwa_mem_bamstats},
    q{1}, q{Set default for non active associated_recipe parameter} );

is( $active_parameter{markduplicates_picardtools_markduplicates},
    q{1}, q{Did not set default for not defined associated_recipe parameter} );

is( $active_parameter{bcftools_mpileup_filter_variant},
    0, q{Set default for scalar parameter} );

my %expected_resource_indel = (
    q{grch37_dbsnp_-138-.vcf} => q{dbsnp,known=true,training=false,truth=false,prior=2.0},
    q{grch37_mills_and_1000g_-gold_standard_indels-.vcf} =>
      q{mills,known=false,training=true,truth=true,prior=12.0},
);
is_deeply( \%{ $active_parameter{gatk_variantrecalibration_resource_indel} },
    \%expected_resource_indel, 'Set default for hash parameter' );

is( $active_parameter{sv_vcfparser_range_feature_file},
    undef, q{Skipped no default and not mandatory parameter} );

## Given a mandatory and active recipe
trap {
    set_default_to_active_parameter(
        {
            active_parameter_href => \%active_parameter,
            associated_recipes_ref =>
              \@{ $parameter{sv_vcfparser_select_file}{associated_recipe} },
            log            => $log,
            parameter_href => \%parameter,
            parameter_name => q{sv_vcfparser_select_file},
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if mandatory parameter for active recipe is not set} );
like(
    $trap->stderr,
    qr/Supply \s+ '-/xms,
    q{Throw fatal log message if mandatory parameter for active recipe is not set}
);

done_testing();
