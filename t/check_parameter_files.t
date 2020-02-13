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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Active_parameter} => [qw{ check_parameter_files }],
        q{MIP::Io::Read}         => [qw{ read_from_file }],
        q{MIP::Parameter}        => [qw{ get_parameter_attribute }],
        q{MIP::Test::Fixtures}   => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ check_parameter_files };
use MIP::Io::Read qw{ read_from_file };
use MIP::Parameter qw{ get_parameter_attribute };

diag(   q{Test check_parameter_files from Active_parameter.pm v}
      . $MIP::Active_parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

## Given files to check for existence when stored as scalar, array and hash
my @order_parameters = qw{ gatk_baserecalibration_known_sites
  gatk_genotypegvcfs_ref_gvcf
  gatk_variantrecalibration_resource_indel
  human_genome_reference
  sv_vcfparser_select_file
  vcfparser_select_file };

my %active_parameter = (

    # To test array parameter
    gatk_baserecalibration_known_sites =>
      [ catfile( $Bin, qw{ data references grch37_dbsnp_-138-.vcf} ), ],
    gatk_genotypegvcfs_ref_gvcf => q{test_file},

    # To test scalar parameter
    human_genome_reference =>
      catfile( $Bin, qw{data references grch37_homo_sapiens_-d5-.fasta} ),
    mip                       => 1,
    gatk_baserecalibration    => 1,
    gatk_genotypegvcfs        => 1,
    gatk_variantrecalibration => 1,
    sv_vcfparser              => 0,

    # To test hash parameter
    gatk_variantrecalibration_resource_indel => {
        catfile( $Bin, qw{data references grch37_dbsnp_-138-.vcf} ) =>
          q{dbsnp,known=true,training=false,truth=false,prior=2.0},
    },
    sv_vcfparser_select_file => q{test_file},
);
my $consensus_analysis_type = q{wgs};
my %parameter               = read_from_file(
    {
        format => q{yaml},
        path   => catfile( dirname($Bin), qw{ definitions rd_dna_parameters.yaml} ),
    }
);

PARAMETER:
foreach my $parameter_name (@order_parameters) {

    my %attribute = get_parameter_attribute(
        {
            parameter_href => \%parameter,
            parameter_name => $parameter_name,
        }
    );

    check_parameter_files(
        {
            active_parameter_href   => \%active_parameter,
            associated_recipes_ref  => $attribute{associated_recipe},
            build_status            => $attribute{build_file},
            consensus_analysis_type => $consensus_analysis_type,
            parameter_exists_check  => $attribute{exists_check},
            parameter_name          => $parameter_name,
        }
    );
}

is( $active_parameter{gatk_genotypegvcfs_ref_gvcf},
    q{test_file}, q{Returned for not required exome mode parameter} );

is( $active_parameter{sv_vcfparser_select_file},
    q{test_file}, q{Returned for not associated recipe} );

is( $active_parameter{vcfparser_select_file},
    undef, q{Returned for not active parameter} );

done_testing();
