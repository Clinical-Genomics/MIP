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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Analysis}       => [qw{ set_ase_chain_recipes }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Analysis qw{ set_ase_chain_recipes };

diag(   q{Test set_ase_chain_recipes from Analysis.pm v}
      . $MIP::Analysis::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 1, } );

## Given allelic specific expression recipes
my %active_parameter = (
    dna_vcf_reformat       => 1,
    gatk_baserecalibration => 1,
    gatk_haplotypecaller   => 1,
    gatk_splitncigarreads  => 1,
    gatk_variantfiltration => 1,
);

## When no DNA VCF file is supplied
set_ase_chain_recipes( { active_parameter_href => \%active_parameter, } );

## Then turn off DNA VCF reformat recipe
is( $active_parameter{dna_vcf_reformat}, 0, q{Turn off DNA VCF reformat recipe} );

## Reset for next test
$active_parameter{dna_vcf_reformat} = 1;

## When DNA VCF file is supplied
$active_parameter{dna_vcf_file} = catfile(qw{ my dna_variants.vcf });

set_ase_chain_recipes( { active_parameter_href => \%active_parameter, } );

my %expected_recipe_mode = (
    dna_vcf_file           => catfile(qw{ my dna_variants.vcf }),
    dna_vcf_reformat       => 1,
    gatk_baserecalibration => 0,
    gatk_haplotypecaller   => 0,
    gatk_splitncigarreads  => 0,
    gatk_variantfiltration => 0,
);

## Then keep DNA VCF reformat recipe as is and turn off ASE recipes
is_deeply( \%active_parameter, \%expected_recipe_mode, q{Switched modes of ASE recipes} );

done_testing();
