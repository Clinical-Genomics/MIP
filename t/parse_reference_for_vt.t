#! /usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ remove_tree };
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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures
  qw{ test_add_io_for_recipe test_constants test_log test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Parse::Reference} => [qw{ parse_reference_for_vt }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parse::Reference qw{ parse_reference_for_vt };

diag(   q{Test parse_reference_for_vt from Reference.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( {} );

## Given analysis parameters
my $recipe_name      = q{vt_core};
my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
        recipe_name   => $recipe_name,
    }
);
$active_parameter{$recipe_name}                     = 2;
$active_parameter{recipe_core_number}{$recipe_name} = 1;
$active_parameter{recipe_time}{$recipe_name}        = 1;
$active_parameter{outscript_dir} = catdir( $Bin, q{test_vt-core} );
$active_parameter{outdata_dir}   = catdir( $Bin, q{test_vt-core} );
$active_parameter{vt_decompose}  = 1;
$active_parameter{vt_normalize}  = 1;
$active_parameter{decompose_normalize_references} = [q{gatk_varianteval_dbsnp}];
$active_parameter{gatk_variantevalall}            = 1;
$active_parameter{gatk_varianteval_dbsnp} =
  catfile( $Bin, qw{ data references grch37_mills_and_1000g_-gold_standard_indels-.vcf} );

my $case_id = $active_parameter{case_id};

my %file_info = test_mip_hashes(
    {
        mip_hash_name => q{file_info},
        recipe_name   => $recipe_name,
    }
);

my %parameter = test_mip_hashes(
    {
        mip_hash_name => q{recipe_parameter},
        recipe_name   => $recipe_name,
    }
);
$parameter{gatk_varianteval_dbsnp}{associated_recipe} = [q{gatk_variantevalall}];
$parameter{gatk_varianteval_dbsnp}{data_type}         = q{SCALAR};
my %process_return = (
    buffers_ref        => [],
    error_messages_ref => [],
    stderrs_ref        => [],
    stdouts_ref        => [],
    success            => 0,
);
test_constants( { test_process_return_href => \%process_return } );

test_add_io_for_recipe(
    {
        file_info_href => \%file_info,
        id             => $case_id,
        parameter_href => \%parameter,
        recipe_name    => $recipe_name,
        step           => q{vcf},
    }
);

## Then run sub and process references
trap {
    parse_reference_for_vt(
        {
            active_parameter_href => \%active_parameter,
            job_id_href           => {},
            parameter_href        => \%parameter,
        }
    );
};

is( $trap->leaveby, q{return}, q{Run sub} );
like( $trap->stderr, qr/Cannot \s+ detect \s+ /xms, q{Log Vt fail} );
remove_tree( catdir( $Bin, q{test_vt-core} ) );

done_testing();
