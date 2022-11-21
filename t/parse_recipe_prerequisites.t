#! /usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
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
use MIP::Test::Fixtures qw{ test_log test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipe}         => [qw{ parse_recipe_prerequisites }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipe qw{ parse_recipe_prerequisites };

diag(   q{Test parse_recipe_prerequisites from Recipe.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( { no_screen => 1, } );

## Given a recipe name and active parameter hash
my $recipe_name = q{deepvariant};
my %active_parameter =
  test_mip_hashes( { mip_hash_name => q{active_parameter}, recipe_name => $recipe_name, } );
$active_parameter{conda_path} = catdir(qw{path to conda});

## Given a parameter hash
my %parameter =
  test_mip_hashes( { mip_hash_name => q{define_parameter}, recipe_name => $recipe_name, } );

## Given a recipe with a chain, file_tag and outfile_suffix
$parameter{$recipe_name}{chain}          = q{TEST};
$parameter{$recipe_name}{file_tag}       = q{_deepvar};
$parameter{$recipe_name}{outfile_suffix} = q{.g.vcf.gz};

## When parsing recipe prerequisites
my %recipe = parse_recipe_prerequisites(
    {
        active_parameter_href => \%active_parameter,
        parameter_href        => \%parameter,
        recipe_name           => $recipe_name,
    }
);

## Then return recipe prerequisites hash
my %expected_recipe = (
    core_number    => 35,
    file_tag       => q{_deepvar},
    gpu_number     => 1,
    job_id_chain   => q{TEST},
    load_env_ref   => [qw{ source path/to/conda/etc/profile.d/conda.sh ; conda activate test }],
    memory         => 175,
    mode           => 2,
    outfile_suffix => q{.g.vcf.gz},
    time           => 10,
);

is_deeply( \%recipe, \%expected_recipe, q{Got recipe prerequisites hash} );

done_testing();
