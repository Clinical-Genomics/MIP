#!/usr/bin/env perl

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
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Update::Recipes} => [qw{ update_recipe_mode_with_start_with }],
        q{MIP::Test::Fixtures}  => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Update::Recipes qw{ update_recipe_mode_with_start_with };

diag(   q{Test update_recipe_mode_with_start_with from Recipes.pm v}
      . $MIP::Update::Recipes::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %active_parameter = (
    fastqc_ar  => 1,
    bwa_mem    => 2,
    star_aln   => 0,
    multiqc_ar => 2,
);
my @recipes            = qw{ fastqc_ar bwa_mem star_aln multiqc_ar };
my @start_with_recipes = qw{ bwa_mem multiqc_ar };

update_recipe_mode_with_start_with(
    {
        active_parameter_href  => \%active_parameter,
        recipes_ref            => \@recipes,
        start_with_recipes_ref => \@start_with_recipes,
    }
);

is( $active_parameter{fastqc_ar}, 2, q{Udated upstreams dependencies recipe mode} );

is( $active_parameter{bwa_mem}, 1, q{Udated start with recipe mode} );

is( $active_parameter{multiqc_ar}, 1, q{Udated downstream dependencies recipe mode} );

is( $active_parameter{star_aln}, 0, q{Did not update switched off recipe } );

done_testing();
