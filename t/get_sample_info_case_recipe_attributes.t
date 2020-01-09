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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Sample_info}    => [qw{ get_sample_info_case_recipe_attributes }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ get_sample_info_case_recipe_attributes };

diag(   q{Test get_sample_info_case_recipe_attributes from Sample_info.pm v}
      . $MIP::Sample_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a recipe when called with an attribute
my $recipe_name = q{bwa_mem};
my %sample_info = (
    recipe => {
        $recipe_name => {
            outdirectory => q{an_outdir},
            outfile      => q{an_outfile},
            path         => catfile(qw{ an_outdir an_outfile}),
            version      => q{major.minor.patch},
        },
    },
);

my $got_attribute = get_sample_info_case_recipe_attributes(
    {
        attribute        => q{path},
        recipe_name      => $recipe_name,
        sample_info_href => \%sample_info,
    }
);

## Then return attribute for recipe
is( $got_attribute, catfile(qw{ an_outdir an_outfile}), q{Got recipe attribute} );

## Given a recipe when called withiut an attribute
my %got_attribute = get_sample_info_case_recipe_attributes(
    {
        recipe_name      => $recipe_name,
        sample_info_href => \%sample_info,
    }
);

## Then return attribute hash for recipes
is_deeply(
    \%got_attribute,
    \%{ $sample_info{recipe}{$recipe_name} },
    q{Got recipe attribute hash}
);

## Given a recipe when sample_info attribute is undef
delete $sample_info{recipe}{$recipe_name}{path};

my $is_undef = get_sample_info_case_recipe_attributes(
    {
        attribute        => q{path},
        recipe_name      => $recipe_name,
        sample_info_href => \%sample_info,
    }
);

## Then return attribute for recipe
is( $is_undef, undef, q{Return undef for missing attribute} );

done_testing();
