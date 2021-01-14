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
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Sample_info}    => [qw{ get_sample_info_sample_recipe_attributes }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ get_sample_info_sample_recipe_attributes };

diag(   q{Test get_sample_info_sample_recipe_attributes from Sample_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a sample_id, recipe when supplied with an infile and attribute
my $attribute   = q{regexp};
my $infile      = q{an_infile};
my $recipe_name = q{qccollect};
my $sample_id   = q{sample_1};
my %sample_info = (
    sample => {
        $sample_id => {
            recipe => {
                $recipe_name        => { $infile => { $attribute => q{regexp_file}, }, },
                q{recipe_no_infile} => { path    => q{a_path}, },
            },
        },
    },
);

my $got_attribute_value = get_sample_info_sample_recipe_attributes(
    {
        attribute        => $attribute,
        infile           => $infile,
        recipe_name      => $recipe_name,
        sample_id        => $sample_id,
        sample_info_href => \%sample_info,
    }
);

## Then return attribute for recipe
is( $got_attribute_value, q{regexp_file}, q{Got recipe attribute value} );

## Given a sample_id, recipe when supplied with an infile but no attribute
my %got_attribute = get_sample_info_sample_recipe_attributes(
    {
        infile           => $infile,
        recipe_name      => $recipe_name,
        sample_id        => $sample_id,
        sample_info_href => \%sample_info,
    }
);

## Then return attribute hash for recipes for infile
is_deeply(
    \%got_attribute,
    \%{ $sample_info{sample}{$sample_id}{recipe}{$recipe_name}{$infile} },
    q{Got recipe attribute hash for infile}
);

# Given a recipe when sample_info attribute is undef
delete $sample_info{sample}{$sample_id}{recipe}{$recipe_name}{$infile}{$attribute};

my $is_undef = get_sample_info_sample_recipe_attributes(
    {
        attribute        => $attribute,
        infile           => $infile,
        recipe_name      => $recipe_name,
        sample_id        => $sample_id,
        sample_info_href => \%sample_info,
    }
);

## Then return attribute for recipe
is( $is_undef, undef, q{Return undef for missing attribute} );

## Given a recipe with no infile

%got_attribute = get_sample_info_sample_recipe_attributes(
    {
        infile           => q{path},
        recipe_name      => q{recipe_no_infile},
        sample_id        => $sample_id,
        sample_info_href => \%sample_info,
    }
);

## Then return attribute hash for recipes
is_deeply(
    \%got_attribute,
    \%{ $sample_info{sample}{$sample_id}{recipe}{recipe_no_infile} },
    q{Got recipe attribute hash for recipe}
);

done_testing();
