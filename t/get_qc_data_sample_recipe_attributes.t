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
        q{MIP::Qc_data}        => [qw{ get_qc_data_sample_recipe_attributes }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qc_data qw{ get_qc_data_sample_recipe_attributes };

diag(   q{Test get_qc_data_sample_recipe_attributes from Qc_data.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a sample_id, recipe when supplied with an infile and attribute
my $attribute   = q{gender};
my $infile      = q{an_infile};
my $recipe_name = q{chanjo_sexcheck};
my $sample_id   = q{sample_1};
my %qc_data     = (
    sample => {
        $sample_id => {
            $infile => {
                $recipe_name => { $attribute => q{male}, },
            },
        },
    },
);

my $got_attribute_value = get_qc_data_sample_recipe_attributes(
    {
        attribute    => $attribute,
        infile       => $infile,
        recipe_name  => $recipe_name,
        sample_id    => $sample_id,
        qc_data_href => \%qc_data,
    }
);

## Then return attribute for recipe
is( $got_attribute_value, q{male}, q{Got recipe attribute value} );

## Given a sample_id, recipe when supplied with an infile but no attribute
my %got_attribute = get_qc_data_sample_recipe_attributes(
    {
        infile       => $infile,
        recipe_name  => $recipe_name,
        sample_id    => $sample_id,
        qc_data_href => \%qc_data,
    }
);

## Then return attribute hash for recipes
is_deeply(
    \%got_attribute,
    \%{ $qc_data{sample}{$sample_id}{$infile}{$recipe_name} },
    q{Got recipe attribute hash}
);

# Given a recipe when qc_data attribute is undef
delete $qc_data{sample}{$sample_id}{$infile}{$recipe_name}{$attribute};

my $is_undef = get_qc_data_sample_recipe_attributes(
    {
        attribute    => $attribute,
        infile       => $infile,
        recipe_name  => $recipe_name,
        sample_id    => $sample_id,
        qc_data_href => \%qc_data,
    }
);

## Then return attribute for recipe
is( $is_undef, undef, q{Return undef for missing attribute} );
done_testing();
