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
        q{MIP::Qc_data}        => [qw{ get_qc_recipe_data }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qc_data qw{ get_qc_recipe_data };

diag(   q{Test get_qc_recipe_data from Qc_data.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a recipe when supplied with an header element index and attribute
my $attribute      = q{gender};
my $recipe_name    = q{chanjo_sexcheck};
my %qc_recipe_data = ( $recipe_name => { $attribute => [qw{ male }], }, );

my $got_attribute_value = get_qc_recipe_data(
    {
        attribute           => $attribute,
        recipe_name         => $recipe_name,
        qc_recipe_data_href => \%qc_recipe_data,
        qc_header_index     => 0,
    }
);

## Then return attribute for recipe
is( $got_attribute_value, q{male}, q{Got qc recipe data attribute value} );

## Given a recipe when supplied with an attribute
my @got_attributes = get_qc_recipe_data(
    {
        attribute           => $attribute,
        recipe_name         => $recipe_name,
        qc_recipe_data_href => \%qc_recipe_data,
    }
);

## Then return attribute array for recipe
is_deeply(
    \@got_attributes,
    \@{ $qc_recipe_data{$recipe_name}{$attribute} },
    q{Got qc recipe data attribute array}
);

## Given a recipe when supplied with no attribute and header element index
my %got_attribute = get_qc_recipe_data(
    {
        recipe_name         => $recipe_name,
        qc_recipe_data_href => \%qc_recipe_data,
    }
);

## Then return attribute array for recipe
is_deeply(
    \%got_attribute,
    \%{ $qc_recipe_data{$recipe_name} },
    q{Got qc recipe data attribute hash}
);
done_testing();
