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
        q{MIP::Qcc_regexp}     => [qw{ get_qcc_regexp_recipe_attribute }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qcc_regexp qw{ get_qcc_regexp_recipe_attribute };

diag(   q{Test get_qcc_regexp_recipe_attribute from Qcc_regexp.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a recipe when called with attribute
my $attribute       = q{version};
my %qcc_regexp_href = ( fastqc => { version => q{a_version}, } );
my $recipe_name     = q{fastqc};

my $got_attribute = get_qcc_regexp_recipe_attribute(
    {
        attribute       => $attribute,
        qcc_regexp_href => \%qcc_regexp_href,
        recipe_name     => $recipe_name,
    }
);

## Then return attribute for recipe
is( $got_attribute, q{a_version}, q{Got recipe attribute} );

## Given a recipe when called without attribute
my %attribute = get_qcc_regexp_recipe_attribute(
    {
        qcc_regexp_href => \%qcc_regexp_href,
        recipe_name     => $recipe_name,
    }
);

my %expected_attribute = ( version => q{a_version} );

## Then return attribute hash for recipe
is_deeply( \%attribute, \%expected_attribute, q{Got recipe attribute hash} );

done_testing();
