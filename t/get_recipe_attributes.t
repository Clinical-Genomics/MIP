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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Get::Parameter} => [qw{ get_recipe_attributes }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::Parameter qw{ get_recipe_attributes };

diag(   q{Test get_recipe_attributes from Parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a program parameter
my %parameter   = test_mip_hashes( { mip_hash_name => q{define_parameter}, } );
my $recipe_name = q{bwa_mem};

my %rec_atr = get_recipe_attributes(
    {
        parameter_href => \%parameter,
        recipe_name    => $recipe_name,
    }
);

## Then return all program attributes
is_deeply( \%{ $parameter{$recipe_name} }, \%rec_atr, q{Got attributes} );

## Given a program parameter attribute
my $attribute     = q{outfile_suffix};
my $got_attribute = get_recipe_attributes(
    {
        attribute      => $attribute,
        parameter_href => \%parameter,
        recipe_name    => $recipe_name,
    }
);

## Then return single attribute
is( $got_attribute, q{.bam}, q{Got attribute} );

done_testing();
