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
        q{MIP::Qc_data}        => [qw{ get_qc_data_case_recipe_attributes }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qc_data qw{ get_qc_data_case_recipe_attributes };

diag(   q{Test get_qc_data_case_recipe_attributes from Qc_data.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a recioe with attributes array values
my @metrics     = qw{ ADM1059A1:2 ADM1059A2:1 ADM1059A3:0 };
my $recipe_name = q{plink_sexcheck};
my %qc_data = ( recipe => { plink_sexcheck => { sample_sexcheck => [@metrics], }, }, );

my $got_attributes_ref = get_qc_data_case_recipe_attributes(
    {
        attribute    => q{sample_sexcheck},
        qc_data_href => \%qc_data,
        recipe_name  => $recipe_name,
    }
);

## Then return attributes_ref
is_deeply( \@{$got_attributes_ref}, \@metrics, q{Got attributes_ref} );

## Given a recipe when calling with no attributes
my %got_attribute = get_qc_data_case_recipe_attributes(
    {
        qc_data_href => \%qc_data,
        recipe_name  => $recipe_name,
    }
);

## Then return attribute(s_ref)
is_deeply(
    \%got_attribute,
    \%{ $qc_data{recipe}{plink_sexcheck} },
    q{Got attribute_href}
);

done_testing();
