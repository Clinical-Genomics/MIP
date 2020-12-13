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
use MIP::Constants qw{ $COLON $COMMA $SPACE $UNDERSCORE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Qc_data}        => [qw{ add_qc_data_evaluation_info }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qc_data qw{ add_qc_data_evaluation_info };

diag(   q{Test add_qc_data_evaluation_info from Qc_data.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given value to add to qc data hash
my $metric = q{duplicates};
my %qc_data;
my $qc_metric_value = 1;
my $recipe_name     = q{plink_relation_check};
my $value           = $recipe_name . $UNDERSCORE . $metric . $COLON . $qc_metric_value;

add_qc_data_evaluation_info(
    {
        qc_data_href => \%qc_data,
        recipe_name  => $recipe_name,
        value        => $value,
    }
);

## Then values should be added to hash on evaluation level
is_deeply( \@{ $qc_data{evaluation}{$recipe_name} },
    [$value], q{Added values for evaluation level} );

done_testing();
