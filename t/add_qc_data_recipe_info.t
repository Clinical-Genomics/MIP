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
        q{MIP::Qc_data}        => [qw{ add_qc_data_recipe_info }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qc_data qw{ add_qc_data_recipe_info };

diag(   q{Test add_qc_data_recipe_info from Qc_data.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given key and value to add to qc data hash when sample and infile
my $key    = q{sample_sexcheck};
my $infile = q{an_infile};
my %qc_data;
my $recipe_name = q{plink_gender_check};
my $sample_id   = q{ADM1059A1};
my $value       = q{ADM1059A1:PASS};

add_qc_data_recipe_info(
    {
        key          => $key,
        infile       => $infile,
        qc_data_href => \%qc_data,
        recipe_name  => $recipe_name,
        sample_id    => $sample_id,
        value        => $value,
    }
);

## Then values should be added to hash on sample level
is_deeply( \@{ $qc_data{sample}{$sample_id}{$infile}{$recipe_name}{$key} },
    [$value], q{Added values on sample level} );

## Given key and value to add when case level
add_qc_data_recipe_info(
    {
        key          => $key,
        qc_data_href => \%qc_data,
        recipe_name  => $recipe_name,
        value        => $value,
    }
);

## Then values should be added to hash on sample level
is_deeply( \@{ $qc_data{recipe}{$recipe_name}{$key} },
    [$value], q{Added values on case level} );

done_testing();
