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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Qc_data}        => [qw{ set_qc_data_recipe_info }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qc_data qw{ set_qc_data_recipe_info };

diag(   q{Test set_qc_data_recipe_info from Qc_data.pm v}
      . $MIP::Qc_data::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a key value pair
my $key    = q{greeting};
my $infile = q{an_infile};
my %qc_data;
my $recipe_name = q{japan};
my $sample_id   = q{sample_1};
my $value       = q{konnichi wa};

set_qc_data_recipe_info(
    {
        key          => $key,
        qc_data_href => \%qc_data,
        recipe_name  => $recipe_name,
        value        => $value,
    }
);

## Then set info on case level
is( $qc_data{recipe}{$recipe_name}{$key}, $value, q{Set case level recipe info} );

## Given a sample id and infile
set_qc_data_recipe_info(
    {
        key          => $key,
        infile       => $infile,
        qc_data_href => \%qc_data,
        recipe_name  => $recipe_name,
        sample_id    => $sample_id,
        value        => $value,
    }
);

## Then set info for sample id and infile level
is( $qc_data{sample}{$sample_id}{$infile}{$recipe_name}{$key},
    $value, q{Set sample id and infile level recipe info} );

## Given a sample id and recipe name
set_qc_data_recipe_info(
    {
        key          => $key,
        qc_data_href => \%qc_data,
        recipe_name  => $recipe_name,
        sample_id    => $sample_id,
        value        => $value,
    }
);

## Then set info for sample id and recipe level
is( $qc_data{sample}{$sample_id}{$recipe_name}{$key},
    $value, q{Set sample id and recipe level recipe info} );

## Given a sample id and key
set_qc_data_recipe_info(
    {
        key          => $key,
        qc_data_href => \%qc_data,
        sample_id    => $sample_id,
        value        => $value,
    }
);

## Then set info for sample id and recipe level
is( $qc_data{sample}{$sample_id}{$key},
    $value, q{Set sample id and key level recipe info} );

done_testing();
