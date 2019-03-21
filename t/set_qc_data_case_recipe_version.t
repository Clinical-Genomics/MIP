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
use Modern::Perl qw{ 2014 };
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
        q{MIP::Qc_data}        => [qw{ set_qc_data_case_recipe_version }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qc_data qw{ set_qc_data_case_recipe_version };

diag(   q{Test set_qc_data_case_recipe_version from Qc_data.pm v}
      . $MIP::Qc_data::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a case level recipe and no version
my %qc_data;
my $recipe_name = q{bwa_mem};
set_qc_data_case_recipe_version(
    {
        qc_data_href => \%qc_data,
        recipe_name  => $recipe_name,
        version      => undef,
    }
);

## Then skip setting version in qc_data for recipe
is( $qc_data{recipe}{$recipe_name}{version},
    undef, q{Skip setting case levele recipe version} );

## Given a case level recipe and version
my $version = q{1.0.0};
set_qc_data_case_recipe_version(
    {
        qc_data_href => \%qc_data,
        recipe_name  => $recipe_name,
        version      => $version,
    }
);

## Then set version in qc_data for recipe
is( $qc_data{recipe}{$recipe_name}{version}, $version,
    q{Set case levele recipe version} );

done_testing();
