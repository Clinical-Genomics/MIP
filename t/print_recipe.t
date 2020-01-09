#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
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
use Test::Trap;

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
        q{MIP::Get::Analysis}  => [qw{ print_recipe }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::Analysis qw{ print_recipe };

diag(   q{Test print_recipe from Analysis.pm v}
      . $MIP::Get::Analysis::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given the option to not print recipes
my %parameter = ( bwa_mem => { type => q{recipe} } );
my $return    = print_recipe(
    {
        define_parameters_files_ref =>
          [ catfile( $Bin, qw{ data test_data define_parameters.yaml } ) ],
        parameter_href    => \%parameter,
        print_recipe_mode => 1,
    }
);

## Do not print
is( $return, undef, q{Do not print} );

## Given a recipe and to print
my %active_parameter = ( print_recipe => 1 );

trap {
    print_recipe(
        {
            define_parameters_files_ref =>
              [ catfile( $Bin, qw{ data test_data define_parameters.yaml } ) ],
            parameter_href    => \%parameter,
            print_recipe      => $active_parameter{print_recipe},
            print_recipe_mode => 1,
        }
    )
};

## Then write bwa_mem recipe and mode
my ( $recipe, $recipe_mode ) = $trap->stdout =~ qr{ (bwa_mem)\s+(1)}sxm;
is( $recipe, q{bwa_mem}, q{Printed recipe} );

is( $recipe_mode, 1, q{Printed correct recipe mode} );

done_testing();
