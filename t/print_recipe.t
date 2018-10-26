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
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

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

my %parameter = ( bwa_mem => { type => q{recipe} } );

my @printed_recipes = print_recipe(
    {
        parameter_href    => \%parameter,
        print_recipe_mode => 1,
        define_parameters_files_ref =>
          [ catfile( $Bin, qw{ data test_data define_parameters.yaml } ) ],
    }
);

is( scalar @printed_recipes, 1, q{Did not print rio block: bamcalibrationblock} );

my @recipe_mode = split $SPACE, $printed_recipes[0];

is( $recipe_mode[1], 1, q{Printed correct recipe mode} );

is( $printed_recipes[0], q{bwa_mem 1}, q{Printed recipe} );

done_testing();
