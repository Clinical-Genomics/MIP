#! /usr/bin/env perl

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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $EMPTY_STR $SPACE };
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
        q{MIP::Environment::Container} => [qw{ get_recipe_executable_bind_path }],
        q{MIP::Test::Fixtures}         => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Container qw{ get_recipe_executable_bind_path };

diag(   q{Test get_recipe_executable_bind_path from Container.pm v}
      . $MIP::Environment::Container::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a recipe name
my $recipe_name = q{bwa_mem};

## Given a recipe bind path
my %active_parameter = (
    recipe_bind_path => { $recipe_name => [qw{ a_dir }], },
    temp_directory   => q{a_temp_dir},
);

## Given a cache and recipe executables
my %parameter = (
    cache   => { recipe              => [$recipe_name], },
    bwa_mem => { program_executables => [qw{ bwa }], },
);

## When getting recipe executable bind paths
my %recipe_executable_bind_path = get_recipe_executable_bind_path(
    {
        active_parameter_href => \%active_parameter,
        parameter_href        => \%parameter,
    }
);

my $xdg_runtime_dir =
  q{a_temp_dir} . $COLON . catfile( $EMPTY_STR, qw{ run user }, q{$(id -u)} );
my %expected_recipe_executable_bind_paths = ( bwa => [ qw{ a_dir }, $xdg_runtime_dir ] );

## Then return recipe executable bind paths
is_deeply(
    \%recipe_executable_bind_path,
    \%expected_recipe_executable_bind_paths,
    q{Got recipe executable bind paths}
);

done_testing();
