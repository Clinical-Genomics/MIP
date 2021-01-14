#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp qw{ tempdir };
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
use MIP::Constants qw{ $COMMA $DOT $SPACE $UNDERSCORE };
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Script::Setup_script} => [qw{ build_script_directories_and_paths }],
        q{MIP::Test::Fixtures}       => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Script::Setup_script qw{ build_script_directories_and_paths };

diag(   q{Test build_script_directories_and_paths from Setup_script.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 1, } );
my $test_dir = tempdir( CLEANUP => 1, );

## Given recipe dirs and recipe name
my $directory_id  = q{dir};
my $outdata_dir   = catdir( $test_dir, q{outdata_dir} );
my $outscript_dir = catdir( $test_dir, q{outscript_dir} );
my $recipe_data_directory_path;
my $recipe_directory = q{recipe_dir};
my $recipe_mode      = 1;
my $recipe_name      = q{bwa_mem};

## When recipe mode is one

( my ( $file_info_path, $file_path_prefix ), $recipe_data_directory_path ) =
  build_script_directories_and_paths(
    {
        directory_id               => $directory_id,
        outdata_dir                => $outdata_dir,
        outscript_dir              => $outscript_dir,
        recipe_data_directory_path => $recipe_data_directory_path,
        recipe_directory           => $recipe_directory,
        recipe_mode                => $recipe_mode,
        recipe_name                => $recipe_name,
    }
  );

my $expected_file_info_path =
  catfile( $outdata_dir, $directory_id, $recipe_directory, q{info},
    $recipe_name . $UNDERSCORE . $directory_id . $DOT );

## Then return the file info path to write stdout and stderr to
is( $file_info_path, $expected_file_info_path, q{Built file_info_path} );

my $expected_file_path_prefix = catfile( $outscript_dir, $directory_id, $recipe_directory,
    $recipe_name . $UNDERSCORE . $directory_id . $DOT );

## Then return the file path prefix to write executable scripts to
is( $file_path_prefix, $expected_file_path_prefix, q{Built file_path_prefix} );

my $expected_recipe_directory = catdir( $outdata_dir, $directory_id, $recipe_directory );

## Then return the directory to write data to
is( $recipe_data_directory_path, $expected_recipe_directory,
    q{Built recipe_data_directory_path} );

## When recipe mode is two
$recipe_mode = 2;

( $file_info_path, $file_path_prefix, $recipe_data_directory_path ) =
  build_script_directories_and_paths(
    {
        directory_id               => $directory_id,
        outdata_dir                => $outdata_dir,
        outscript_dir              => $outscript_dir,
        recipe_data_directory_path => $recipe_data_directory_path,
        recipe_directory           => $recipe_directory,
        recipe_mode                => $recipe_mode,
        recipe_name                => $recipe_name,
    }
  );

$expected_file_info_path =
  catfile( $outdata_dir, $directory_id, $recipe_directory, q{info},
    q{dry_run} . $UNDERSCORE . $recipe_name . $UNDERSCORE . $directory_id . $DOT );

## Then return the file info path to write stdout and stderr to
is( $file_info_path, $expected_file_info_path, q{Built dry run file_info_path} );

$expected_file_path_prefix = catfile( $outscript_dir, $directory_id, $recipe_directory,
    q{dry_run} . $UNDERSCORE . $recipe_name . $UNDERSCORE . $directory_id . $DOT );

## Then return the file path prefix to write executable scripts to
is( $file_path_prefix, $expected_file_path_prefix, q{Built dry run file_path_prefix} );

done_testing();
