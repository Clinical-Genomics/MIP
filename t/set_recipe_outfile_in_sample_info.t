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
        q{MIP::Sample_info}    => [qw{ set_recipe_outfile_in_sample_info }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw(set_recipe_outfile_in_sample_info);

diag(   q{Test set_recipe_outfile_in_sample_info from Sample_info.pm v}
      . $MIP::Sample_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Init hash
my %sample_info;

# Test variables
my $test_recipe_name = q{test_recipe};
my $directory        = q{test_directory};
my $outfile          = q{test.yaml};
my $path             = catfile( $directory, $outfile );
my $version          = q{1.0.1};

## Family level
set_recipe_outfile_in_sample_info(
    {
        sample_info_href => \%sample_info,
        recipe_name      => $test_recipe_name,
        outdirectory     => $directory,
        outfile          => $outfile,
        path             => $path,
        version          => $version,
    }
);

## Test
is( exists $sample_info{recipe}{$test_recipe_name}, 1, q{Created case level hash key} );

is( $sample_info{recipe}{$test_recipe_name}{outdirectory},
    $directory, q{Assigned correct value to case level outdirectory} );

is( $sample_info{recipe}{$test_recipe_name}{outfile},
    $outfile, q{Assigned correct value to case level outfile} );

is( $sample_info{recipe}{$test_recipe_name}{path},
    $path, q{Assigned correct value to case level path} );

is( $sample_info{recipe}{$test_recipe_name}{version},
    $version, q{Assigned correct value to case level version} );

## Sample level, without infile
my $sample_id = q{test_sample_id};
my $infile    = q{test_infile};

set_recipe_outfile_in_sample_info(
    {
        sample_info_href => \%sample_info,
        sample_id        => $sample_id,
        recipe_name      => $test_recipe_name,
        outdirectory     => $directory,
        outfile          => $outfile,
        path             => $path,
        version          => $version,
    }
);
## Test

is( $sample_info{sample}{$sample_id}{recipe}{$test_recipe_name}{path},
    $path, q{Assigned correct value to sample level path} );

my %test_no_infile = (
    outdirectory => q{Value to sample level outdirectory not assigned},
    outfile      => q{Value to sample level outfile not assigned},
    path         => q{Value to sample level path not assigned},
    version      => q{Value to sample level version not assigned},
);

while ( my ( $parameter, $test_comment ) = each %test_no_infile ) {
    my $test_result =
      $sample_info{sample}{$sample_id}{recipe}{$test_recipe_name}{$infile}{$parameter};
    is( $test_result, undef, $test_comment );
}

## Sample level, with infile
set_recipe_outfile_in_sample_info(
    {
        sample_info_href => \%sample_info,
        sample_id        => $sample_id,
        infile           => $infile,
        recipe_name      => $test_recipe_name,
        outdirectory     => $directory,
        outfile          => $outfile,
        path             => $path,
        version          => $version,
    }
);

## Test
is( exists $sample_info{sample}{$sample_id}{recipe}{$test_recipe_name}{$infile},
    1, q{Created sample level hash key} );

is( $sample_info{sample}{$sample_id}{recipe}{$test_recipe_name}{$infile}{outdirectory},
    $directory, q{Assigned correct value to sample level outdirectory} );

is( $sample_info{sample}{$sample_id}{recipe}{$test_recipe_name}{$infile}{outfile},
    $outfile, q{Assigned correct value to sample level outfile} );

is( $sample_info{sample}{$sample_id}{recipe}{$test_recipe_name}{$infile}{path},
    $path, q{Assigned correct value to sample level path} );

is( $sample_info{sample}{$sample_id}{recipe}{$test_recipe_name}{$infile}{version},
    $version, q{Assigned correct value to sample level version} );

done_testing();
