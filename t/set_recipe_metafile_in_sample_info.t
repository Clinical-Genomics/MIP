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
        q{MIP::Sample_info}    => [qw{ set_recipe_metafile_in_sample_info }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ set_recipe_metafile_in_sample_info };

diag(   q{Test set_recipe_metafile_in_sample_info from Sample_info.pm v}
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
my $metafile         = q{test_metafile};
my $directory        = q{test_directory};
my $file             = q{test.yaml};
my $path             = catfile( $directory, $file );
my $version          = q{1.0.1};
my $processed_by     = q{picard_markduplicates};

## Family level
set_recipe_metafile_in_sample_info(
    {
        sample_info_href => \%sample_info,
        recipe_name      => $test_recipe_name,
        metafile_tag     => $metafile,
        directory        => $directory,
        file             => $file,
        path             => $path,
        version          => $version,
        processed_by     => $processed_by,
    }
);

## Test
is( exists $sample_info{recipe}{$test_recipe_name}{$metafile},
    1, q{Created case level hash key} );

is( $sample_info{recipe}{$test_recipe_name}{$metafile}{directory},
    $directory, q{Assigned correct value to case level directory} );

is( $sample_info{recipe}{$test_recipe_name}{$metafile}{file},
    $file, q{Assigned correct value to case level file} );

is( $sample_info{recipe}{$test_recipe_name}{$metafile}{path},
    $path, q{Assigned correct value to case level path} );

is( $sample_info{recipe}{$test_recipe_name}{$metafile}{version},
    $version, q{Assigned correct value to case level version} );

is( $sample_info{recipe}{$test_recipe_name}{$metafile}{processed_by},
    $processed_by, q{Assigned correct value to case level processed_by} );

## Sample level
my $sample_id = q{test_sample_id};
my $infile    = q{test_infile};

set_recipe_metafile_in_sample_info(
    {
        sample_info_href => \%sample_info,
        sample_id        => $sample_id,
        infile           => $infile,
        recipe_name      => $test_recipe_name,
        metafile_tag     => $metafile,
        directory        => $directory,
        file             => $file,
        path             => $path,
        version          => $version,
        processed_by     => $processed_by,
    }
);

## Test
is(
    exists $sample_info{sample}{$sample_id}{recipe}{$test_recipe_name}
      {$infile}{$metafile},
    1,
    q{Created sample level hash key}
);

is(
    $sample_info{sample}{$sample_id}{recipe}{$test_recipe_name}
      {$infile}{$metafile}{directory},
    $directory, q{Assigned correct value to sample level directory}
);

is(
    $sample_info{sample}{$sample_id}{recipe}{$test_recipe_name}{$infile}{$metafile}{file},
    $file, q{Assigned correct value to sample level file}
);

is(
    $sample_info{sample}{$sample_id}{recipe}{$test_recipe_name}{$infile}{$metafile}{path},
    $path, q{Assigned correct value to sample level path}
);

is(
    $sample_info{sample}{$sample_id}{recipe}{$test_recipe_name}
      {$infile}{$metafile}{version},
    $version, q{Assigned correct value to sample level version}
);

is(
    $sample_info{sample}{$sample_id}{recipe}{$test_recipe_name}
      {$infile}{$metafile}{processed_by},
    $processed_by, q{Assigned correct value to sample level processed_by}
);

done_testing();
