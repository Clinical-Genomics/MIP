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
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Sample_info} => [qw{ set_processing_metafile_in_sample_info }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ set_processing_metafile_in_sample_info };

diag(   q{Test set_processing_metafile_in_sample_info from Sample_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Init hash
my %sample_info;

# Test variables
my $metafile_tag = q{most_complete_bam};
my $metafile     = q{test_metafile};
my $directory    = q{test_directory};
my $path         = catfile( $directory, $metafile );

## Family level
set_processing_metafile_in_sample_info(
    {
        metafile_tag     => $metafile_tag,
        path             => $path,
        sample_info_href => \%sample_info,
    }
);

## Test
is( exists $sample_info{$metafile_tag}, 1, q{Created case level hash key} );

is( $sample_info{$metafile_tag}{path},
    $path, q{Assigned correct value to case level path} );

## Sample level
my $sample_id = q{test_sample_id};

set_processing_metafile_in_sample_info(
    {
        metafile_tag     => $metafile_tag,
        path             => $path,
        sample_id        => $sample_id,
        sample_info_href => \%sample_info,
    }
);

## Test
is( exists $sample_info{sample}{$sample_id}{$metafile_tag},
    1, q{Created sample level hash key} );

is( $sample_info{sample}{$sample_id}{$metafile_tag}{path},
    $path, q{Assigned correct value to sample level path} );

done_testing();
