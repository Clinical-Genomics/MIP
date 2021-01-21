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
use MIP::Constants qw{ $COMMA $FORWARD_SLASH $SPACE };

## Constants
Readonly my $CHR_23 => q{23};

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::File_info} => [qw{ set_io_files }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ set_io_files };

diag(   q{Test set_io_files from File_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a chain id
my $chain_id = 1;

## Given contigs and contig_files
my @contigs = ( 1 .. $CHR_23, qw{ X Y MT } );
my @contig_files =
  map { catfile( qw{ a temp dir }, q{file.} . $_ . q{.bam} ) } @contigs;

## Given a file_info hash
my %file_info;

## Given file contig hash with contig keys and file paths
my %file_contig;
CONTIG:
while ( my ( $index, $contig ) = each @contigs ) {

    $file_contig{$contig} = $contig_files[$index];
}

## Given file_paths and suffixes
my @file_paths =
  ( catfile(qw{ a test dir file_1.1.txt}), catfile(qw{ a test dir file_2.txt}), );
my @file_path_prefixes =
  ( catfile(qw{ a test dir file_1}), catfile(qw{ a test dir file_2}), );
my @file_suffixes = qw{.1.txt .txt};

## Given temporary file paths
my @temp_file_paths =
  ( catfile(qw{ a temp dir file_1.1.txt}), catfile(qw{ a temp dir file_2.txt}), );
my @temp_file_path_prefixes =
  ( catfile(qw{ a temp dir file_1}), catfile(qw{ a temp dir file_2}), );

## Given a sample id
my $id = q{sample_1};

## Given a input stream direction
my $stream = q{in};

## Given a temp dir
my $temp_directory = catfile(qw{ a temp dir});

## Given a recipe name
my $recipe_name = q{bwa_mem};

## When setting io files
set_io_files(
    {
        chain_id       => $chain_id,
        id             => $id,
        file_paths_ref => \@file_paths,
        file_info_href => \%file_info,
        recipe_name    => $recipe_name,
        stream         => $stream,
        temp_directory => $temp_directory,
    }
);

## Then set io features for $chain_id and input stream
is(
    $file_info{io}{$chain_id}{$id}{$recipe_name}{$stream}{dir_path},
    catfile(qw{ a test dir }) . $FORWARD_SLASH,
    q{Set file path}
);

is(
    $file_info{io}{$chain_id}{$id}{$recipe_name}{$stream}{dir_path_prefix},
    catfile(qw{ a test dir }),
    q{Set dir path prefix}
);

is_deeply(
    $file_info{io}{$chain_id}{$id}{$recipe_name}{$stream}{file_names},
    [qw{ file_1.1.txt file_2.txt }],
    q{Set file name}
);

is_deeply( \@{ $file_info{io}{$chain_id}{$id}{$recipe_name}{$stream}{file_name_prefixes} },
    [qw{ file_1 file_2 }], q{Set file name prefixes} );

is_deeply( \@{ $file_info{io}{$chain_id}{$id}{$recipe_name}{$stream}{file_paths} },
    \@file_paths, q{Set file paths} );

is_deeply( \@{ $file_info{io}{$chain_id}{$id}{$recipe_name}{$stream}{file_path_prefixes} },
    \@file_path_prefixes, q{Set file path prefixes} );

is_deeply( \@{ $file_info{io}{$chain_id}{$id}{$recipe_name}{$stream}{file_suffixes} },
    \@file_suffixes, q{Set file suffixes} );

is( $file_info{io}{$chain_id}{$id}{$recipe_name}{$stream}{file_suffix},
    q{.txt}, q{Set file suffix} );

is( $file_info{io}{$chain_id}{$id}{$recipe_name}{$stream}{file_constant_suffix},
    undef, q{Did not set file constant suffix} );

## Then set io features for $chain_id and temp stream
my $temp_stream = q{temp};

is(
    $file_info{io}{$chain_id}{$id}{$recipe_name}{$temp_stream}{dir_path},
    catfile(qw{ a temp dir }) . $FORWARD_SLASH,
    q{Set file path for temp stream}
);

is(
    $file_info{io}{$chain_id}{$id}{$recipe_name}{$temp_stream}{dir_path_prefix},
    catfile(qw{ a temp dir }),
    q{Set dir path prefix for temp stream}
);

is_deeply(
    $file_info{io}{$chain_id}{$id}{$recipe_name}{$temp_stream}{file_names},
    [qw{ file_1.1.txt file_2.txt }],
    q{Set file name for temp stream}
);

is_deeply( \@{ $file_info{io}{$chain_id}{$id}{$recipe_name}{$temp_stream}{file_name_prefixes} },
    [qw{ file_1 file_2 }], q{Set file name prefixes for temp stream} );

is_deeply( \@{ $file_info{io}{$chain_id}{$id}{$recipe_name}{$temp_stream}{file_paths} },
    \@temp_file_paths, q{Set file paths for temp stream} );

is_deeply( \@{ $file_info{io}{$chain_id}{$id}{$recipe_name}{$temp_stream}{file_path_prefixes} },
    \@temp_file_path_prefixes, q{Set file path prefixes for temp stream} );

is_deeply( \@{ $file_info{io}{$chain_id}{$id}{$recipe_name}{$temp_stream}{file_suffixes} },
    \@file_suffixes, q{Set file suffixes for temp stream} );

is( $file_info{io}{$chain_id}{$id}{$recipe_name}{$temp_stream}{file_suffix},
    q{.txt}, q{Set file suffix for temp stream} );

is( $file_info{io}{$chain_id}{$id}{$recipe_name}{$temp_stream}{file_constant_suffix},
    undef, q{Did not set file constant suffix for temp stream} );

## Given constant file features
@file_paths =
  ( catfile(qw{ a test dir file.fastq.gz}), catfile(qw{ a test dir file.fastq.gz}), );

## When setting io files
set_io_files(
    {
        chain_id       => $chain_id,
        id             => $id,
        file_paths_ref => \@file_paths,
        file_info_href => \%file_info,
        recipe_name    => $recipe_name,
        stream         => $stream,
        temp_directory => $temp_directory,
    }
);

## Then set file constant suffix for chain_id and stream
is( $file_info{io}{$chain_id}{$id}{$recipe_name}{$stream}{file_constant_suffix},
    q{.fastq.gz}, q{Set file constant suffix} );

is( $file_info{io}{$chain_id}{$id}{$recipe_name}{$temp_stream}{file_constant_suffix},
    q{.fastq.gz}, q{Set file constant suffix for temp stream} );

is( $file_info{io}{$chain_id}{$id}{$recipe_name}{$stream}{file_name_prefix},
    q{file}, q{Set file constant file name prefix} );

is( $file_info{io}{$chain_id}{$id}{$recipe_name}{$temp_stream}{file_name_prefix},
    q{file}, q{Set file constant file name prefix for temp stream} );

is(
    $file_info{io}{$chain_id}{$id}{$recipe_name}{$stream}{dir_path_prefix},
    catfile(qw{ a test dir }),
    q{Set file constant dir path}
);

is(
    $file_info{io}{$chain_id}{$id}{$recipe_name}{$temp_stream}{dir_path_prefix},
    catfile(qw{ a temp dir }),
    q{Set file constant dir path for temp stream}
);

## Given contig files

## When setting io files using contig files
set_io_files(
    {
        chain_id       => $chain_id,
        id             => $id,
        file_paths_ref => \@contig_files,
        file_info_href => \%file_info,
        recipe_name    => $recipe_name,
        stream         => $stream,
        temp_directory => $temp_directory,
    }
);

## Then set file_path_href to file_contigs hash
is_deeply( \%{ $file_info{io}{$chain_id}{$id}{$recipe_name}{$stream}{file_path_href} },
    \%{file_contig}, q{Set file path hash} );
done_testing();
