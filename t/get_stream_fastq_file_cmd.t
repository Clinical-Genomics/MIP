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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Fastq} => [qw{ get_stream_fastq_file_cmd }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Fastq qw{ get_stream_fastq_file_cmd };

diag(   q{Test get_stream_fastq_file_cmd from Fastq.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given an infile dir
my $infile_dir = catdir( $Bin, qw{ data 643594-miptest test_data ADM1059A1 fastq  } );

## Given a compressed and uncompressed fastq file
my @infiles =
  qw{ 7_161011_HHJJCCCXY_ADM1059A1_NAATGCGC_1.fastq 7_161011_HHJJCCCXY_ADM1059A1_NAATGCGC_1.fastq.gz };

## Given fastq file paths
my @fastq_files = map { catfile( $infile_dir, $_ ) } @infiles;

## When building the command to stream a chunk from each fastq file
my @read_file_chunk_cmds = get_stream_fastq_file_cmd( { fastq_files_ref => \@fastq_files, } );

## Then commands streaming a chunk of uncompressed fastq file is returned
like( $read_file_chunk_cmds[0], qr/[<(]cat/sxm, q{Got cmd for uncompressed fastq file} );

## Then commands streaming a chunk of compressed fastq files is returned
like( $read_file_chunk_cmds[1], qr/[<(]gzip/sxm, q{Got cmd for compressed fastq file} );

done_testing();
