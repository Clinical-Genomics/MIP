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
use MIP::Test::Fixtures qw{ test_log test_standard_cli };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

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
        q{MIP::File::Format::Dragen} => [qw{ create_dragen_fastq_list_sample_id }],
        q{MIP::Test::Fixtures}       => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Dragen qw{ create_dragen_fastq_list_sample_id };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test create_dragen_fastq_list_sample_id from Dragen.pm v}
      . $MIP::File::Format::Dragen::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create temporary directories and paths
my $test_dir        = File::Temp->newdir();
my $fastq_file_path = catfile( $test_dir, q{test_dragen_fastq.csv} );

my $log = test_log( { no_screen => 1, } );

## Given headers and fastq info
my @fastq_list_headers      = qw{ RGID RGSM RGLB Lane Read1File Read2File };
my @dragen_fastq_list_lines = qw{ agct sample-1 lib-1 lane-1 fastq-1 fastq-2 };

create_dragen_fastq_list_sample_id(
    {
        fastq_list_lines_ref => \@dragen_fastq_list_lines,
        fastq_list_file_path => $fastq_file_path,
        log                  => $log,
    }
);

## Then dragen fastq list file should have been created
ok( -e $fastq_file_path, q{Created dragen fastq list file} );

## Given a dragen fastq list file
# Read the created file to check content

# Create anonymous filehandle
my $filehandle = IO::Handle->new();

open $filehandle, q{<}, $fastq_file_path
  or $log->logdie(qq{Can't open $fastq_file_path: $ERRNO });

my @lines;
LINE:
while (<$filehandle>) {

    chomp;
    push @lines, split $COMMA;
}
close $filehandle;
my @expected_lines = ( @fastq_list_headers, @dragen_fastq_list_lines );
is_deeply( \@lines, \@expected_lines, q{Correctly wrote file content} );

done_testing();
