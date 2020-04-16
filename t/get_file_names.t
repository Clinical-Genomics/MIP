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
        q{MIP::File::Path}     => [qw{ get_file_names }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Path qw{ get_file_names };

diag(   q{Test get_file_names from File.pm v}
      . $MIP::File::Path::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given an infile directory, when applying rules
my $infile_directory =
  catfile( $Bin, qw{ data 643594-miptest test_data ADM1059A1 fastq } );

my @infiles = get_file_names(
    {
        file_directory   => $infile_directory,
        rule_name        => q{*.fastq*},
        rule_skip_subdir => q{original_fastq_files},
    }
);

my @expected_files = qw{ 1_161011_TestFilev2_ADM1059A1_TCCGGAGA_1.fastq.gz
  1_161011_TestFilev2_ADM1059A1_TCCGGAGA_2.fastq.gz
  2_161011_TestFilev2-Interleaved_ADM1059A1_TCCGGAGA_1.fastq.gz
  2_161011_TestFilev2_ADM1059A1_TCCGGAGA_1.fastq.gz
  7_161011_HHJJCCCXY_ADM1059A1_NAATGCGC_1.fastq.gz
  8_161011_HHJJCCCXY_ADM1059A1_NAATGCGC_1.fastq.gz };

## Then skip sub dir
is_deeply( \@infiles, \@expected_files, q{Found all files when skipping sub dir} );

## Given an infile directory, when applying rule file name
@infiles = get_file_names(
    {
        file_directory => $infile_directory,
        rule_name      => q{*.fastq*},
    }
);

push @expected_files, q{test.fastq.gz};

## Then find all files recursively
is_deeply( \@infiles, \@expected_files, q{Found all files recursively} );

## Given an infile directory, when applying no rules
@infiles = get_file_names( { file_directory => $infile_directory, } );

my @expected_file_objects = qw{ fastq
  1_161011_TestFilev2_ADM1059A1_TCCGGAGA_1.fastq.gz
  1_161011_TestFilev2_ADM1059A1_TCCGGAGA_2.fastq.gz
  2_161011_TestFilev2-Interleaved_ADM1059A1_TCCGGAGA_1.fastq.gz
  2_161011_TestFilev2_ADM1059A1_TCCGGAGA_1.fastq.gz
  7_161011_HHJJCCCXY_ADM1059A1_NAATGCGC_1.fastq.gz
  8_161011_HHJJCCCXY_ADM1059A1_NAATGCGC_1.fastq.gz
  original_fastq_files
  test.fastq.gz
  test.txt };

## Then find all files and dir recursively
is_deeply( \@infiles, \@expected_file_objects, q{Found all files and dirs recursively} );

done_testing();
