#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
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
use MIP::Constants qw{ $COMMA $SPACE $UNDERSCORE };
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
        q{MIP::Fastq}          => [qw{ define_mip_fastq_file_features }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Fastq qw{ define_mip_fastq_file_features };

diag(   q{Test define_mip_fastq_file_features from Fastq.pm v}
      . $MIP::Fastq::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given fastq file name info
my $date               = q{150703};
my $direction          = 1;
my $flowcell           = q{Undetermined-flow-rider};
my $index              = q{ATCG};
my $lane               = 1;
my $original_file_name = q{file-1};
my $sample_id          = q{sample-1};

## Define file formats
my ( $mip_file_format, $mip_file_format_with_direction,
    $original_file_name_prefix, $run_barcode )
  = define_mip_fastq_file_features(
    {
        date               => $date,
        direction          => $direction,
        flowcell           => $flowcell,
        index              => $index,
        lane               => $lane,
        original_file_name => $original_file_name,
        sample_id          => $sample_id,
    }
  );
my $expected_file_format = join $UNDERSCORE,
  ( $sample_id, $date, $flowcell, $index, q{lane} . $lane );
my $expected_file_format_with_direction = join $UNDERSCORE,
  ( $sample_id, $date, $flowcell, $index, q{lane} . $lane, $direction );
my $expected_original_file_name_prefix = substr $original_file_name, 0,
  index $original_file_name, q{.fastq};
my $expected_run_barcode = join $UNDERSCORE, ( $date, $flowcell, $lane, $index );

## Then mip file format should be return
is( $mip_file_format, $expected_file_format, q{Define mip file format} );

## Then mip file format including direction should be return
is(
    $mip_file_format_with_direction,
    $expected_file_format_with_direction,
    q{Define mip file format with direction}
);

## Then mip original file name prefix should be return
is(
    $original_file_name_prefix,
    $expected_original_file_name_prefix,
    q{Define original file name prefix}
);

## Then mip run barcode should be return
is( $run_barcode, $expected_run_barcode, q{Define run bar code} );

done_testing();
