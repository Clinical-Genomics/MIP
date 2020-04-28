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
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $DATE => q{161011};

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Fastq}          => [qw{ parse_fastq_infiles_format }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Fastq qw{ parse_fastq_infiles_format };

diag(   q{Test parse_fastq_infiles_format from Fastq.pm v}
      . $MIP::Fastq::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given compressed file, when following MIP filename convention
my $file_name = q{1_161011_TestFilev2_ADM1059A1_TCCGGAGA_1.fastq.gz};

## Parse infile according to filename convention
my %infile_info = parse_fastq_infiles_format( { file_name => $file_name, } );

my %expected_features = (
    date             => $DATE,
    direction        => 1,
    flowcell         => q{TestFilev2},
    infile_sample_id => q{ADM1059A1},
    index            => q{TCCGGAGA},
    lane             => 1,
);

## Then return all features from filename
is_deeply( \%infile_info, \%expected_features,
    q{Compressed file follows file convention} );

## Given compressed file, when followinf MIP filename convention
$file_name = q{1_161011_TestFilev2_ADM1059A1_TCCGGAGA_1.fastq};

## When parsing infile according to filename convention
%infile_info = parse_fastq_infiles_format( { file_name => $file_name, } );

## Then return all features from filename
is_deeply( \%infile_info, \%expected_features,
    q{Uncompressed file follows file convention} );

## Given file, when not following MIP filename convention
$file_name = q{TestFilev2_ADM1059A1_TCCGGAGA_1.fastq};

## When parsing infile according to filename convention
my @responces = trap { parse_fastq_infiles_format( { file_name => $file_name, } ) };

## Then return no features from filename
is( $responces[0], undef, q{Return undef for file not following file convention} );
like(
    $trap->stderr,
    qr/Missing \s+ file \s+ name \s+ feature/xms,
    q{Throw warning for missing features}
);

done_testing();
