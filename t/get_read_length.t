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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Fastq}          => [qw{ get_read_length }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Fastq qw{ get_read_length };

diag(   q{Test get_read_length from Fastq.pm v}
      . $MIP::Fastq::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $READ_LENGTH => 151;

## Given a compressed fastq file
my $directory = catdir( $Bin, qw{ data 643594-miptest test_data ADM1059A1 fastq } );

my $file = q{1_161011_TestFilev2_ADM1059A1_TCCGGAGA_1.fastq.gz};

# For compressed file
my $read_file_command = q{gzip -d -c};

my $read_length = get_read_length(
    {
        file_path         => catfile( $directory, $file ),
        read_file_command => $read_file_command,
    }
);

## Then read length should be returned
is( $read_length, $READ_LENGTH, q{Return read length} );

done_testing();
