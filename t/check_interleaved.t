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
use Test::Trap qw{ :stderr:output(systemsafe) };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Check::File}    => [qw{ check_interleaved }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::File qw{ check_interleaved };

diag(   q{Test check_interleaved from File.pm v}
      . $MIP::Check::File::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

## Given interleaved file, when casava version less than 1.4
my $directory = catdir( $Bin, qw{ data 643594-miptest test_data ADM1059A1 fastq } );
my $interleaved_file  = q{2_161011_TestFilev2-Interleaved_ADM1059A1_TCCGGAGA_1.fastq.gz};
my $read_file_command = q{gzip -d -c};

my @returns = trap {
    check_interleaved(
        {
            file_path         => catfile( $directory, $interleaved_file ),
            read_file_command => $read_file_command,
        }
    )
};

## Then return true if detected interleaved
ok( $returns[0], q{Detected interleaved casava < 1.4 fastq file} );

## Given file, when casava version < 1.4 without dash in instrument id
$interleaved_file = q{2_161011_TestFilev2_ADM1059A1_TCCGGAGA_1.fastq.gz};

@returns = trap {
    check_interleaved(
        {
            file_path         => catfile( $directory, $interleaved_file ),
            read_file_command => $read_file_command,
        }
    )
};

## Then return true if detected interleaved
ok( $returns[0],
    q{Detected interleaved casava < 1.4 fastq file without dash in instrument id} );

## Given file, when casava version >= 1.8
$interleaved_file = q{8_161011_HHJJCCCXY_ADM1059A1_NAATGCGC_1.fastq.gz};

@returns = trap {
    check_interleaved(
        {
            file_path         => catfile( $directory, $interleaved_file ),
            read_file_command => $read_file_command,
        }
    )
};

## Then return true if detected interleaved
ok( $returns[0], q{Detected interleaved casava >= 1.8 fastq file} );

## Given file, when casava version >= 1.8 without dash in instrument id
$interleaved_file = q{7_161011_HHJJCCCXY_ADM1059A1_NAATGCGC_1.fastq.gz};

@returns = trap {
    check_interleaved(
        {
            file_path         => catfile( $directory, $interleaved_file ),
            read_file_command => $read_file_command,
        }
    )
};

## Then return true if detected interleaved
ok( $returns[0],
    q{Detected interleaved casava >= 1.8 fastq file without dash in instrument id} );

## Given wrong read direction file
$directory        = catdir( $Bin, qw{ data 643594-miptest test_data bad_input } );
$interleaved_file = q{3_161011_TestFilev2-Interleaved_ADM1059A1_TCCGGAGA_1.fastq.gz};

trap {
    check_interleaved(
        {
            file_path         => catfile( $directory, $interleaved_file ),
            read_file_command => $read_file_command,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if interleaved file format cannot be found} );
like(
    $trap->stderr,
    qr/Read \s+ direction \s+ is/xms,
    q{Throw fatal log message if interleaved file format cannot be found}
);

## Given malformed file
$directory = catdir( $Bin, qw{ data 643594-miptest } );
my $file = q{643594-miptest_pedigree.yaml};
$read_file_command = q{cat};

trap {
    check_interleaved(
        {
            file_path         => catfile( $directory, $file ),
            read_file_command => $read_file_command,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if malformed file} );
like(
    $trap->stderr,
    qr/Malformed \s+ fastq \s+ file/xms,
    q{Throw fatal log message if malformed file}
);

done_testing();
