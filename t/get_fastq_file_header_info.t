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
use MIP::Constants qw{ $COMMA $NEWLINE $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

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
        q{MIP::Get::File}      => [qw{ get_fastq_file_header_info }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::File qw{ get_fastq_file_header_info };

diag(   q{Test get_fastq_file_header_info from File.pm v}
      . $MIP::Get::File::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $LANE             => 7;
Readonly my $LANE_V_1_8       => 8;
Readonly my $RUN_NUMBER       => 41;
Readonly my $RUN_NUMBER_V_1_8 => 255;
Readonly my $SPACE            => q{ };
Readonly my $TILE             => 1110;
Readonly my $TILE_V_1_8       => 1101;
Readonly my $X_POS            => 21_797;
Readonly my $X_POS_V_1_8      => 6066;
Readonly my $Y_POS            => 11_330;
Readonly my $Y_POS_V_1_8      => 1063;

## Creates log object
my $log = test_log( {} );

## Given file, when format is lesser than Casava 1.4
my $directory = catfile( $Bin, qw{ data 643594-miptest test_data ADM1059A1 fastq} );
my $file_format_illumina_v1_4 = q{1_161011_TestFilev2_ADM1059A1_TCCGGAGA_1.fastq.gz};
my $read_file_command         = q{gzip -d -c};

my %fastq_info_header = get_fastq_file_header_info(
    {
        file_path         => catfile( $directory, $file_format_illumina_v1_4 ),
        log               => $log,
        read_file_command => $read_file_command,
    }
);
my %expected_header_element_v1_4 = (
    instrument_id => q{@ST-E00266},
    run_number    => $RUN_NUMBER,
    flowcell      => q{H2V7FCCXX},
    lane          => $LANE,
    tile          => $TILE,
    x_pos         => $X_POS,
    y_pos         => $Y_POS,
    direction     => 1,
);

## Then return expected_header_elements for < v1.4
is_deeply(
    \%fastq_info_header,
    \%expected_header_element_v1_4,
    q{Found Illumina header format < v1.4}
);

## Given file, when format is greather than Casava 1.8
my $file_format_illumina_v1_8 = q{8_161011_HHJJCCCXY_ADM1059A1_NAATGCGC_1.fastq.gz};

%fastq_info_header = get_fastq_file_header_info(
    {
        file_path         => catfile( $directory, $file_format_illumina_v1_8 ),
        log               => $log,
        read_file_command => $read_file_command,
    }
);
my %expected_header_element_v1_8 = (
    instrument_id => q{@ST-E00269},
    run_number    => $RUN_NUMBER_V_1_8,
    flowcell      => q{HHJJCCCXY},
    lane          => $LANE_V_1_8,
    tile          => $TILE_V_1_8,
    x_pos         => $X_POS_V_1_8,
    y_pos         => $Y_POS_V_1_8,
    direction     => 1,
    filtered      => q{N},
    control_bit   => 0,
    index         => q{NAATGCGC},
);

## Then return expected_header_elements for > v1.8
is_deeply(
    \%fastq_info_header,
    \%expected_header_element_v1_8,
    q{Found Illumina header format > v1.8}
);

## Given file, when format is bad
$directory = catfile( $Bin, qw{ data 643594-miptest } );
my $file_bad_format = q{643594-miptest_pedigree.yaml};

trap {
    get_fastq_file_header_info(
        {
            file_path         => catfile( $directory, $file_bad_format ),
            log               => $log,
            read_file_command => q{cat},
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if format cannot be found} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if format cannot be found } );

done_testing();
