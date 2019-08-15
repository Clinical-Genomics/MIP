#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
use FindBin qw{ $Bin };
use Getopt::Long;
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
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $COMMA            => q{,};
Readonly my $NEWLINE          => qq{\n};
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

### User Options
GetOptions(

    # Display help text
    q{h|help} => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },

    # Display version number
    q{v|version} => sub {
        done_testing();
        say {*STDOUT} $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $VERSION . $NEWLINE;
        exit;
    },
    q{vb|verbose} => $VERBOSE,
  )
  or (
    done_testing(),
    help(
        {
            USAGE     => $USAGE,
            exit_code => 1,
        }
    )
  );

BEGIN {

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Log::MIP_log4perl} => [qw{ initiate_logger }],
        q{MIP::Script::Utils}     => [qw{ help }],
    );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Get::File});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
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

## Create temp logger
my $test_dir      = File::Temp->newdir();
my $test_log_path = catfile( $test_dir, q{test.log} );

## Creates log object
my $log = initiate_logger(
    {
        file_path => $test_log_path,
        log_name  => q{TEST},
    }
);

## Given file, when format is lesser than Casava 1.4
my $directory = catfile( $Bin, qw{ data 643594-miptest test_data ADM1059A1 fastq} );
my $file_format_illumina_v1_4 = q{1_161011_TestFilev2_ADM1059A1_TCCGGAGA_1.fastq.gz};
my $read_file_command         = q{zcat};

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

######################
####SubRoutines#######
######################

sub build_usage {

## Function  : Build the USAGE instructions
## Returns   :
## Arguments : $program_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $program_name;

    my $tmpl = {
        program_name => {
            default     => basename($PROGRAM_NAME),
            store       => \$program_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help     Display this help message
    -v/--version  Display version
END_USAGE
}
