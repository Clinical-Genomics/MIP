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
Readonly my $COMMA   => q{,};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

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
    my @modules = (q{MIP::Check::File});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
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

## Given interleaved file, when casava version less than 1.4
my $sample_id = q{ADM1059A1};
my $directory = catdir( $Bin, qw{ data 643594-miptest test_data ADM1059A1 fastq } );
my $interleaved_file  = q{2_161011_TestFilev2-Interleaved_ADM1059A1_TCCGGAGA_1.fastq.gz};
my $read_file_command = q{zcat};

my @returns = trap {
    check_interleaved(
        {
            file_path         => catfile( $directory, $interleaved_file ),
            log               => $log,
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
            log               => $log,
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
            log               => $log,
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
            log               => $log,
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
            log               => $log,
            read_file_command => $read_file_command,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if interleaved file format cannot be found} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if interleaved file format cannot be found} );

## Given malformed file
$directory = catdir( $Bin, qw{ data 643594-miptest } );
my $file = q{643594-miptest_pedigree.yaml};
$read_file_command = q{cat};

trap {
    check_interleaved(
        {
            file_path         => catfile( $directory, $file ),
            log               => $log,
            read_file_command => $read_file_command,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if malformed file} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message if malformed file} );

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
