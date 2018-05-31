#!/usr/bin/env perl

use 5.018;
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
use Modern::Perl qw{ 2014 };
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
        say {*STDOUT} $NEWLINE
          . basename($PROGRAM_NAME)
          . $SPACE
          . $VERSION
          . $NEWLINE;
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
    my @modules = (q{MIP::Parse::File});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Parse::File qw{ parse_fastq_infiles };

diag(   q{Test parse_fastq_infiles from File.pm v}
      . $MIP::Parse::File::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create temp logger
my $test_dir = File::Temp->newdir();
my $test_log_path = catfile( $test_dir, q{test.log} );

## Creates log object
my $log = initiate_logger(
    {
        file_path => $test_log_path,
        log_name  => q{TEST},
    }
);

## Given compressed file, when proper data
my %active_parameter = ( sample_ids => [qw{ ADM1059A1 }] );
my %file_info;
my %indir_path = ( ADM1059A1 =>
      catdir( $Bin, qw{ data 643594-miptest test_data ADM1059A1 fastq} ), );
my %infile =
  ( ADM1059A1 => [qw{ 1_161011_TestFilev2_ADM1059A1_TCCGGAGA_1.fastq.gz }], );
my %infile_both_strands_prefix;
my %infile_lane_prefix;
my %lane;
my %sample_info;

my $is_file_uncompressed = parse_fastq_infiles(
    {
        active_parameter_href           => \%active_parameter,
        file_info_href                  => \%file_info,
        indir_path_href                 => \%indir_path,
        infile_both_strands_prefix_href => \%infile_both_strands_prefix,
        infile_href                     => \%infile,
        infile_lane_prefix_href         => \%infile_lane_prefix,
        lane_href                       => \%lane,
        log                             => $log,
        sample_info_href                => \%sample_info,
    }
);

## Then return undef
is( $is_file_uncompressed, undef, q{No files uncompressed} );

## Given uncompressed file
push @{ $active_parameter{sample_ids} }, q{ADM1059A2};
push @{ $infile{ADM1059A2} },
  qw{ 1_161011_TestFilev2_ADM1059A2_CGCTCATT_1.fastq ADM1059A2.fastq.gz };
$indir_path{ADM1059A2} =
  catdir( $Bin, qw{ data 643594-miptest test_data bad_input } );

$is_file_uncompressed = parse_fastq_infiles(
    {
        active_parameter_href           => \%active_parameter,
        file_info_href                  => \%file_info,
        indir_path_href                 => \%indir_path,
        infile_both_strands_prefix_href => \%infile_both_strands_prefix,
        infile_href                     => \%infile,
        infile_lane_prefix_href         => \%infile_lane_prefix,
        lane_href                       => \%lane,
        log                             => $log,
        sample_info_href                => \%sample_info,
    }
);

## Then return true
ok( $is_file_uncompressed,
    q{Files uncompressed and got run info from headers} );

## Given file, when no sample_id in file name
push @{ $active_parameter{sample_ids} }, q{ADM1059A3};
push @{ $infile{ADM1059A3} },            qw{ 643594-miptest_pedigree.yaml };
$indir_path{ADM1059A3} = catdir( $Bin, qw{ data 643594-miptest test_data } );

trap {
    parse_fastq_infiles(
        {
            active_parameter_href           => \%active_parameter,
            file_info_href                  => \%file_info,
            indir_path_href                 => \%indir_path,
            infile_both_strands_prefix_href => \%infile_both_strands_prefix,
            infile_href                     => \%infile,
            infile_lane_prefix_href         => \%infile_lane_prefix,
            lane_href                       => \%lane,
            log                             => $log,
            sample_info_href                => \%sample_info,
        }
      )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if sample_id in file name cannot be found} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if sample_id in file name cannot be found} );

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
