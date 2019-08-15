#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir };
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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $COMMA   => q{,};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };
Readonly my $DATE    => q{161011};

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
    my %perl_module = ( q{MIP::Script::Utils} => [qw{ help }], );

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

use MIP::Parse::File qw{ parse_fastq_infiles_format };

diag(   q{Test parse_fastq_infiles_format from File.pm v}
      . $MIP::Parse::File::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given compressed file, when following MIP filename convention
my $file_name = q{1_161011_TestFilev2_ADM1059A1_TCCGGAGA_1.fastq.gz};

## Parse infile according to filename convention
my %infile_info = parse_fastq_infiles_format( { file_name => $file_name, } );

my %expected_features = (
    lane             => 1,
    date             => $DATE,
    flowcell         => q{TestFilev2},
    infile_sample_id => q{ADM1059A1},
    index            => q{TCCGGAGA},
    direction        => 1,
);

## Then return all features from filename
is_deeply( \%infile_info, \%expected_features,
    q{Compressed file follows file convention} );

## Given compressed file, when followinf MIP filename convention
$file_name = q{1_161011_TestFilev2_ADM1059A1_TCCGGAGA_1.fastq};

## Parse infile according to filename convention
%infile_info = parse_fastq_infiles_format( { file_name => $file_name, } );

## Then return all features from filename
is_deeply( \%infile_info, \%expected_features,
    q{Uncompressed file follows file convention} );

## Given file, when not following MIP filename convention
$file_name = q{TestFilev2_ADM1059A1_TCCGGAGA_1.fastq};

## Parse infile according to filename convention
%infile_info = parse_fastq_infiles_format( { file_name => $file_name, } );

## Then return no features from filename
is( keys %infile_info, 0, q{Return 0 features for file not following file convention} );

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
