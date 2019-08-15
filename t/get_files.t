#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
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
    my @modules = (q{MIP::Get::File});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Get::File qw{ get_files };

diag(   q{Test get_files from File.pm v}
      . $MIP::Get::File::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given an infile directory, when applying rules
my $infile_directory =
  catfile( $Bin, qw{ data 643594-miptest test_data ADM1059A1 fastq } );

my @infiles = get_files(
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
@infiles = get_files(
    {
        file_directory => $infile_directory,
        rule_name      => q{*.fastq*},
    }
);

push @expected_files, q{test.fastq.gz};

## Then find all files recursively
is_deeply( \@infiles, \@expected_files, q{Found all files recursively} );

## Given an infile directory, when applying no rules
@infiles = get_files( { file_directory => $infile_directory, } );

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
