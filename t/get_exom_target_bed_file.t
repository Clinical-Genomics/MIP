#!/usr/bin/env perl

use Modern::Perl qw{ 2014 };
use warnings qw{ FATAL utf8 };
use autodie;
use 5.018;
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use FindBin qw{ $Bin };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catdir catfile };
use Getopt::Long;
use Test::More;
use File::Temp;

## CPANM
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };
use MIP::Log::MIP_log4perl qw{ initiate_logger };

our $USAGE = build_usage( {} );

##Constants
Readonly my $COMMA      => q{,};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

my $VERBOSE = 1;
our $VERSION = q{1.0.0};

###User Options
GetOptions(

    # Display help text
    'h|help' => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },

    # Display version number
    'v|version' => sub {
        done_testing();
        say {*STDOUT} $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $VERSION,
          $NEWLINE;
        exit;
    },
    'vb|verbose' => $VERBOSE,
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
    my %perl_module;

    $perl_module{q{MIP::Script::Utils}}     = [qw{ help }];
    $perl_module{q{MIP::Log::MIP_log4perl}} = [qw{ initiate_logger }];

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load } . $module;
    }

    ## Modules
    my @modules = (q{MIP::Get::File});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load } . $module;
    }
}

use MIP::Get::File qw{ get_exom_target_bed_file };

diag(   q{Test get_exom_target_bed_file from File.pm v}
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
my $log           = initiate_logger(
    {
        file_path => $test_log_path,
        log_name  => q{MIP},
    }
);

## Base arguments
my $sample_1         = q{sample_1};
my $sample_2         = q{sample_2};
my $sample_3         = q{sample_3};
my %exome_target_bed = (
    q{exome_target_file.bed}   => $sample_1 . $COMMA . $sample_2,
    q{exome_target_file_2.bed} => $sample_3,
);

my $exome_target_bed_file = get_exom_target_bed_file(
    {
        exome_target_bed_href => \%exome_target_bed,
        sample_id             => $sample_1,
        log                   => $log,
    }
);

is( q{exome_target_file.bed}, $exome_target_bed_file,
    q{Get exom target bed file for sample_1} );

$exome_target_bed_file = get_exom_target_bed_file(
    {
        exome_target_bed_href => \%exome_target_bed,
        sample_id             => $sample_1,
        log                   => $log,
        file_ending           => q{.interval_list},
    }
);

is( q{exome_target_file.bed.interval_list},
    $exome_target_bed_file,
    q{Get exom target bed.interval_list file for sample_1} );

$exome_target_bed_file = get_exom_target_bed_file(
    {
        exome_target_bed_href => \%exome_target_bed,
        sample_id             => $sample_3,
        log                   => $log,
    }
);

is( q{exome_target_file_2.bed}, $exome_target_bed_file,
    q{Get exom target bed file for sample_3} );

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

##build_usage

##Function : Build the USAGE instructions
##Returns  : ""
##Arguments: $program_name
##         : $program_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $program_name;

    my $tmpl = {
        program_name => {
            default     => basename($PROGRAM_NAME),
            strict_type => 1,
            store       => \$program_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
