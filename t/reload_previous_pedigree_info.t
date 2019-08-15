#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
use FindBin qw{ $Bin };
use Getopt::Long;
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
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $COMMA             => q{,};
Readonly my $NEWLINE           => qq{\n};
Readonly my $SPACE             => q{ };
Readonly my $EXPECTED_COVERAGE => 150;

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
    my @modules = (q{MIP::File::Format::Pedigree});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::File::Format::Pedigree qw{ reload_previous_pedigree_info };

diag(   q{Test reload_previous_pedigree_info from Pedigree.pm v}
      . $MIP::File::Format::Pedigree::VERSION
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

my %sample_info = (
    sample => {
        sample_1 => {
            analysis_type => q{wes},
            phenotype     => q{affected},
            sex           => q{female},
        },
        sample_3 => {
            analysis_type => q{wts},
            phenotype     => q{unknown},
            sex           => q{other},
        },
        sample_4 => {
            analysis_type => q{wgs},
            phenotype     => q{unknown},
            sex           => q{unknown},
        },
    },
);

my %active_parameter =
  ( sample_info_file => catfile( $Bin, qw{data test_data file_does_not_exists} ), );

reload_previous_pedigree_info(
    {
        log                   => $log,
        sample_info_href      => \%sample_info,
        sample_info_file_path => $active_parameter{sample_info_file},
    }
);

is( $sample_info{sample}{sample_1}{analysis_type}, q{wes}, q{No qc_sample_info to load} );

## Test for previous info to be reloaded
%active_parameter =
  ( sample_info_file => catfile( $Bin, qw{data test_data qc_sample_info.yaml} ), );

reload_previous_pedigree_info(
    {
        log                   => $log,
        sample_info_href      => \%sample_info,
        sample_info_file_path => $active_parameter{sample_info_file},
    }
);

is( $sample_info{sample}{sample_1}{analysis_type}, q{wes}, q{Updated parameter} );

is( $sample_info{sample}{sample_1}{capture_kit},
    q{agilent_sureselect.v5},
    q{Reloaded sample parameter not present in current analysis} );

my %reloaded_sample_2 = (
    analysis_type     => q{wgs},
    capture_kit       => q{agilent_sureselect.v5},
    expected_coverage => $EXPECTED_COVERAGE,
    father            => 0,
    mother            => 0,
    phenotype         => q{unaffected},
    sample_id         => q{sample_2},
    sample_name       => q{sample_2},
    sex               => q{male},
);

while ( my ( $key, $value ) = each %reloaded_sample_2 ) {

    is( $sample_info{sample}{sample_2}{$key},
        $value, q{Reloaded previous parameter: } . $key );
}

is( $sample_info{sample}{sample_4}{sex}, q{unknown}, q{Keept new parameter} );

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
