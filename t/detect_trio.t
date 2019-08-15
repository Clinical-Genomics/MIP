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
use Params::Check qw{ allow check last_error };
use open qw{ :encoding(UTF-8) :std };
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
    my @modules = (q{MIP::File::Format::Pedigree});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::File::Format::Pedigree qw{ detect_trio };

diag(   q{Test detect_trio from PEDIGREE.pm v}
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

## Given only single sample
my %active_parameter = ( sample_ids => [qw{ sample_1 }], );
my %sample_info;

my $is_trio = detect_trio(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        sample_info_href      => \%sample_info,
    }
);

## Then do not detect trio
is( $is_trio, undef, q{Single sample - did not detect trio} );

## Given sample info when 3 children are present and no parents
%active_parameter = ( sample_ids => [qw{ child_1 child_2 child_3 }], );
%sample_info      = (
    sample => {
        child_1 => {
            father => 0,
            mother => 0,
        },
        child_2 => {
            father => 0,
            mother => 0,
        },
        child_3 => {
            father => 0,
            mother => 0,
        },
    },
);

$is_trio = detect_trio(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        sample_info_href      => \%sample_info,
    }
);

## Then do not detect trio
is( $is_trio, undef, q{Three children - did not detect trio} );

## Given more samples than a trio, when one child has parents, but not in analysis
%sample_info = (
    sample => {
        child_1 => {
            father => 1,
            mother => 2,
        },
        child_2 => {
            father => 0,
            mother => 0,
        },
        child_3 => {
            father => 0,
            mother => 0,
        },
    },
);

$is_trio = detect_trio(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        sample_info_href      => \%sample_info,
    }
);

## Then do not detect trio
is( $is_trio, undef,
    q{Three children where one is trio but not in analysis - did not detect trio} );

## Given more samples than a trio, when one child has parents in analysis
%active_parameter = ( sample_ids => [qw{ child_1 child_2 child_3 father_1 mother_1 }], );

%sample_info = (
    sample => {
        child_1 => {
            father => q{father_1},
            mother => q{mother_2},
        },
        child_2 => {
            father => 0,
            mother => 0,
        },
        child_3 => {
            father => 0,
            mother => 0,
        },
    },
);

$is_trio = detect_trio(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        sample_info_href      => \%sample_info,
    }
);

## Then do not detect trio due to to many samples
is( $is_trio, undef,
    q{Three children where one is trio in analysis - did not detect trio} );

## Given a trio, when correct number of samples in analysis
%active_parameter = ( sample_ids => [qw{ child_1 father_1 mother_1 }], );

%sample_info = (
    sample => {
        child_1 => {
            father => q{father_1},
            mother => q{mother_1},
        },
        father_1 => {
            father => 0,
            mother => 0,
        },
        mother_1 => {
            father => 0,
            mother => 0,
        },
    },
);

$is_trio = detect_trio(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        sample_info_href      => \%sample_info,
    }
);

## Then do detect trio
is( $is_trio, 1, q{Trio - did detect trio} );
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
