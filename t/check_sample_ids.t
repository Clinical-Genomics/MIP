#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
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
    my @modules = (q{MIP::Check::Parameter});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Check::Parameter qw{ check_sample_ids };

diag(   q{Test check_sample_ids from Parameter.pm v}
      . $MIP::Check::Parameter::VERSION
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

## Given case_id, when no sample_ids
my %active_parameter = ( case_id => q{case-1}, );

trap {
    check_sample_ids(
        {
            case_id        => $active_parameter{case_id},
            log            => $log,
            sample_ids_ref => \@{ $active_parameter{sample_ids} },
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if the path cannot be found} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message} );

## Given case_id, when equals sample_ids
%active_parameter = (
    case_id    => q{sample-1},
    sample_ids => [qw{ sample-1 sample-2 }],
);

trap {
    check_sample_ids(
        {
            case_id        => $active_parameter{case_id},
            log            => $log,
            sample_ids_ref => \@{ $active_parameter{sample_ids} },
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if the path cannot be found} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message} );

## Given case_id, when duplicate sample_ids
%active_parameter = (
    case_id    => q{case-1},
    sample_ids => [qw{ sample-1 sample-2 sample-1 }],
);

trap {
    check_sample_ids(
        {
            case_id        => $active_parameter{case_id},
            log            => $log,
            sample_ids_ref => \@{ $active_parameter{sample_ids} },
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if the path cannot be found} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message} );

## Given case_id, when sample_ids contain "-"
%active_parameter = (
    case_id    => q{case-1},
    sample_ids => [qw{ sample_1 sample-2 sample-3 }],
);

trap {
    check_sample_ids(
        {
            case_id        => $active_parameter{case_id},
            log            => $log,
            sample_ids_ref => \@{ $active_parameter{sample_ids} },
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if the path cannot be found} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message} );

## Given case_id, when sample_ids are correct
%active_parameter = (
    case_id    => q{case-1},
    sample_ids => [qw{ sample-1 sample-2 sample-3 }],
);

my $is_ok = check_sample_ids(
    {
        case_id        => $active_parameter{case_id},
        log            => $log,
        sample_ids_ref => \@{ $active_parameter{sample_ids} },
    }
);

## Then  and throw FATAL log message
ok( $is_ok, q{Sample ids are ok} );

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
