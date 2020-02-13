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

use MIP::Check::Parameter qw{ check_sample_id_in_hash_parameter };

diag(   q{Test check_sample_id_in_hash_parameter from Parameter.pm v}
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

my %active_parameter = ( sample_ids => [qw{ sample-1 sample-2 }], );
my %parameter;

## Given no hash parameter
my $is_skipped = check_sample_id_in_hash_parameter(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        parameter_names_ref   => [qw{ analysis_type }],
        parameter_href        => \%parameter,
        sample_ids_ref        => \@{ $active_parameter{sample_ids} },
    }
);

## Then we have nothing to test
ok( $is_skipped, q{Skipped parameter} );

## Given a hash parameter, when undef sample_ids for parameter
%active_parameter = (
    analysis_type => undef,
    sample_ids    => [qw{ sample-1 sample-2 }],
);

$is_skipped = check_sample_id_in_hash_parameter(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        parameter_names_ref   => [qw{ analysis_type }],
        parameter_href        => \%parameter,
        sample_ids_ref        => \@{ $active_parameter{sample_ids} },
    }
);
## Then we still have nothing to test
ok( $is_skipped, q{Skipped undef sample id for parameter} );

## Given a hash parameter, when parameter is mandatory
%active_parameter = (
    analysis_type => { q{sample-1} => q{wgs}, },
    sample_ids    => [qw{ sample-1 sample-2 }],
);

trap {
    check_sample_id_in_hash_parameter(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_names_ref   => [qw{ analysis_type }],
            parameter_href        => \%parameter,
            sample_ids_ref        => \@{ $active_parameter{sample_ids} },
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if mandatory parameter} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message for mandatory parameter} );

## Given a hash parameter, when parameter is not mandatory
%active_parameter = (
    analysis_type => { q{sample-1} => q{wgs}, },
    sample_ids    => [qw{ sample-1 sample-2 }],
);

%parameter = ( analysis_type => { mandatory => q{no}, }, );

$is_skipped = check_sample_id_in_hash_parameter(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        parameter_names_ref   => [qw{ analysis_type }],
        parameter_href        => \%parameter,
        sample_ids_ref        => \@{ $active_parameter{sample_ids} },
    }
);

## Then we can skip the check
ok( $is_skipped, q{Skipped check for not mandatory parameter} );

## Given a hash parameter, when for ok parameters
%active_parameter = (
    analysis_type => {
        q{sample-1} => q{wgs},
        q{sample-2} => q{wes},
    },
    sample_ids => [qw{ sample-1 sample-2 }],
);

%parameter = ( analysis_type => { mandatory => q{yes}, }, );

my $is_ok = check_sample_id_in_hash_parameter(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        parameter_names_ref   => [qw{ analysis_type }],
        parameter_href        => \%parameter,
        sample_ids_ref        => \@{ $active_parameter{sample_ids} },
    }
);

## Then we can skip the check
ok( $is_ok, q{Passed check for hash parameter} );
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
