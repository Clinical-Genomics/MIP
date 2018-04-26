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
use Params::Check qw{ allow check last_error };
use open qw{ :encoding(UTF-8) :std };
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
    my @modules = (q{MIP::Check::Path});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Check::Path qw{ check_command_in_path };

diag(   q{Test check_command_in_path from Path.pm v}
      . $MIP::Check::Path::VERSION
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

my %active_parameter = (
    conda_path => catfile( $Bin, qw{ data modules miniconda } ),
    pbwa_mem   => 0,
);
my %parameter;

## Given switched off active parameter, when no parameter defined
my $return = check_command_in_path(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        parameter_href        => \%parameter,
    }
);

## Then return undef
is( $return, undef, q{Skip check when no parameter defined in parameter hash} );

## Given switched off active parameter and defined program parameter, when no defined
## program_name_path
%parameter = (
    pbwa_mem => { type => q{program}, },

);
$return = check_command_in_path(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        parameter_href        => \%parameter,
    }
);

## Then return undef
is( $return, undef,
    q{Skip check when no program_name_path defined in parameter hash} );

## Given switched on active parameter, defined program and program_name_path parameter, when program is in path and executable
%active_parameter = (
    conda_path => catfile( $Bin, qw{ data modules miniconda } ),
    pbwa_mem   => 1,
);

%parameter = (
    pbwa_mem => {
        program_name_path => [qw{ bwa }],
        type              => q{program},
    },
);

trap {
    check_command_in_path(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
        }
      )
};

## Then INFO message should broadcast
like( $trap->stderr, qr/INFO/xms,
    q{Found bin and executable: Throw INFO log message} );

## Given switched on active parameter, defined program and program_name_path parameter, when program is in path and executable
%active_parameter = (
    conda_path => catfile( $Bin, qw{ data modules miniconda } ),
    pbwa_mem   => 1,
);

%parameter = (
    pbwa_mem => {
        program_name_path => [qw{ no_binary }],
        type              => q{program},
    },
);

trap {
    check_command_in_path(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
        }
      )
};

## Then FATAL message should broadcast
like( $trap->stderr, qr/FATAL/xms,
    q{No bin and executable - Throw FATAL log message} );

## Given switched on active parameter, defined program and program_name_path
## parameter, when program is in env path and executable for
## module source environment command
%active_parameter = (
    conda_path => catfile( $Bin, qw{ data modules miniconda } ),
    module_source_environment_command => {
        prankvariant => q{source activate test_env_1},
    },
    prankvariant => 1,
);

%parameter = (
    prankvariant => {
        program_name_path => [qw{ genmod }],
        type              => q{program},
    },
);

trap {
    check_command_in_path(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
        }
      )
};

## Then INFO message should broadcast
like( $trap->stderr, qr/INFO/xms,
    q{Found bin and executable module source env cmd: Throw INFO log message} );

## Given switched on active parameter, defined program and program_name_path
## parameter, when program is in env path and not executable for
## module source environment command
%active_parameter = (
    conda_path => catfile( $Bin, qw{ data modules miniconda } ),
    module_source_environment_command => {
        prankvariant => q{source activate test_env_1},
    },
    prankvariant => 1,
);

%parameter = (
    prankvariant => {
        program_name_path => [qw{ no_binary }],
        type              => q{program},
    },
);

trap {
    check_command_in_path(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
        }
      )
};

## Then FATAL message should broadcast
like( $trap->stderr, qr/FATAL/xms,
    q{No bin and executable in module source env cmd - Throw FATAL log message}
);

## Given switched on active parameter, when program is in env path and
## executable for program source environment command
%active_parameter = (
    conda_path => catfile( $Bin, qw{ data modules miniconda } ),
    program_source_environment_command => {
        genmod => q{source activate test_env_1},
    },
);

trap {
    check_command_in_path(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
        }
      )
};

## Then INFO message should broadcast
like( $trap->stderr, qr/INFO/xms,
q{Found bin and executable in program source env cmd: Throw INFO log message}
);

## Given switched on active parameter, defined program and program_name_path
## parameter, when program is in env path and not executable for
## program source environment command
%active_parameter = (
    conda_path => catfile( $Bin, qw{ data modules miniconda } ),
    program_source_environment_command => {
        not_a_binary => q{source activate test_env_1},
    },
);

trap {
    check_command_in_path(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
        }
      )
};

## Then FATAL message should broadcast
like( $trap->stderr, qr/FATAL/xms,
    q{No bin and executable in program source env cmd - Throw FATAL log message}
);

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
