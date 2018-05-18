#!/usr/bin/env perl

use 5.022;
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
    my @modules = (q{MIP::Recipes::Install::Conda});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Set::Parameter qw{ set_conda_env_names_and_paths };

diag(   q{Test set_conda_env_names_and_paths from MIP::Set::Parameter.pm v}
      . $MIP::Set::Parameter::VERSION
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

## Given a test hash with a defined name for MIP main env (emip) and two sub
## environments whereof one is named.
my %parameter = (
    conda_dir_path   => catdir( $Bin, qw{ data modules miniconda } ),
    environment_name => {
        emip       => q{Test},
        test_env_1 => undef,
        test_env_2 => q{test_2},
    },
    installations => [qw{ emip test_env_1 test_env_2 }],
);

## When the subroutine is being ran
set_conda_env_names_and_paths(
    {
        log            => $log,
        parameter_href => \%parameter,
    }
);
## Then:
## The unnamed sub environment shall be named with the main environment name
## followed by the sub environment
is( $parameter{environment_name}{test_env_1},
    q{Test_test_env_1}, q{Set undefined environment name} );
## The named sub environment shall retain it's name
is( $parameter{environment_name}{test_env_2},
    q{test_2}, q{Leave defined environment name} );
## The path to MIP's conda env shall be set
is(
    $parameter{emip}{conda_prefix_path},
    catdir( $parameter{conda_dir_path}, qw{ envs Test } ),
    q{Set main environment specific conda paths}
);
## The path to the MIP's sub environment shall be set
is(
    $parameter{test_env_2}{conda_prefix_path},
    catdir( $parameter{conda_dir_path}, qw{ envs test_2 } ),
    q{Set sub environment specific conda paths}
);

## Given that MIP's main and sub environments are unnamed
$parameter{environment_name}{emip}       = undef;
$parameter{environment_name}{test_env_1} = undef;

## When the subroutine is being ran
trap {
    set_conda_env_names_and_paths(
        {
            log            => $log,
            parameter_href => \%parameter,
        }
      )
};

## Then:
## The name of sub environemnt shall be the environemnt
is( $parameter{environment_name}{test_env_1},
    q{test_env_1},
    q{Set undefined environment name when emip env name is undefined} );
## The conda path for MIP's main environment shall be conda's base environment
is(
    $parameter{emip}{conda_prefix_path},
    $parameter{conda_dir_path},
    q{Set main environment to base }
);
## The Log shall warn
like( $trap->stderr, qr/WARN/xms,
    q{Warn when installing in conda's base environment} );

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
