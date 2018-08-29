#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use Getopt::Long;
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie;
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.1';

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
    my %perl_module = ( q{MIP::Script::Utils} => [qw{ help }], );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Get::Parameter});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Get::Parameter qw{ get_module_parameters };

diag(   q{Test get_module_parameters from Parameter.pm v}
      . $MIP::Get::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $program_name = q{chanjo_sexcheck};

my %active_parameter = (
    module_time                      => { chanjo_sexcheck => 1, },
    module_core_number               => { chanjo_sexcheck => 1, },
    source_main_environment_commands => [ qw{ source activate mip }, ],
);

my ( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
    {
        active_parameter_href => \%active_parameter,
        program_name          => $program_name,
    }
);

is( $core_number, 1, q{Got module core_number} );

is( $time, 1, q{Got module time} );

is_deeply( \@source_environment_cmds, [qw{source activate mip}],
    q{Got main source environment command} );

## Test module specific source command
%active_parameter = (
    module_time        => { chanjo_sexcheck => 1, },
    module_core_number => { chanjo_sexcheck => 1, },
    module_source_environment_command =>
      { chanjo_sexcheck => [ qw{ source activate python_v3.6_tools }, ], },
);

( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
    {
        active_parameter_href => \%active_parameter,
        program_name          => $program_name,
    }
);

is_deeply(
    \@source_environment_cmds,
    [qw{ source activate python_v3.6_tools }],
    q{Got module source environment commmand}
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
            strict_type => 1,
            store       => \$program_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
