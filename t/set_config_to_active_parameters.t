#!/usr/bin/env perl

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use Getopt::Long;
use Params::Check qw{ check allow last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };
use 5.026;

## CPANM
use autodie;
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $SPACE   => q{ };
Readonly my $NEWLINE => qq{\n};
Readonly my $COMMA   => q{,};

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
    my @modules = (q{MIP::Set::Parameter});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Set::Parameter qw{ set_config_to_active_parameters };

diag(   q{Test set_config_to_active_parameters from Parameter.pm v}
      . $MIP::Set::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %config = (
    hash   => { scalar => q{config_scalar_1}, },
    array  => [qw{ config_array_1 config_array_2 }],
    scalar => q{config_scalar_2},
);

my %active_parameter;

set_config_to_active_parameters(
    {
        config_parameter_href => \%config,
        active_parameter_href => \%active_parameter,
    }
);
## Tests

is( $active_parameter{hash}{scalar}, q{config_scalar_1}, q{Set hash from config} );

is( $active_parameter{array}[0], q{config_array_1}, q{Set array from config} );

is( $active_parameter{scalar}, q{config_scalar_2}, q{Set scalar from config} );

%active_parameter = (
    hash   => { scalar => q{active_scalar_1} },
    array  => [qw{ active_array_1 active_array_2 }],
    scalar => q{active_scalar_2},
);

set_config_to_active_parameters(
    {
        config_parameter_href => \%config,
        active_parameter_href => \%active_parameter,
    }
);

is( $active_parameter{hash}{scalar},
    q{active_scalar_1}, q{Did not set hash from config} );

is( $active_parameter{array}[0], q{active_array_1}, q{Did not set array from config} );

is( $active_parameter{scalar}, q{active_scalar_2}, q{Did not set scalar from config} );

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
