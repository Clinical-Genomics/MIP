#!/usr/bin/env perl

use Modern::Perl qw{2014};
use warnings qw{FATAL utf8};
use autodie;
use 5.018;    #Require at least perl 5.18
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{-no_match_vars};
use Params::Check qw{check allow last_error};

use FindBin qw{$Bin};    #Find directory of script
use File::Basename qw{dirname basename};
use File::Spec::Functions qw{catdir};
use Getopt::Long;
use Test::More;
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use Script::Utils qw{help};

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $SPACE   => q{ };
Readonly my $NEWLINE => qq{\n};

###User Options
GetOptions(
    q{h|help} => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },    #Display help text
    q{v|version} => sub {
        done_testing();
        say {*STDOUT} $NEWLINE, basename($PROGRAM_NAME),
          $SPACE, $VERSION, $NEWLINE;
        exit;
    },    #Display version number
    q{vb|verbose} => $VERBOSE,
  )
  or (
    done_testing(),
    Script::Utils::help(
        {
            USAGE     => $USAGE,
            exit_code => 1,
        }
    )
  );

BEGIN {

### Check all internal dependency modules and imports
##Modules with import
    my %perl_module;

    $perl_module{'Script::Utils'} = [qw{help}];

  PERL_MODULES:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load } . $module;
    }

## Modules
    my @modules = (q{MIP::Update::Programs});

  MODULES:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load } . $module;
    }
}

use MIP::Update::Programs qw{update_program_mode_with_dry_run_all};

diag(
    q{Test update_program_mode_with_dry_run_all }
      . $MIP::Update::Programs::VERSION
      . q{, Perl}
      . $PERL_VERSION,
    $EXECUTABLE_NAME
);

### No update of program parameters

my $simulation_mode = 0;

my @programs = qw{pfastqc pbwa_mem ppeddy};

my %active_parameter = (
    pfastqc  => 0,
    pbwa_mem => 1,
    ppeddy   => 2,
);

update_program_mode_with_dry_run_all(
    {
        active_parameter_href => \%active_parameter,
        programs_ref          => \@programs,
        dry_run_all           => $simulation_mode,
    }
);

is( $active_parameter{pfastqc}, 0, q{No update pfastqc} );

is( $active_parameter{pbwa_mem}, 1, q{No update pbwa_mem} );

is( $active_parameter{ppeddy}, 2, q{No update ppeddy} );

### Update of program parameters

# Set simulation mode
$simulation_mode = 1;

update_program_mode_with_dry_run_all(
    {
        active_parameter_href => \%active_parameter,
        programs_ref          => \@programs,
        dry_run_all           => $simulation_mode,
    }
);

is( $active_parameter{pfastqc}, 0, q{No update pfastqc} );

is( $active_parameter{pbwa_mem}, 2, q{Update pbwa_mem to simulation mode} );

is( $active_parameter{ppeddy}, 2, q{No update ppeddy} );

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

    check( $tmpl, $arg_href, 1 ) or croak qw(Could not parse arguments!);

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
