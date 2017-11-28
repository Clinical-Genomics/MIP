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
use File::Spec::Functions qw{ catdir };
use Getopt::Long;
use Test::More;

## CPANM
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
    my %perl_module;

    $perl_module{q{MIP::Script::Utils}} = [qw{ help }];

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Update::Programs});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Update::Programs qw{ update_program_mode };

diag(   q{Test update_program_mode from Programs.pm v}
      . $MIP::Update::Programs::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my @programs         = qw{ cnvnator delly_call delly_reformat tiddit };
my %active_parameter = (
    pmanta          => 1,
    pdelly_call     => 1,
    pdelly_reformat => 1,
    pcnvnator       => 1,
    ptiddit         => 1,
);

my @warning_msgs = update_program_mode(
    {
        active_parameter_href   => \%active_parameter,
        programs_ref            => \@programs,
        consensus_analysis_type => q{wgs},
    }
);

is( @warning_msgs, 0, q{No updates to programs mode} );

@warning_msgs = update_program_mode(
    {
        active_parameter_href   => \%active_parameter,
        programs_ref            => \@programs,
        consensus_analysis_type => q{wes},
    }
);
## Alias
my $manta_mode          = $active_parameter{pmanta};
my $delly_call_mode     = $active_parameter{pdelly_call};
my $delly_reformat_mode = $active_parameter{pdelly_reformat};
my $tiddit_mode         = $active_parameter{ptiddit};
my $cnvnator_mode       = $active_parameter{pcnvnator};

## Test program mode updates and warnings
is( $manta_mode,          1, q{Updated programs mode for manta} );
is( $delly_call_mode,     0, q{Updated programs mode for delly_call} );
is( $delly_reformat_mode, 0, q{Updated programs mode for delly_reformat} );
is( $tiddit_mode,         0, q{Updated programs mode for tiddit} );
is( $cnvnator_mode,       0, q{Updated programs mode for cnvnator} );
isnt( @warning_msgs, 0, q{Generated warning message} );

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
