#!/usr/bin/env perl

use 5.026;
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
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;

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
    my @modules = (q{MIP::Get::Parameter});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Get::Parameter qw{ get_capture_kit };

diag(   q{Test get_capture_kit from Parameter.pm v}
      . $MIP::Get::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Test defined switch but no user info - default ie. "latest" capture kit
my $capture_kit = q{latest};

my %parameter = (
    supported_capture_kit => {
        q{agilent_sureselect.v5} =>
          q{genome_reference_source_version_agilent_sureselect_targets_-v5-.bed},
        latest =>
          q{genome_reference_source_version_agilent_sureselect_targets_cre_-v1-.bed},
    },
);

my %user_supply_switch = ( exome_target_bed => 0, );

my $exome_target_bed_file = get_capture_kit(
    {
        capture_kit                    => $capture_kit,
        supported_capture_kit_href     => $parameter{supported_capture_kit},
        user_supplied_parameter_switch => $user_supply_switch{exome_target_bed},
    }
);
is(
    $exome_target_bed_file,
    q{genome_reference_source_version_agilent_sureselect_targets_cre_-v1-.bed},
    q{Got latest default capture kit}
);

## Test defined switch but no user info - no default capture kit
$capture_kit           = q{agilent_sureselect.v4};
$exome_target_bed_file = get_capture_kit(
    {
        capture_kit                    => $capture_kit,
        supported_capture_kit_href     => $parameter{supported_capture_kit},
        user_supplied_parameter_switch => $user_supply_switch{exome_target_bed},
    }
);
is( $exome_target_bed_file, q{agilent_sureselect.v4}, q{Got supplied capture kit} );

## Test undefined switch but no user info - default capture kit
$capture_kit        = q{agilent_sureselect.v5};
%user_supply_switch = ();

$exome_target_bed_file = get_capture_kit(
    {
        capture_kit                    => $capture_kit,
        supported_capture_kit_href     => $parameter{supported_capture_kit},
        user_supplied_parameter_switch => $user_supply_switch{exome_target_bed},
    }
);
is(
    $exome_target_bed_file,
    q{genome_reference_source_version_agilent_sureselect_targets_-v5-.bed},
    q{Got capture kit with undef user switch}
);

## ## Test undefined switch but no user info - no default capture kit
$capture_kit           = q{agilent_sureselect.v4};
$exome_target_bed_file = get_capture_kit(
    {
        capture_kit                    => $capture_kit,
        supported_capture_kit_href     => $parameter{supported_capture_kit},
        user_supplied_parameter_switch => $user_supply_switch{exome_target_bed},
    }
);
is( $exome_target_bed_file, q{agilent_sureselect.v4},
    q{Got supplied capture kit with undef user switch} );

## Test defined switch and user info - Get no cpature kit
%user_supply_switch = ( exome_target_bed => 1, );

$exome_target_bed_file = get_capture_kit(
    {
        capture_kit                    => $capture_kit,
        supported_capture_kit_href     => $parameter{supported_capture_kit},
        user_supplied_parameter_switch => $user_supply_switch{exome_target_bed},
    }
);
is( $exome_target_bed_file, undef, q{Did not get any capture kit} );

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
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
