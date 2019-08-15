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
our $VERSION = '1.0.2';

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
    my @modules = (q{MIP::Set::Pedigree});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Set::Pedigree qw{ set_pedigree_capture_kit_info };

diag(   q{Test set_pedigree_capture_kit_info from Pedigree.pm v}
      . $MIP::Set::Pedigree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %active_parameter;

my %parameter = (
    supported_capture_kit => {
        default => {
            q{agilent_sureselect.v4} =>
              q{genome_reference_source_version_agilent_sureselect_targets_-v4-.bed},
            q{agilent_sureselect.v5} =>
              q{genome_reference_source_version_agilent_sureselect_targets_-v5-.bed},
        },
    },
);

my %pedigree = (
    case    => q{case_1},
    samples => [
        {
            analysis_type => q{wes},
            capture_kit   => q{agilent_sureselect.v5},
            father        => 0,
            mother        => 0,
            phenotype     => q{affected},
            sample_id     => q{sample_1},
            sex           => q{female},
        },
        {
            analysis_type => q{wgs},
            capture_kit   => q{agilent_sureselect.v4},
            father        => 0,
            mother        => 0,
            phenotype     => q{unaffected},
            sample_id     => q{sample_2},
            sex           => q{male},
        },
        {
            analysis_type => q{wts},
            capture_kit   => q{agilent_sureselect.v5},
            father        => 0,
            mother        => 0,
            phenotype     => q{unknown},
            sample_id     => q{sample_3},
            sex           => q{other},
        },
        {
            analysis_type => q{wgs},
            father        => q{sample_1},
            mother        => q{sample_2},
            phenotype     => q{unknown},
            sample_id     => q{sample_4},
            sex           => q{unknown},
        },
    ],
);

my %sample_info = (
    sample => {
        sample_1 => {
            analysis_type     => q{wes},
            expected_coverage => 30,
            capture_kit       => q{agilent_sureselect.v5},
        },
        sample_2 => {
            analysis_type     => q{wes},
            expected_coverage => 30,
            capture_kit       => q{agilent_sureselect.v4},
        },
        sample_3 => {
            analysis_type     => q{wes},
            expected_coverage => 30,
            capture_kit       => q{agilent_sureselect.v5},
        },
    },
);

my %user_supply_switch = ( exome_target_bed => 0, );
set_pedigree_capture_kit_info(
    {
        active_parameter_href   => \%active_parameter,
        parameter_href          => \%parameter,
        pedigree_href           => \%pedigree,
        sample_info_href        => \%sample_info,
        user_supply_switch_href => \%user_supply_switch,
    }
);

my $capture_kit_string = $active_parameter{exome_target_bed}
  {q{genome_reference_source_version_agilent_sureselect_targets_-v5-.bed}};
is( $capture_kit_string, q{sample_1,sample_3}, q{Set sample_ids for capture kit 1} );

$capture_kit_string = $active_parameter{exome_target_bed}
  {q{genome_reference_source_version_agilent_sureselect_targets_-v4-.bed}};
is( $capture_kit_string, q{sample_2}, q{Set sample_ids for capture kit 2} );

## Test unknown capture kit
%active_parameter = ();
%sample_info      = (
    sample => {
        sample_1 => {
            analysis_type     => q{wes},
            expected_coverage => 30,
            capture_kit       => q{unknown_capture_kit},
        },
    },
);

set_pedigree_capture_kit_info(
    {
        active_parameter_href   => \%active_parameter,
        parameter_href          => \%parameter,
        pedigree_href           => \%pedigree,
        sample_info_href        => \%sample_info,
        user_supply_switch_href => \%user_supply_switch,
    }
);
$capture_kit_string = $active_parameter{exome_target_bed}{q{unknown_capture_kit}};
is( $capture_kit_string, q{sample_1}, q(Unknown capture kit) );

## Test no capture kit
%active_parameter = ();
%sample_info      = ();

set_pedigree_capture_kit_info(
    {
        active_parameter_href   => \%active_parameter,
        parameter_href          => \%parameter,
        pedigree_href           => \%pedigree,
        sample_info_href        => \%sample_info,
        user_supply_switch_href => \%user_supply_switch,
    }
);

is( $active_parameter{exome_target_bed},
    undef, q(No capture kit from cmd, config or pedigree) );

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
