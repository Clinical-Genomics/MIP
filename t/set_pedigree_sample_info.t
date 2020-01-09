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
    my @modules = (q{MIP::Set::Pedigree});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Set::Pedigree qw{ set_pedigree_sample_info };

diag(   q{Test set_pedigree_sample_info from Pedigree.pm v}
      . $MIP::Set::Pedigree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $TOTAL_SAMPLE_IDS => 4;
Readonly my $SAMPLE_3_INDEX   => 3;

my %active_parameter;

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
        },
    },
);

my @user_input_sample_ids;
my %user_supply_switch = ( sample_ids => 0 );

my @pedigree_sample_ids = set_pedigree_sample_info(
    {
        active_parameter_href     => \%active_parameter,
        pedigree_href             => \%pedigree,
        sample_info_href          => \%sample_info,
        user_supply_switch_href   => \%user_supply_switch,
        user_input_sample_ids_ref => \@user_input_sample_ids,
    }
);

## Return of sample_ids
is( scalar @pedigree_sample_ids, $TOTAL_SAMPLE_IDS, q{Found all sample_ids} );

## Transfer to active parameters
is( scalar @{ $active_parameter{sample_ids} },
    $TOTAL_SAMPLE_IDS, q{Set all sample_ids to active parameter} );

## Sample level addition to sample info
foreach my $pedigree_sample_href ( @{ $pedigree{samples} } ) {

    # Alias
    my $sample_id = $pedigree_sample_href->{sample_id};

    foreach my $key ( keys %{$pedigree_sample_href} ) {

        is( $sample_info{sample}{$sample_id}{$key}, $pedigree_sample_href->{$key},
                q{Set key: }
              . $key
              . q{ to '}
              . $pedigree_sample_href->{$key}
              . q{' for '}
              . $sample_id
              . q{'} );
    }
}

## User input
@user_input_sample_ids = qw{ sample_1 sample_2 };
%user_supply_switch    = ( sample_ids => 1 );
%active_parameter      = ();
%sample_info           = ();

@pedigree_sample_ids = set_pedigree_sample_info(
    {
        active_parameter_href     => \%active_parameter,
        pedigree_href             => \%pedigree,
        sample_info_href          => \%sample_info,
        user_supply_switch_href   => \%user_supply_switch,
        user_input_sample_ids_ref => \@user_input_sample_ids,
    }
);

## Sample level addition to sample info for user supplied sample ids
foreach my $key ( keys %{ $sample_info{sample}{sample_1} } ) {

    is( $sample_info{sample}{sample_1}{$key}, $pedigree{samples}[0]{$key},
            q{Set key: }
          . $key
          . q{ to '}
          . $pedigree{samples}[0]{$key}
          . q{' for ' sample_1 '} );
}

## Sample level addition to sample info for not supplied sample ids
is( $sample_info{sample}{sample_3}{analysis_type}, undef,
        q{Did not set key: analysis_type to '}
      . $pedigree{samples}[$SAMPLE_3_INDEX]{analysis_type}
      . q{' for ' sample_3 '} );

is( $active_parameter{sample_ids}, undef, q{Did not set sample_ids to active parameter} );

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
