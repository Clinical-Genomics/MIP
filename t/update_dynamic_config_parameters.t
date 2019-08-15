#!/usr/bin/env perl

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catdir catfile };
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
    my @modules = (q{MIP::Update::Parameters});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Update::Parameters qw{ update_dynamic_config_parameters };

diag(   q{Test update_dynamic_config_parameters from Update::Parameters.pm v}
      . $MIP::Update::Parameters::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %active_parameter = (
    cluster_constant_path  => catfile(qw{ root dir_1 dir_2 }),
    analysis_constant_path => q{analysis},
    case_id                => q{case_1},
    pedigree_file =>
      catfile(qw{ cluster_constant_path! case_id! case_id!_pedigree.yaml }),
    sample_info_file => catfile(
        qw{ cluster_constant_path! case_id! analysis_constant_path! case_id!_qc_sample_info.yaml }
    ),
);

my @order_parameters = qw{ pedigree_file sample_info_file };

## Loop through all parameters and update info
PARAMETER:
foreach my $parameter_name (@order_parameters) {

    ## Updates the active parameters to particular user/cluster for dynamic config parameters following specifications. Leaves other entries untouched.
    update_dynamic_config_parameters(
        {
            active_parameter_href => \%active_parameter,
            parameter_name        => $parameter_name,
        }
    );
}
my $updated_pedigree_file = catfile(qw{ root dir_1 dir_2 case_1 case_1_pedigree.yaml });
is( $active_parameter{pedigree_file},
    $updated_pedigree_file, q{Updated pedigree file path} );

my $updated_sample_info_file =
  catfile(qw{ root dir_1 dir_2 case_1 analysis case_1_qc_sample_info.yaml });
is( $active_parameter{sample_info_file},
    $updated_sample_info_file, q{Updated sample_info_file path} );

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
