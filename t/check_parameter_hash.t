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
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
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
    my @modules = (q{MIP::Check::Parameter});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Check::Parameter qw{ check_parameter_hash };
use MIP::File::Format::Yaml qw{ load_yaml };

diag(   q{Test check_parameter_hash from Parameter.pm v}
      . $MIP::Check::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $definitions_file =
  catfile( $Bin, qw{ data test_data define_parameters.yaml } );

## Loads a YAML file into an arbitrary hash and returns it.
my %parameter = load_yaml( { yaml_file => $definitions_file, } );

## Load mandatory keys and values for parameters
my %mandatory_key = load_yaml(
    {
        yaml_file => catfile(
            dirname($Bin), qw{ definitions mandatory_parameter_keys.yaml }
        ),
    }
);

## Load non mandatory keys and values for parameters
my %non_mandatory_key = load_yaml(
    {
        yaml_file => catfile(
            dirname($Bin), qw{ definitions non_mandatory_parameter_keys.yaml }
        ),

    }
);

my $error_msg = check_parameter_hash(
    {
        parameter_href         => \%parameter,
        mandatory_key_href     => \%mandatory_key,
        non_mandatory_key_href => \%non_mandatory_key,
        file_path              => $definitions_file,
    }
);

is( $error_msg, undef, q{No errors detected} );

### Introduce errors
## Test if mandatory key is missing
$parameter{WRONG}{values} = [qw{ WRONG_ARRAY }];

$error_msg = check_parameter_hash(
    {
        parameter_href         => \%parameter,
        mandatory_key_href     => \%mandatory_key,
        non_mandatory_key_href => \%non_mandatory_key,
        file_path              => $definitions_file,
    }
);

isnt( $error_msg, undef, q{Mandatory key is missing and detected} );

## Reset
delete $parameter{WRONG};

## Add wrong mandatory data_type type
$parameter{mip}{associated_program} = q{SHOULD_BE_ARRAY};

$error_msg = check_parameter_hash(
    {
        parameter_href         => \%parameter,
        mandatory_key_href     => \%mandatory_key,
        non_mandatory_key_href => \%non_mandatory_key,
        file_path              => $definitions_file,
    }
);

isnt( $error_msg, undef, q{Wrong mandatory data_type type detected} );

## Reset
$parameter{mip}{associated_program} = [qw{ ARRAY }];

## Add wrong mandatory data_type value
$parameter{mip}{data_type} = q{NOT_EXACTLY_SCALAR};

$error_msg = check_parameter_hash(
    {
        parameter_href         => \%parameter,
        mandatory_key_href     => \%mandatory_key,
        non_mandatory_key_href => \%non_mandatory_key,
        file_path              => $definitions_file,
    }
);

isnt( $error_msg, undef, q{Wrong mandatory data_type value detected} );

## Reset
$parameter{mip}{data_type} = q{SCALAR};

## Add wrong non mandatory type
$parameter{email_types}{mandatory} = [qw{ NOT_EXACTLY_A_SCALAR }];

$error_msg = check_parameter_hash(
    {
        parameter_href         => \%parameter,
        mandatory_key_href     => \%mandatory_key,
        non_mandatory_key_href => \%non_mandatory_key,
        file_path              => $definitions_file,
    }
);

isnt( $error_msg, undef, q{Wrong non mandatory data_type detected} );

## Reset
$parameter{email_types}{mandatory} = q{SCALAR};

## Add wrong non mandatory value
$parameter{email_types}{mandatory} = q{NEITHER_yes_OR_no};

$error_msg = check_parameter_hash(
    {
        parameter_href         => \%parameter,
        mandatory_key_href     => \%mandatory_key,
        non_mandatory_key_href => \%non_mandatory_key,
        file_path              => $definitions_file,
    }
);

isnt( $error_msg, undef, q{Wrong non mandatory value detected} );

## Reset
$parameter{email_types}{mandatory} = q{SCALAR};

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
