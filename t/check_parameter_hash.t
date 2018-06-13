#!/usr/bin/env perl

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catfile catdir };
use FindBin qw{ $Bin };
use Getopt::Long;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };
use 5.018;

## CPANM
use Modern::Perl qw{ 2014 };
use autodie;
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = 1.0.0;

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

## Given valid input
my $is_ok = check_parameter_hash(
    {
        parameter_href         => \%parameter,
        mandatory_key_href     => \%mandatory_key,
        non_mandatory_key_href => \%non_mandatory_key,
        file_path              => $definitions_file,
    }
);

## Then return true
ok( $is_ok, q{No errors in parameters detected} );

## Given wrong data type, when SCALAR
$parameter{family_id}{data_type} = [q{wrong_data_type}];

trap {
    check_parameter_hash(
        {
            parameter_href         => \%parameter,
            mandatory_key_href     => \%mandatory_key,
            non_mandatory_key_href => \%non_mandatory_key,
            file_path              => $definitions_file,
        }
      )
};

## Then throw FATAL log message and croak
like( $trap->stderr, qr/ARRAY/xms,
    q{Throw fatal log message for wrong data type - array} );

## Reset parameter
$parameter{family_id}{data_type} = q{SCALAR};

## Given wrong data type, when ARRAY or HASH
$parameter{family_id}{associated_program} = q{not_an_array};

trap {
    check_parameter_hash(
        {
            parameter_href         => \%parameter,
            mandatory_key_href     => \%mandatory_key,
            non_mandatory_key_href => \%non_mandatory_key,
            file_path              => $definitions_file,
        }
      )
};

## Then throw FATAL log message and croak
like( $trap->stderr, qr/SCALAR/xms,
    q{Throw fatal log message for wrong data type - scalar} );

## Reset parameter
$parameter{family_id}{associated_program} = [qw{ mip }];

## Given not allowed value
$parameter{family_id}{data_type} = q{not_valid_value};

trap {
    check_parameter_hash(
        {
            parameter_href         => \%parameter,
            mandatory_key_href     => \%mandatory_key,
            non_mandatory_key_href => \%non_mandatory_key,
            file_path              => $definitions_file,
        }
      )
};

## Then throw FATAL log message and croak
like(
    $trap->stderr,
    qr/Found\s+illegal\s+value/xms,
    q{Throw fatal log message for illegal value}
);

## Given a missing mandatory key
delete $parameter{family_id}{data_type};

trap {
    check_parameter_hash(
        {
            parameter_href         => \%parameter,
            mandatory_key_href     => \%mandatory_key,
            non_mandatory_key_href => \%non_mandatory_key,
            file_path              => $definitions_file,
        }
      )
};

## Then throw FATAL log message and croak
like(
    $trap->stderr,
    qr/Missing\s+mandatory\s+key/xms,
    q{Throw fatal log message for missing mandatory key}
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
