#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.04;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Parameter}          => [qw{ check_parameter_hash }],
        q{MIP::File::Format::Yaml} => [qw{ load_yaml }],
        q{MIP::Test::Fixtures}     => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parameter qw{ check_parameter_hash };
use MIP::File::Format::Yaml qw{ load_yaml };

diag(   q{Test check_parameter_hash from Parameter.pm v}
      . $MIP::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $definitions_file = catfile( $Bin, qw{ data test_data define_parameters.yaml } );

## Loads a YAML file into an arbitrary hash and returns it.
my %parameter = load_yaml( { yaml_file => $definitions_file, } );

## Load required keys and values for parameters
my %required = load_yaml(
    {
        yaml_file => catfile( dirname($Bin), qw{ definitions required_parameters.yaml } ),
    }
);

## Load non required keys and values for parameters
my %not_required = load_yaml(
    {
        yaml_file =>
          catfile( dirname($Bin), qw{ definitions not_required_parameters.yaml } ),

    }
);

## Given valid input
my $is_ok = check_parameter_hash(
    {
        file_path         => $definitions_file,
        not_required_href => \%not_required,
        parameter_href    => \%parameter,
        required_href     => \%required,
    }
);

## Then return true
ok( $is_ok, q{No errors in parameters detected} );

## Given wrong data type, when SCALAR
$parameter{case_id}{data_type} = [q{wrong_data_type}];

trap {
    check_parameter_hash(
        {
            file_path         => $definitions_file,
            not_required_href => \%not_required,
            parameter_href    => \%parameter,
            required_href     => \%required,
        }
    )
};

## Then throw FATAL log message and croak
like( $trap->stderr, qr/ARRAY/xms,
    q{Throw fatal log message for wrong data type - array} );

## Reset parameter
$parameter{case_id}{data_type} = q{SCALAR};

## Given wrong data type, when ARRAY or HASH
$parameter{case_id}{associated_recipe} = q{not_an_array};

trap {
    check_parameter_hash(
        {
            file_path         => $definitions_file,
            not_required_href => \%not_required,
            parameter_href    => \%parameter,
            required_href     => \%required,
        }
    )
};

## Then throw FATAL log message and croak
like( $trap->stderr, qr/SCALAR/xms,
    q{Throw fatal log message for wrong data type - scalar} );

## Reset parameter
$parameter{case_id}{associated_recipe} = [qw{ mip }];

## Given not allowed value
$parameter{case_id}{data_type} = q{not_valid_value};

trap {
    check_parameter_hash(
        {
            file_path         => $definitions_file,
            not_required_href => \%not_required,
            parameter_href    => \%parameter,
            required_href     => \%required,
        }
    )
};

## Then throw FATAL log message and croak
like(
    $trap->stderr,
    qr/Found\s+illegal\s+value/xms,
    q{Throw fatal log message for illegal value}
);

## Given a missing required key
delete $parameter{case_id}{data_type};

trap {
    check_parameter_hash(
        {
            file_path         => $definitions_file,
            not_required_href => \%not_required,
            parameter_href    => \%parameter,
            required_href     => \%required,
        }
    )
};

## Then throw FATAL log message and croak
like(
    $trap->stderr,
    qr/Missing\s+required\s+key/xms,
    q{Throw fatal log message for missing required key}
);

done_testing();
