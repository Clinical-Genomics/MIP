package MIP::File::Format::Parameter;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ parse_definition_file };
}

## Constants
Readonly my $SPACE => q{ };

sub parse_definition_file {

## Function : Parse and check the definition parameters file
## Returns  : %parameter
## Arguments: $define_parameters_path            => File defining the supported parameters
##          : $non_mandatory_parameter_keys_path => Non mandatory keys
##          : $mandatory_parameter_keys_path     => Mandatory keys

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $define_parameters_path;
    my $non_mandatory_parameter_keys_path;
    my $mandatory_parameter_keys_path;

    my $tmpl = {
        define_parameters_path => {
            defined     => 1,
            required    => 1,
            store       => \$define_parameters_path,
            strict_type => 1,
        },
        non_mandatory_parameter_keys_path => {
            defined     => 1,
            required    => 1,
            store       => \$non_mandatory_parameter_keys_path,
            strict_type => 1,
        },
        mandatory_parameter_keys_path => {
            defined     => 1,
            required    => 1,
            store       => \$mandatory_parameter_keys_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Parameter qw{  check_parameter_hash };
    use MIP::File::Format::Yaml qw{ load_yaml };

    ## Loads a YAML file into an arbitrary hash and returns it.
    my %parameter = load_yaml( { yaml_file => $define_parameters_path, } );

    ## Load mandatory keys and values for parameters
    my %mandatory_key = load_yaml(
        {
            yaml_file => $mandatory_parameter_keys_path,
        }
    );

    ## Load non mandatory keys and values for parameters
    my %non_mandatory_key = load_yaml(
        {
            yaml_file => $non_mandatory_parameter_keys_path,
        }
    );

    ## Eval parameter hash
    check_parameter_hash(
        {
            file_path              => $define_parameters_path,
            non_mandatory_key_href => \%non_mandatory_key,
            mandatory_key_href     => \%mandatory_key,
            parameter_href         => \%parameter,
        }
    );

    return %parameter;
}

1;
