package MIP::Get::Parameter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      get_recipe_resources
      get_recipe_attributes
    };
}

sub get_recipe_attributes {

## Function : Return recipe attributes
## Returns  : $attribute | %attribute
## Arguments: $attribute      => Attribute key
##          : $parameter_href => Holds all parameters
##          : $recipe_name    => Recipe name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $attribute;
    my $parameter_href;
    my $recipe_name;

    my $tmpl = {
        attribute => {
            store       => \$attribute,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        recipe_name => {
            required    => 1,
            defined     => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Get attribute value
    if ( defined $attribute && $attribute ) {

        return $parameter_href->{$recipe_name}{$attribute};
    }

    ## Get recipe attribute hash
    return %{ $parameter_href->{$recipe_name} };
}

sub get_recipe_resources {

## Function : Return recipe resources
## Returns  : $recipe_resource | %recipe_resource
## Arguments: $active_parameter_href => The active parameters for this analysis hash {REF}
##          : $recipe_name           => Recipe name
##          : $recipe_resource       => Recipe parameter key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $recipe_name;
    my $resource;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        resource => {
            allow       => [qw{ core_number gpu_number load_env_ref memory time }],
            store       => \$resource,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ get_package_env_cmds };
    use MIP::Environment::Cluster qw{ check_recipe_memory_allocation };

    ## Initilize variable
    my @source_environment_cmds = get_package_env_cmds(
        {
            active_parameter_href => $active_parameter_href,
            package_name          => $recipe_name,
        }
    );

    my $core_number    = $active_parameter_href->{recipe_core_number}{$recipe_name};
    my $process_memory = $active_parameter_href->{recipe_memory}{$recipe_name};
    my $memory;

    ## Multiply memory with processes that are to be launched in the recipe
    if ( $process_memory and $core_number ) {
        $memory = $process_memory * $core_number;
    }
    ## Set default recipe memory allocation if it hasn't been specified
    elsif ( not $process_memory and $core_number ) {
        $memory = $core_number * $active_parameter_href->{core_ram_memory};
    }
    elsif ( not $process_memory and not $core_number ) {
        $memory = $active_parameter_href->{core_ram_memory};
    }
    else {
        $memory = $process_memory;
    }

    check_recipe_memory_allocation(
        {
            node_ram_memory          => $active_parameter_href->{node_ram_memory},
            recipe_memory_allocation => $memory,
        }
    );

    my %recipe_resource = (
        core_number  => $core_number,
        gpu_number   => $active_parameter_href->{recipe_gpu_number}{$recipe_name},
        load_env_ref => \@source_environment_cmds,
        memory       => $memory,
        time         => $active_parameter_href->{recipe_time}{$recipe_name},
    );

    ## Return specified recipe resource
    if ( defined $resource && $resource ) {
        return $recipe_resource{$resource};
    }

    ## Return recipe resource hash
    return %recipe_resource;

}

1;
