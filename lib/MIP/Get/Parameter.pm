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

## CPANM
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      get_gatk_intervals
      get_install_parameter_attribute
      get_package_source_env_cmds
      get_program_version
      get_recipe_resources
      get_recipe_attributes
    };
}

## Constants
Readonly my $TWO => 2;

sub get_gatk_intervals {

## Function : Generate and return interval hash
## Returns  : %gatk_intervals
## Arguments: $analysis_type         => Analysis type
##          : $contigs_ref           => Contigs to split in file
##          : $exome_target_bed_href => Exome target bed files lnked to sample ids
##          : $file_ending           => File ending to add {Optional}
##          : $filehandle            => Filehandle to write to
##          : $log                   => Log
##          : $max_cores_per_node    => Maximum core per node
##          : $outdirectory          => Outdirectory
##          : $reference_dir         => MIP reference directory
##          : $sample_id             => Sample_id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_type;
    my $contigs_ref;
    my $exome_target_bed_href;
    my $file_ending;
    my $filehandle;
    my $log;
    my $outdirectory;
    my $reference_dir;
    my $sample_id;

    ## Default(s)
    my $max_cores_per_node;

    my $tmpl = {
        analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$analysis_type,
            strict_type => 1,
        },
        contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$contigs_ref,
            strict_type => 1,
        },
        exome_target_bed_href => {
            store => \$exome_target_bed_href,
        },
        file_ending => {
            store       => \$file_ending,
            strict_type => 1,
        },
        filehandle => { store => \$filehandle, },
        log        => {
            store => \$log,
        },
        max_cores_per_node => {
            allow       => qr/ ^\d+$ /sxm,
            default     => 1,
            store       => \$max_cores_per_node,
            strict_type => 1,
        },
        outdirectory => {
            store       => \$outdirectory,
            strict_type => 1,
        },
        reference_dir => {
            store       => \$reference_dir,
            strict_type => 1,
        },
        sample_id => {
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ get_exome_target_bed_file };
    use MIP::Contigs qw{ generate_contig_interval_file };

    ## Store gatk interval
    my %gatk_intervals;

    if ( $analysis_type eq q{wes} ) {

        my $exome_target_bed_file = get_exome_target_bed_file(
            {
                exome_target_bed_href => $exome_target_bed_href,
                file_ending           => $file_ending,
                sample_id             => $sample_id,
            }
        );

        ## Generate contig specific interval_list and return hash with paths
        %gatk_intervals = generate_contig_interval_file(
            {
                contigs_ref           => $contigs_ref,
                exome_target_bed_file => $exome_target_bed_file,
                filehandle            => $filehandle,
                max_process_number    => $max_cores_per_node,
                outdirectory          => $outdirectory,
                reference_dir         => $reference_dir,
            }
        );
    }
    else {
        ## Key-value pairs are identical for WGS/WTS
        %gatk_intervals = map { $_ => [$_] } @{$contigs_ref};
    }

    return %gatk_intervals;
}

sub get_install_parameter_attribute {

## Function : Return parameter attribute from hash
## Returns  :
## Arguments: $parameter_href => Holds all parameters {REF}
##          : $parameter_name => Name of key to return

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $parameter_name;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        parameter_name => {
            store       => \$parameter_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Make sure that the supplied key exists
    croak(qq{Could not find parameter_name key: '$parameter_name' in hash})
      if ( not exists $parameter_href->{$parameter_name} );

    ## Hash attribute
    if ( ref $parameter_href->{$parameter_name} eq q{HASH} ) {

        return %{ $parameter_href->{$parameter_name} };
    }
    ## ARRAY attribute
    if ( ref $parameter_href->{$parameter_name} eq q{ARRAY} ) {

        return @{ $parameter_href->{$parameter_name} };
    }
    ## Scalar attribute
    return $parameter_href->{$parameter_name};
}

sub get_package_source_env_cmds {

## Function : Get package source environment commands
## Returns  : @source_environment_cmds
## Arguments: $active_parameter_href => The active parameters for this analysis hash {REF}
##          : $package_name          => Package name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $package_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        package_name => {
            defined     => 1,
            required    => 1,
            store       => \$package_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ get_package_env_attributes };
    use MIP::Environment::Manager qw{ get_env_method_cmds };

    ## Initilize variable
    my @source_environment_cmds;

    my ( $env_name, $env_method ) = get_package_env_attributes(
        {
            load_env_href => $active_parameter_href->{load_env},
            package_name  => $package_name,
        }
    );

    ## Could not find recipe within env
    if ( not $env_name ) {

        ## Fall back to MIPs MAIN env
        ( $env_name, $env_method ) = get_package_env_attributes(
            {
                load_env_href => $active_parameter_href->{load_env},
                package_name  => q{mip},
            }
        );
    }

    ## Get env load command
    my @env_method_cmds = get_env_method_cmds(
        {
            action     => q{load},
            env_name   => $env_name,
            env_method => $env_method,
        }
    );
    push @source_environment_cmds, @env_method_cmds;

    return @source_environment_cmds;
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

    use MIP::Environment::Cluster qw{ check_recipe_memory_allocation };

    ## Initilize variable
    my @source_environment_cmds = get_package_source_env_cmds(
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
