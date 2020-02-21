package MIP::Active_parameter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile splitpath };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw { any };

## MIPs lib/
use MIP::Constants qw{ $COMMA $DOT $LOG_NAME $PIPE $SINGLE_QUOTE $SPACE $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.15;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_parameter_files
      check_recipe_mode
      get_not_allowed_temp_dirs
      get_package_env_attributes
      get_user_supplied_pedigree_parameter
      parse_program_executables
      parse_recipe_resources
      set_binary_path
      set_default_analysis_type
      set_default_conda_path
      set_default_human_genome
      set_default_infile_dirs
      set_default_parameter
      set_default_pedigree_fam_file
      set_default_program_test_file
      set_default_reference_dir
      set_default_reference_info_file
      set_default_store_file
      set_default_temp_directory
      set_default_uninitialized_parameter
      set_default_vcfparser_select_file
      set_exome_target_bed
      set_gender_sample_ids
      set_parameter_reference_dir_path
      set_pedigree_sample_id_parameter
      set_recipe_resource
      update_recipe_mode_for_start_with_option
      update_reference_parameters
      update_to_absolute_path
    };
}

sub check_parameter_files {

## Function : Checks that files/directories files exists
## Returns  :
## Arguments: $active_parameter_href   => Holds all set parameter for analysis
##          : $associated_recipes_ref  => The parameters recipe(s) {REF}
##          : $build_status            => Build status of parameter
##          : $case_id                 => Case_id
##          : $consensus_analysis_type => Consensus analysis type for checking e.g. WGS specific files
##          : $parameter_exists_check  => Check if intendend file exists in reference directory
##          : $parameter_name          => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $associated_recipes_ref;
    my $build_status;
    my $consensus_analysis_type;
    my $parameter_exists_check;
    my $parameter_name;

    ## Default(s)
    my $case_id;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        associated_recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$associated_recipes_ref,
            strict_type => 1,
        },
        build_status => {
            store       => \$build_status,
            strict_type => 1,
        },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
            strict_type => 1,
        },
        consensus_analysis_type => {
            store       => \$consensus_analysis_type,
            strict_type => 1,
        },
        parameter_exists_check => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_exists_check,
            strict_type => 1,
        },
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path
      qw{ check_filesystem_objects_existance check_filesystem_objects_and_index_existance };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %only_wgs = ( gatk_genotypegvcfs_ref_gvcf => 1, );

    ## Do nothing since parameter is not required unless exome mode is enabled
    return
      if ( exists $only_wgs{$parameter_name}
        && $consensus_analysis_type =~ / wgs /xsm );

    ## Check all recipes that use parameter
  ASSOCIATED_RECIPE:
    foreach my $associated_recipe ( @{$associated_recipes_ref} ) {

        ## Active associated recipe
        my $associated_recipe_name = $active_parameter_href->{$associated_recipe};

        ## Active parameter
        my $active_parameter = $active_parameter_href->{$parameter_name};

        ## Only check active associated recipes parameters
        next ASSOCIATED_RECIPE if ( not $associated_recipe_name );

        ## Only check active parameters
        next ASSOCIATED_RECIPE if ( not defined $active_parameter );

        if ( ref $active_parameter eq q{ARRAY} ) {

            ## Get path for array elements
          PATH:
            foreach my $path ( @{ $active_parameter_href->{$parameter_name} } ) {

                check_filesystem_objects_and_index_existance(
                    {
                        is_build_file  => $build_status,
                        object_name    => $path,
                        object_type    => $parameter_exists_check,
                        parameter_name => $parameter_name,
                        path           => $path,
                    }
                );
            }
            return;
        }
        elsif ( ref $active_parameter eq q{HASH} ) {

            ## Get path for hash keys
          PATH:
            for my $path ( keys %{ $active_parameter_href->{$parameter_name} } ) {

                check_filesystem_objects_and_index_existance(
                    {
                        is_build_file  => $build_status,
                        object_name    => $path,
                        object_type    => $parameter_exists_check,
                        parameter_name => $parameter_name,
                        path           => $path,
                    }
                );
            }
            return;
        }

        ## File
        my $path = $active_parameter_href->{$parameter_name};

        check_filesystem_objects_and_index_existance(
            {
                is_build_file  => $build_status,
                object_name    => $path,
                object_type    => $parameter_exists_check,
                parameter_name => $parameter_name,
                path           => $path,
            }
        );
        return;
    }
    return;
}

sub check_recipe_mode {

## Function : Check correct value for recipe mode in MIP
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_href        => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Parameter qw{ get_cache };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Set allowed values
    my %is_allowed = map { $_ => 1 } ( 0 .. 2 );

    my @recipes = get_cache(
        {
            parameter_href => $parameter_href,
            parameter_name => q{recipe},
        }
    );

  RECIPE:
    foreach my $recipe (@recipes) {

        my $err_msg = q{Recipe: } . $recipe . q{ does not exist in %active_parameters};
        croak($err_msg) if ( not exists $active_parameter_href->{$recipe} );

        ## Unpack
        my $recipe_mode = $active_parameter_href->{$recipe};

        next RECIPE if ( $is_allowed{$recipe_mode} );

        #If not an allowed value in active parameters
        $log->fatal(
            $SINGLE_QUOTE
              . $active_parameter_href->{$recipe}
              . q{' Is not an allowed mode for recipe '--}
              . $recipe
              . q{'. Set to: }
              . join $PIPE,
            ( sort keys %is_allowed )
        );
        exit 1;
    }
    return 1;
}

sub get_not_allowed_temp_dirs {

## Function : Get paths that should not be set by mistake to temp_dir
## Returns  : @not_allowed_temp_dirs
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @is_not_allowed_temp_dirs = (
        $active_parameter_href->{cluster_constant_path},
        $active_parameter_href->{outdata_dir},
        $active_parameter_href->{outscript_dir},
        $active_parameter_href->{reference_dir},
    );
    return @is_not_allowed_temp_dirs;
}

sub get_package_env_attributes {

## Function : Get environment name and method for package (recipe, program or MIP)
## Returns  : $env_name, $env_method or "undef"
## Arguments: $load_env_href => Load env hash defining environments for packages {REF}
##          : $package_name  => Package name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $load_env_href;
    my $package_name;

    my $tmpl = {
        load_env_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$load_env_href,
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

  ENV:
    foreach my $env_name ( keys %{$load_env_href} ) {

        next ENV if ( not exists $load_env_href->{$env_name}{$package_name} );

        ## Found package_name within env
        ## Unpack
        my $env_method = $load_env_href->{$env_name}{method};
        return $env_name, $env_method;
    }
    return;
}

sub get_user_supplied_pedigree_parameter {

## Function : Detect if user supplied info on parameters otherwise collected from pedigree
## Returns  : %is_user_supplied - Hash where 1=user input and 0=no user input
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Define what should be checked
    my %is_user_supplied = (
        analysis_type         => 0,
        dna_sample_id         => 0,
        exome_target_bed      => 0,
        expected_coverage     => 0,
        sample_ids            => 0,
        supported_capture_kit => 0,
        time_point            => 0,
    );

    ## Detect user supplied info
  USER_PARAMETER:
    foreach my $parameter ( keys %is_user_supplied ) {

        ## If hash and supplied
        if ( ref $active_parameter_href->{$parameter} eq q{HASH}
            && keys %{ $active_parameter_href->{$parameter} } )
        {

            $is_user_supplied{$parameter} = 1;
        }
        elsif ( ref $active_parameter_href->{$parameter} eq q{ARRAY}
            && @{ $active_parameter_href->{$parameter} } )
        {
            ## If array and supplied
            $is_user_supplied{$parameter} = 1;
        }
        elsif ( defined $active_parameter_href->{$parameter}
            and not ref $active_parameter_href->{$parameter} )
        {

            ## If scalar and supplied
            $is_user_supplied{$parameter} = 1;
        }
    }
    return %is_user_supplied;
}

sub parse_program_executables {

## Function : Checking commands in your path and executable
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_href        => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ set_binary_path };
    use MIP::Environment::Path qw{ check_binary_in_path };

  PARAMETER:
    foreach my $parameter_name ( keys %{$active_parameter_href} ) {

        ## Only check path(s) for parameters with "type" key
        next PARAMETER
          if ( not exists $parameter_href->{$parameter_name}{type} );

        ## Only check path(s) for parameters with type value eq "recipe"
        next PARAMETER
          if ( not $parameter_href->{$parameter_name}{type} eq q{recipe} );

        ## Only check path(s) for active recipes
        next PARAMETER if ( not $active_parameter_href->{$parameter_name} );

        ## Alias
        my $program_executables_ref =
          \@{ $parameter_href->{$parameter_name}{program_executables} };

      PROGRAM:
        foreach my $program ( @{$program_executables_ref} ) {

            my $binary_path = check_binary_in_path(
                {
                    active_parameter_href => $active_parameter_href,
                    binary                => $program,
                    program_name          => $parameter_name,
                }
            );

            ## Set to use downstream
            set_binary_path(
                {
                    active_parameter_href => $active_parameter_href,
                    binary                => $program,
                    binary_path           => $binary_path,
                }
            );
        }
    }
    return;
}

sub parse_recipe_resources {

## Function : Check core number and memory requested against environment provisioned
## Returns  : 1
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Cluster
      qw{ check_max_core_number check_recipe_memory_allocation };

    ## Check that the recipe core number do not exceed the maximum per node
    foreach my $recipe_name ( keys %{ $active_parameter_href->{recipe_core_number} } ) {

        ## Limit number of cores requested to the maximum number of cores available per node
        $active_parameter_href->{recipe_core_number}{$recipe_name} =
          check_max_core_number(
            {
                max_cores_per_node => $active_parameter_href->{max_cores_per_node},
                core_number_requested =>
                  $active_parameter_href->{recipe_core_number}{$recipe_name},
            }
          );
    }

    ## Check that the recipe memory do not exceed the maximum per node
    ## Limit recipe_memory to node max memory if required
    foreach my $recipe_name ( keys %{ $active_parameter_href->{recipe_memory} } ) {

        $active_parameter_href->{recipe_memory}{$recipe_name} =
          check_recipe_memory_allocation(
            {
                node_ram_memory => $active_parameter_href->{node_ram_memory},
                recipe_memory_allocation =>
                  $active_parameter_href->{recipe_memory}{$recipe_name},
            }
          );
    }

    return 1;
}

sub set_binary_path {

## Function : Set binary path to active parameters
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $binary                => Binary to set
##          : $binary_path           => Path to binary

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $binary;
    my $binary_path;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        binary => { defined => 1, required => 1, store => \$binary, strict_type => 1, },
        binary_path => { required => 1, store => \$binary_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( not defined $binary_path );

    $active_parameter_href->{binary_path}{$binary} = $binary_path;
    return;
}

sub set_default_analysis_type {

## Function : Set default analysis type to active parameters
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    map { $active_parameter_href->{analysis_type}{$_} = q{wgs} }
      @{ $active_parameter_href->{sample_ids} };
    return;
}

sub set_default_conda_path {

## Function : Set default conda path to active parameters
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $conda_path            => Conda bin file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $conda_path;
    my $bin_file;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        bin_file => {
            default     => q{conda},
            store       => \$bin_file,
            strict_type => 1,
        },
        conda_path => { defined => 1, required => 1, store => \$conda_path, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Path qw{ get_conda_path };

    ## Set conda path
    $active_parameter_href->{$conda_path} =
      get_conda_path( { bin_file => $bin_file, } );

    if (   not $active_parameter_href->{$conda_path}
        or not -d $active_parameter_href->{$conda_path} )
    {

        croak(q{Failed to find default conda path});
    }
    return;
}

sub set_default_human_genome {

## Function : Set default human genome reference to active parameters
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $parameter_name        => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        parameter_name => { defined => 1, required => 1, store => \$parameter_name, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Now we now what human genome reference to build from
    $active_parameter_href->{$parameter_name} =
      $active_parameter_href->{human_genome_reference};

    return;
}

sub set_default_infile_dirs {

## Function : Set default infile dirs to active parameters
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Build default for infile_dirs
  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        if ( not exists $active_parameter_href->{analysis_type}{$sample_id} ) {

            set_default_analysis_type(
                {
                    active_parameter_href => $active_parameter_href,
                }
            );
        }
        my $path = catfile(
            $active_parameter_href->{cluster_constant_path},
            $active_parameter_href->{case_id},
            $active_parameter_href->{analysis_type}{$sample_id},
            $sample_id,
            q{fastq}
        );

        $active_parameter_href->{infile_dirs}{$path} = $sample_id;
    }
    return;
}

sub set_default_pedigree_fam_file {

## Function : Set default pedigree_fam_file to active parameters
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set pedigree fam file
    $active_parameter_href->{pedigree_fam_file} = catfile(
        $active_parameter_href->{outdata_dir},
        $active_parameter_href->{case_id},
        $active_parameter_href->{case_id} . $DOT . q{fam}
    );
    return;
}

sub set_default_program_test_file {

## Function : Set default path to file with program test commands
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( $active_parameter_href->{program_test_file} );

    $active_parameter_href->{program_test_file} =
      catfile( $Bin, qw{templates program_test_cmds.yaml } );

    return;
}

sub set_default_reference_dir {

## Function : Set default reference dir to active parameters
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set reference dir to current working dir
    $active_parameter_href->{reference_dir} = cwd();
    return;
}

sub set_default_reference_info_file {

## Function : Set default reference_info_file
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set reference info file
    $active_parameter_href->{reference_info_file} =
      catfile( $active_parameter_href->{outdata_dir}, q{reference_info.yaml} );
    return;
}

sub set_default_store_file {

## Function : Set default store_file to active parameters
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set store file
    $active_parameter_href->{store_file} = catfile( $active_parameter_href->{outdata_dir},
        $active_parameter_href->{case_id} . $UNDERSCORE . q{deliverables.yaml} );
    return;
}

sub set_default_temp_directory {

## Function : Set default temp directory to active parameters
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Mip download
    if ( exists $active_parameter_href->{download_pipeline_type} ) {

        $active_parameter_href->{temp_directory} =
          catfile( cwd(), qw{ mip_download $SLURM_JOB_ID } );
        return;
    }

    ## Mip analyse
    $active_parameter_href->{temp_directory} =
      catfile( $active_parameter_href->{outdata_dir}, q{$SLURM_JOB_ID} );

    return;
}

sub set_default_uninitialized_parameter {

## Function : Initiate hash keys for install
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $parameter_name        => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        parameter_name => { defined => 1, required => 1, store => \$parameter_name, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( exists $active_parameter_href->{$parameter_name} );

    $active_parameter_href->{$parameter_name} = [];

    return;
}

sub set_default_vcfparser_select_file {

## Function : Set default vcfparser select file to active parameters
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $parameter_name        => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        parameter_name => { defined => 1, required => 1, store => \$parameter_name, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Build default for (sv_)vcfparser select file
    my $path = catfile(
        $active_parameter_href->{cluster_constant_path},
        $active_parameter_href->{case_id},
        q{gene_panels.bed}
    );

    $active_parameter_href->{$parameter_name} = $path;
    return;
}

sub set_exome_target_bed {

## Function : Set exome target bed parameter in active_parameter hash
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $exome_target_bed_file => Exome target bed file to set
##          : $sample_id_string      => Sample id string to attach to exome bed file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $exome_target_bed_file;
    my $sample_id_string;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        exome_target_bed_file => {
            defined     => 1,
            required    => 1,
            store       => \$exome_target_bed_file,
            strict_type => 1,
        },
        sample_id_string => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id_string,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Add sample_ids as string to exome_target_bed_file
    $active_parameter_href->{exome_target_bed}{$exome_target_bed_file} =
      $sample_id_string;

    return;
}

sub set_default_parameter {

## Function : Set default parameter in active_parameter hash
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_name        => Parameter name to set default for
##          : $parameter_default     => Parameter default value to set (scalar|array_ref|hash_ref)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_name;
    my $parameter_default;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
        parameter_default => {
            defined  => 1,
            required => 1,
            store    => \$parameter_default,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Add default value for parameter
    ## Can be scalar|array_ref|hash_ref
    $active_parameter_href->{$parameter_name} = $parameter_default;

    return;
}

sub set_gender_sample_ids {

## Function : Set gender sample_ids of the current analysis and count each instance
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Sample_info qw{ get_pedigree_sample_id_attributes };

    my %gender_map = (
        1      => \&_set_male_gender,
        2      => \&_set_female_gender,
        female => \&_set_female_gender,
        male   => \&_set_male_gender,
        other  => \&_set_other_and_male_gender,
    );

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        my $sex = get_pedigree_sample_id_attributes(
            {
                attribute        => q{sex},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        if ( exists $gender_map{$sex} ) {

            $gender_map{$sex}->(
                {
                    active_parameter_href => $active_parameter_href,
                    sample_id             => $sample_id
                }
            );
            next SAMPLE_ID;
        }
        ## Set as other and male
        _set_other_and_male_gender(
            {
                active_parameter_href => $active_parameter_href,
                sample_id             => $sample_id
            }
        );

    }
    return;
}

sub set_parameter_reference_dir_path {

## Function : Set path for supplied reference(s) associated with parameter that should
##          : reside in the mip reference directory to full path.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_name        => Parameter to update

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Unpack
    my $reference_dir = $active_parameter_href->{reference_dir};

    # $parameter can be array_ref, hash_ref, point to file or undef
    my $parameter = $active_parameter_href->{$parameter_name};

    return if ( not defined $parameter );

    if ( ref $parameter eq q{ARRAY} ) {

      FILE:
        foreach my $file ( @{$parameter} ) {

            ## Split to restate
            my ( $volume, $directory, $file_name ) = splitpath($file);

            ## Update original element - works since array_ref
            ## Parameter element now stores path instead of file
            $file = catfile( $reference_dir, $file_name );
        }
        return;
    }
    elsif ( ref $parameter eq q{HASH} ) {

      FILE:
        foreach my $file ( keys %{$parameter} ) {

            ## Split to restate
            my ( $volume, $directory, $file_name ) = splitpath($file);

            ## Update original key with path and add potential annotation key
            ## by deleting original value (returns value deleted)
            my $path = catfile( $reference_dir, $file_name );
            $active_parameter_href->{$parameter_name}{$path} =
              delete $active_parameter_href->{$parameter_name}{$file};
        }
        return;
    }

    ## File
    ## Split to restate
    my ( $volume, $directory, $file_name ) =
      splitpath( $active_parameter_href->{$parameter_name} );

    ## Restate to allow for changing mip reference directory between runs
    $active_parameter_href->{$parameter_name} = $file_name;

    ## Update original value
    my $path = catfile( $reference_dir, $active_parameter_href->{$parameter_name} );
    $active_parameter_href->{$parameter_name} = $path;

    return;
}

sub set_pedigree_sample_id_parameter {

## Function : Set pedigree parameter in active_parameter hash
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $pedigree_key          => Pedigree key to set
##          : $pedigree_value        => Pedigree value to set
##          : $sample_id             => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $pedigree_key;
    my $pedigree_value;
    my $sample_id;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        pedigree_key => {
            defined     => 1,
            required    => 1,
            store       => \$pedigree_key,
            strict_type => 1,
        },
        pedigree_value => {
            defined     => 1,
            required    => 1,
            store       => \$pedigree_value,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Add value for sample_id using pedigree info
    $active_parameter_href->{$pedigree_key}{$sample_id} = $pedigree_value;

    return;
}

sub set_recipe_resource {

## Function : Set recipe resource allocation for specific recipe(s)
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %set_hash_key_map = (
        set_recipe_core_number => q{recipe_core_number},
        set_recipe_time        => q{recipe_time},
        set_recipe_memory      => q{recipe_memory},
    );

  HASH_KEY:
    while ( my ( $set_hash_key, $target_hash_key ) = each %set_hash_key_map ) {

      RECIPE:
        while ( my ( $recipe, $core_number ) =
            each %{ $active_parameter_href->{$set_hash_key} } )
        {

            $active_parameter_href->{$target_hash_key}{$recipe} = $core_number;
        }
    }
    return;
}

sub update_recipe_mode_for_start_with_option {

## Function : Update recipe mode depending on recipe to start with
## Returns  :
## Arguments: $active_parameter_href  => The active parameters for this analysis hash {REF}
##          : $recipes_ref            => Recipes in MIP
##          : $start_with_recipes_ref => Recipes to run

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $recipes_ref;
    my $start_with_recipes_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$recipes_ref,
            strict_type => 1,
        },
        start_with_recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$start_with_recipes_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Run mode
    my $run_mode = 1;

    ## Dry run mode
    my $dry_run_mode = 2;

  RECIPE:
    foreach my $recipe_name ( @{$recipes_ref} ) {

        next RECIPE if ( not $active_parameter_href->{$recipe_name} );

        ## If recipe is uppstream of start with recipe
        if ( not any { $_ eq $recipe_name } @{$start_with_recipes_ref} ) {

            ## Change recipe mode to dry run
            $active_parameter_href->{$recipe_name} = $dry_run_mode;
            next RECIPE;
        }
        ## Recipe or downstream dependency recipe
        # Change recipe mode to active
        $active_parameter_href->{$recipe_name} = $run_mode;
    }
    return;
}

sub update_reference_parameters {

## Function : Update reference parameters with mip_reference directory path
## Returns  :
## Arguments: $active_parameter_href  => Holds all set parameter for analysis
##          : $associated_recipes_ref => The parameters recipe(s) {REF}
##          : $parameter_name         => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $associated_recipes_ref;
    my $parameter_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        associated_recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$associated_recipes_ref,
            strict_type => 1,
        },
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Check all recipes that use parameter
  ASSOCIATED_RECIPE:
    foreach my $associated_recipe ( @{$associated_recipes_ref} ) {

        my $recipe_name = $active_parameter_href->{$associated_recipe};

        ## Only check active recipes parameters
        next ASSOCIATED_RECIPE if ( not $recipe_name );

        ## Update path for supplied reference(s) associated with
        ## parameter that should reside in the mip reference directory to full path
        set_parameter_reference_dir_path(
            {
                active_parameter_href => $active_parameter_href,
                parameter_name        => $parameter_name,
            }
        );

        ## Only need to perform update once per parameter
        return;
    }
    return;
}

sub update_to_absolute_path {

## Function : Change relative path to absolute path in active_parameters for
##            parameters with key value pair "update_path: absolute path" in parameter hash
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_href        => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw{ get_absolute_path };
    use MIP::Parameter qw{ set_cache };

    ## Adds dynamic aggregate information from definitions to parameter hash
    # Collect all path that should be made absolute
    set_cache(
        {
            aggregates_ref => [q{update_path:absolute_path}],
            parameter_href => $parameter_href,
        }
    );

  DYNAMIC_PARAMETER:
    foreach my $parameter_name ( @{ $parameter_href->{cache}{absolute_path} } ) {

        next DYNAMIC_PARAMETER
          if ( not exists $active_parameter_href->{$parameter_name} );

        next DYNAMIC_PARAMETER
          if ( not defined $active_parameter_href->{$parameter_name} );

        if ( ref $active_parameter_href->{$parameter_name} eq q{ARRAY} ) {

          VALUE:
            foreach my $parameter_value ( @{ $active_parameter_href->{$parameter_name} } )
            {

                $parameter_value = get_absolute_path(
                    {
                        parameter_name => $parameter_name,
                        path           => $parameter_value,
                    }
                );
            }
            next DYNAMIC_PARAMETER;
        }
        if ( ref $active_parameter_href->{$parameter_name} eq q{HASH} ) {

            ## Alias
            my $parameter_name_href = $active_parameter_href->{$parameter_name};

            ## Cannot use each since we are updating key within loop
          KEY:
            foreach my $key ( keys %{$parameter_name_href} ) {

                ## Return absolute path for supplied key path or croaks and exists if path does not exists
                my $updated_key = get_absolute_path(
                    {
                        parameter_name => $parameter_name,
                        path           => $key,
                    }
                );

                $parameter_name_href->{$updated_key} =
                  delete $parameter_name_href->{$key};
            }
            next DYNAMIC_PARAMETER;
        }
        ## Scalar

        $active_parameter_href->{$parameter_name} = get_absolute_path(
            {
                parameter_name => $parameter_name,
                path           => $active_parameter_href->{$parameter_name},
            }
        );
    }
    return;
}

sub _set_female_gender {

## Function : Set female gender sample_ids of the current analysis and count each instance
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $sample_id             => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_id;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    $active_parameter_href->{found_female}++;
    push @{ $active_parameter_href->{gender}{females} }, $sample_id;
    return;
}

sub _set_male_gender {

## Function : Set male gender sample_ids of the current analysis and count each instance
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $sample_id             => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_id;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    $active_parameter_href->{found_male}++;
    push @{ $active_parameter_href->{gender}{males} }, $sample_id;
    return;
}

sub _set_other_and_male_gender {

## Function : Set other gender sample_ids of the current analysis and count each instance
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $sample_id             => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_id;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Include since it might be male to enable analysis of Y. For WGS estimation of gender
    ## will be performed from fastq reads
    $active_parameter_href->{found_male}++;

    $active_parameter_href->{found_other}++;
    push @{ $active_parameter_href->{gender}{others} }, $sample_id;
    return;
}

1;
