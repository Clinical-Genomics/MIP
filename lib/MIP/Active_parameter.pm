package MIP::Active_parameter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile splitpath };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw { any };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $COLON $COMMA $DOT $LOG_NAME $NEWLINE $PIPE $SINGLE_QUOTE $SPACE $TAB $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      add_gender
      add_recipe_bind_paths
      check_load_env_packages
      check_parameter_files
      check_recipe_mode
      check_sample_id_in_hash_parameter
      check_sample_id_in_hash_parameter_path
      get_active_parameter_attribute
      get_binary_path
      get_exome_target_bed_file
      get_matching_values_key
      get_not_allowed_temp_dirs
      get_package_env_attributes
      get_package_env_cmds
      get_recipe_resources
      get_user_supplied_pedigree_parameter
      parse_infiles
      parse_recipe_resources
      parse_vep_plugin
      remove_sample_id_from_gender
      set_conda_paths
      set_default_analysis_type
      set_default_container_config_file
      set_default_human_genome
      set_default_infile_dirs
      set_default_parameter
      set_default_pedigree_fam_file
      set_default_program_test_file
      set_default_reference_dir
      set_default_reference_info_file
      set_default_store_file
      set_default_qccollect_store_metrics_outfile
      set_default_temp_directory
      set_default_transcript_annotation
      set_default_uninitialized_parameter
      set_default_vcfparser_select_file
      set_exome_target_bed
      set_gender_estimation
      set_gender_sample_ids
      set_include_y
      set_load_env_environment
      set_parameter_reference_dir_path
      set_pedigree_sample_id_parameter
      set_programs_for_installation
      set_recipe_mode
      set_recipe_resource
      set_vcfparser_outfile_counter
      update_dynamic_config_parameters
      update_recipe_mode_for_start_with_option
      update_recipe_mode_with_dry_run_all
      update_reference_parameters
      update_to_absolute_path
      update_with_dynamic_config_parameters
      write_references
    };
}

sub add_gender {

## Function : Add sample_id to list of genders of the current analysis
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $gender                => Gender to add to
##          : $sample_id             => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $gender;
    my $sample_id;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        gender => {
            allow       => [qw{ females males others }],
            defined     => 1,
            required    => 1,
            store       => \$gender,
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

    push @{ $active_parameter_href->{gender}{$gender} }, $sample_id;
    return;
}

sub add_recipe_bind_paths {

## Function : Add recipe specific bind paths to defaults
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $export_bind_paths_ref => Bind paths to add to
##          : $recipe_name           => Recipe name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $export_bind_paths_ref;
    my $recipe_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        export_bind_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$export_bind_paths_ref,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( not $active_parameter_href->{recipe_bind_path}{$recipe_name} );

    push @{$export_bind_paths_ref}, @{ $active_parameter_href->{recipe_bind_path}{$recipe_name} };

    return;
}

sub check_load_env_packages {

## Function : Check that package name name are included in MIP as either "mip", "recipe" or "program_executables"
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

    ## Constants
    Readonly my @LOAD_ENV_KEYS => qw{ installation method mip };

    my @program_executables = get_cache(
        {
            parameter_href => $parameter_href,
            parameter_name => q{program_executables},
        }
    );

    my @recipes = get_cache(
        {
            parameter_href => $parameter_href,
            parameter_name => q{recipe},
        }
    );

    ## Allowed packages/keywords in load_env section of config
    my @allowed_packages = ( @program_executables, @LOAD_ENV_KEYS, @recipes, );

  ENV:
    foreach my $env ( keys %{ $active_parameter_href->{load_env} } ) {

      PACKAGE:
        foreach my $package ( keys %{ $active_parameter_href->{load_env}{$env} } ) {

            ## is program executable, installation, method, mip or recipe_name
            next PACKAGE if ( any { $_ eq $package } @allowed_packages );

            my $err_msg =
                q{Could not find load_env package: '}
              . $package
              . q{' in MIP as either recipe or program_executables};
            croak($err_msg);
        }
    }
    return 1;
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

sub check_sample_id_in_hash_parameter {

## Function : Check sample_id provided in hash parameter is included in the
##          : analysis and only represented once
## Returns  : 1
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_href        => Holds all parameters {REF}
##          : $parameter_names_ref   => Parameter name list {REF}
##          : $sample_ids_ref        => Array to loop in for parameter {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_names_ref;
    my $parameter_href;
    my $sample_ids_ref;

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
        parameter_names_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$parameter_names_ref,
            strict_type => 1,
        },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Parameter qw{ get_parameter_attribute };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

  PARAMETER:
    foreach my $parameter_name ( @{$parameter_names_ref} ) {

        ## Skip undef parameters in current analysis
        next PARAMETER
          if ( not defined $active_parameter_href->{$parameter_name} );

      SAMPLE_ID:
        foreach my $sample_id ( @{$sample_ids_ref} ) {

            my $sample_id_value = get_active_parameter_attribute(
                {
                    active_parameter_href => $active_parameter_href,
                    attribute             => $sample_id,
                    parameter_name        => $parameter_name,
                }
            );
            my $is_mandatory = get_parameter_attribute(
                {
                    attribute      => q{mandatory},
                    parameter_href => $parameter_href,
                    parameter_name => $parameter_name,
                }
            );

            ## Check that a value is defined
            next SAMPLE_ID if ( defined $sample_id_value );

            ## Not mandatory - skip
            next PARAMETER
              if ( defined $is_mandatory and $is_mandatory eq q{no} );

            my %parameter_name_hash = get_active_parameter_attribute(
                {
                    active_parameter_href => $active_parameter_href,
                    parameter_name        => $parameter_name,
                }
            );
            my $parameter_name_sample_ids = join $COMMA . $SPACE, ( keys %parameter_name_hash );

            $log->fatal( q{Could not find value for }
                  . $sample_id
                  . q{ for parameter '--}
                  . $parameter_name
                  . $SINGLE_QUOTE
                  . q{. Provided sample_ids for parameter are: }
                  . $parameter_name_sample_ids );
            exit 1;
        }
    }
    return 1;
}

sub check_sample_id_in_hash_parameter_path {

## Function : Check sample_id provided in hash path parameter is included in the
##          : analysis and only represented once
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_names_ref   => Parameter name list {REF}
##          : $sample_ids_ref        => Array to loop in for parameter {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_names_ref;
    my $sample_ids_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        parameter_names_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$parameter_names_ref,
            strict_type => 1,
        },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

  PARAMETER:
    foreach my $parameter_name ( @{$parameter_names_ref} ) {

        # Hash to test duplicate sample_ids later
        my %seen;

        my %parameter_name_hash = get_active_parameter_attribute(
            {
                active_parameter_href => $active_parameter_href,
                parameter_name        => $parameter_name,
            }
        );

      PATH:
        while ( my ( $path, $sample_ids_str ) = each %parameter_name_hash ) {

            ## Get sample ids for parameter from string
            my @parameter_samples = split $COMMA, $sample_ids_str;

          SAMPLE_ID:
            foreach my $sample_id (@parameter_samples) {

                # Increment instance to check duplicates later
                $seen{$sample_id}++;

                ## Check sample_id are unique
                next SAMPLE_ID if ( $seen{$sample_id} < 2 );

                $log->fatal( q{Sample_id: }
                      . $sample_id
                      . q{ is not uniqe in '--}
                      . $parameter_name . q{': }
                      . $path . q{=}
                      . $sample_ids_str, );
                exit 1;
            }
        }

        ## Check all sample ids are present in parameter string
      SAMPLE_ID:
        foreach my $sample_id ( @{$sample_ids_ref} ) {

            ## If sample_id is not present in parameter_name hash
            next SAMPLE_ID if ( exists $seen{$sample_id} );

            my $seen_sample_ids = join $COMMA . $SPACE, ( keys %seen );

            $log->fatal( q{Could not detect }
                  . $sample_id
                  . q{ for '--}
                  . $parameter_name
                  . q{'. Provided sample_ids are: }
                  . $seen_sample_ids, );
            exit 1;
        }
    }
    return 1;
}

sub get_active_parameter_attribute {

## Function : Get active parameter attribute
## Returns  : $attribute
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $attribute             => Attribute to return
##          : $parameter_name        => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $attribute;
    my $parameter_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        attribute => {
            store       => \$attribute,
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

    ## Return entire parameter attribute hash if no specific attribute key supplied
    return %{ $active_parameter_href->{$parameter_name} } if ( not defined $attribute );

    ## Unpack
    my $parameter_attribute = $active_parameter_href->{$parameter_name}{$attribute};

    return if ( not defined $parameter_attribute );

    return @{$parameter_attribute} if ( ref $parameter_attribute eq q{ARRAY} );

    return %{$parameter_attribute} if ( ref $parameter_attribute eq q{HASH} );

    ## Return scalar parameter attribute value
    return $parameter_attribute;
}

sub get_binary_path {

## Function : Get path to binary
## Returns  : $active_parameter_href->{binary_path}{$binary}
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $binary                => Binary name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $binary;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        binary => {
            defined     => 1,
            required    => 1,
            store       => \$binary,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( not defined $active_parameter_href->{binary_path}{$binary} );

    return $active_parameter_href->{binary_path}{$binary};
}

sub get_exome_target_bed_file {

## Function : Get exome_target_bed file for specific sample_id and add file_ending if supplied
## Returns  : $exome_target_bed_file
## Arguments: $exome_target_bed_href => Exome target bed files lnked to sample ids
##          : $file_ending           => File ending to add to file
##          : $sample_id             => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $exome_target_bed_href;
    my $file_ending;
    my $sample_id;

    my $tmpl = {
        exome_target_bed_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$exome_target_bed_href,
            strict_type => 1,
        },
        file_ending => { store => \$file_ending, strict_type => 1, },
        sample_id   => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## To track sample_ids with capture kits
    my %seen;

  BED_FILE:
    while ( my ( $exome_target_bed_file, $sample_id_string ) = each %{$exome_target_bed_href} ) {

        my @capture_kit_samples = split $COMMA, $sample_id_string;

        @seen{@capture_kit_samples}++;

        ## If capture_kit sample_id is associated with exome_target_bed file
        if ( exists $seen{$sample_id} ) {

            ## Reset each iterator for next call
            scalar keys %{$exome_target_bed_href};

            return $exome_target_bed_file if ( not defined $file_ending );

            return $exome_target_bed_file . $file_ending;
        }
    }

    $log->fatal(
        q{Could not detect }
          . $sample_id
          . q{ in '-exome_target_bed' associated files in sub routine get_exome_target_bed_file},
        $NEWLINE
    );
    exit 1;
}

sub get_matching_values_key {

## Function : Return the key if the hash value exists
## Returns  : "key pointing to matched value"
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_name        => MIP parameter name
##          : $query_value           => Value to query in the hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_name;
    my $query_value;

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
        query_value => {
            defined     => 1,
            required    => 1,
            store       => \$query_value,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( not exists $active_parameter_href->{$parameter_name} );

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my @collapsed_values;

    ## Values are now keys and vice versa
    my %reversed = reverse %{ $active_parameter_href->{$parameter_name} };

  PAIR:
    while ( my ( $key, $value ) = each %{ $active_parameter_href->{$parameter_name} } ) {

        next PAIR if ( exists $reversed{$value} and $reversed{$value} eq $key );

        push @collapsed_values, $value;
    }
    if (@collapsed_values) {
        $log->fatal(qq{Found duplicated values in parameter: $parameter_name });
        $log->fatal( q{Duplicated values: } . join $SPACE, @collapsed_values );
        exit 1;
    }

    return $reversed{$query_value} if ( exists $reversed{$query_value} );

    return;
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
        $active_parameter_href->{cluster_constant_path}, $active_parameter_href->{outdata_dir},
        $active_parameter_href->{outscript_dir},         $active_parameter_href->{reference_dir},
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

sub get_package_env_cmds {

## Function : Get package environment commands
## Returns  : @env_method_cmds
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
    my $conda_init_path =
      catfile( $active_parameter_href->{conda_path}, qw{ etc profile.d conda.sh } );
    my @env_method_cmds = get_env_method_cmds(
        {
            action          => q{load},
            conda_init_path => $conda_init_path,
            env_name        => $env_name,
            env_method      => $env_method,
        }
    );
    return @env_method_cmds;
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
            allow       => [qw{ core_number gpu_number load_env_ref memory mode time }],
            store       => \$resource,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ get_package_env_cmds };
    use MIP::Environment::Cluster qw{ check_recipe_memory_allocation };

    ## Initilize variable
    my @environment_cmds = get_package_env_cmds(
        {
            active_parameter_href => $active_parameter_href,
            package_name          => $recipe_name,
        }
    );

    my $core_number     = $active_parameter_href->{recipe_core_number}{$recipe_name};
    my $process_memory  = $active_parameter_href->{recipe_memory}{$recipe_name};
    my $core_ram_memory = $active_parameter_href->{core_ram_memory};

    my $memory =
        ( $process_memory     and $core_number )     ? $process_memory * $core_number
      : ( not $process_memory and $core_number )     ? $core_number * $core_ram_memory
      : ( not $process_memory and not $core_number ) ? $core_ram_memory
      :                                                $process_memory;

    check_recipe_memory_allocation(
        {
            node_ram_memory          => $active_parameter_href->{node_ram_memory},
            recipe_memory_allocation => $memory,
        }
    );

    my %recipe_resource = (
        core_number  => $core_number,
        gpu_number   => $active_parameter_href->{recipe_gpu_number}{$recipe_name},
        load_env_ref => \@environment_cmds,
        memory       => $memory,
        mode         => $active_parameter_href->{$recipe_name},
        time         => $active_parameter_href->{recipe_time}{$recipe_name},
    );

    return $recipe_resource{$resource} if ($resource);

    return %recipe_resource;

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

sub parse_infiles {

## Function : Collects the ".fastq(.gz)" files from the supplied infiles directory. Checks if any files exist.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $file_info_href        => File info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ get_matching_values_key };
    use MIP::File_info qw{ set_infiles };
    use MIP::File::Path qw{ get_file_names };
    use MIP::Validate::Case qw{ check_infiles };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Collect input files governed by sample_ids
  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Return the key if the hash value exists
        my $infile_directory = get_matching_values_key(
            {
                active_parameter_href => $active_parameter_href,
                parameter_name        => q{infile_dirs},
                query_value           => $sample_id,
            }
        );

        _check_infile_directory(
            {
                infile_directory => $infile_directory,
                sample_id        => $sample_id,
            }
        );

        ## Get the file(s) from filesystem
        my @infiles = get_file_names(
            {
                file_directory   => $infile_directory,
                rule_name        => q{*.fastq*},
                rule_skip_subdir => q{original_fastq_files},
            }
        );

        ## Check infiles found and that they contain sample_id
        check_infiles(
            {
                infiles_ref      => \@infiles,
                infile_directory => $infile_directory,
                sample_id        => $sample_id,
            }
        );

        ## Set the infile features i.e. dir and files
        set_infiles(
            {
                file_info_href   => $file_info_href,
                infiles_ref      => \@infiles,
                infile_directory => $infile_directory,
                sample_id        => $sample_id,
            }
        );

        ## Broadcast to user
        $log->info(q{Reads from platform:});
        $log->info( q{Sample id: } . $sample_id );
        $log->info( $TAB . q{Inputfiles:} );

        ## Log each file from platform
      FILE:
        foreach my $file (@infiles) {

            # Indent for visability
            $log->info( $TAB x 2, $file );
        }
    }
    return 1;
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

    use MIP::Environment::Cluster qw{ check_max_core_number check_recipe_memory_allocation };

    ## Check that the recipe core number do not exceed the maximum per node
    foreach my $recipe_name ( keys %{ $active_parameter_href->{recipe_core_number} } ) {

        ## Limit number of cores requested to the maximum number of cores available per node
        $active_parameter_href->{recipe_core_number}{$recipe_name} = check_max_core_number(
            {
                max_cores_per_node    => $active_parameter_href->{max_cores_per_node},
                core_number_requested => $active_parameter_href->{recipe_core_number}{$recipe_name},
            }
        );
    }

    ## Check that the recipe memory do not exceed the maximum per node
    ## Limit recipe_memory to node max memory if required
    foreach my $recipe_name ( keys %{ $active_parameter_href->{recipe_memory} } ) {

        $active_parameter_href->{recipe_memory}{$recipe_name} = check_recipe_memory_allocation(
            {
                node_ram_memory          => $active_parameter_href->{node_ram_memory},
                recipe_memory_allocation => $active_parameter_href->{recipe_memory}{$recipe_name},
            }
        );
    }

    return 1;
}

sub parse_vep_plugin {

## Function : Parse VEP plugin options
## Returns  : 1
## Arguments: $active_parameter_href => Holds all set parameter for analysis
##          : $mip_vep_plugins_ref   => MIP VEP plugin parameter names

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $mip_vep_plugins_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        mip_vep_plugins_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$mip_vep_plugins_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Vep qw{ check_vep_plugin };

  PLUGIN_PARAM:
    foreach my $parameter_name ( @{$mip_vep_plugins_ref} ) {

        check_vep_plugin(
            {
                parameter_name       => $parameter_name,
                vep_plugin_href      => \%{ $active_parameter_href->{$parameter_name} },
                vep_plugins_dir_path => $active_parameter_href->{vep_plugins_dir_path},
            }
        );
    }
    return 1;
}

sub remove_sample_id_from_gender {

## Function : Remove sample_id from gender
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $gender                => Gender to remove sample id from
##          : $sample_id             => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $gender;
    my $sample_id;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        gender => {
            allow       => [qw{ females males others }],
            defined     => 1,
            required    => 1,
            store       => \$gender,
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

    @{ $active_parameter_href->{gender}{$gender} } =
      grep { not /$sample_id/xms } @{ $active_parameter_href->{gender}{$gender} };
    return;
}

sub set_conda_paths {

## Function : Set path to root conda dir and the conda_environment_path based on the root path
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $environment_name      => MIP conda environment name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $environment_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        environment_name => {
            defined     => 1,
            required    => 1,
            store       => \$environment_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Path qw{ get_conda_path };

    ## Get path to conda
    $active_parameter_href->{conda_path} = $active_parameter_href->{conda_path}
      // get_conda_path( {} );

    ## Set path to conda env
    $active_parameter_href->{conda_environment_path} =
      catdir( $active_parameter_href->{conda_path}, q{envs}, $environment_name );

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
    $active_parameter_href->{$parameter_name} = $active_parameter_href->{human_genome_reference};

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
            $sample_id, q{fastq}
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

sub set_default_qccollect_store_metrics_outfile {

## Function : Set default qccollect_store_metrics_outfile to active parameters
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
    $active_parameter_href->{qccollect_store_metrics_outfile} =
      catfile( $active_parameter_href->{outdata_dir},
        $active_parameter_href->{case_id} . $UNDERSCORE . q{metrics_deliverables.yaml} );
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
    if ( exists $active_parameter_href->{download_pipeline} ) {

        $active_parameter_href->{temp_directory} =
          catfile( cwd(), qw{ mip_download $SLURM_JOB_ID } );
        return;
    }

    ## Mip analyse
    $active_parameter_href->{temp_directory} =
      catfile( $active_parameter_href->{outdata_dir}, q{$SLURM_JOB_ID} );

    return;
}

sub set_default_transcript_annotation {

## Function : Set default transcript annotation to active parameters
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

    ## Now we now what transcript annotation reference to build from
    $active_parameter_href->{$parameter_name} = $active_parameter_href->{transcript_annotation};

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
    $active_parameter_href->{exome_target_bed}{$exome_target_bed_file} = $sample_id_string;

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

sub set_load_env_environment {

## Function : Set load_env environment name based on supplied environment name
## Returns  :
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

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments};

    ## Unpack
    my $env_name = $active_parameter_href->{environment_name};

    return 0 if ( not defined $env_name );

    if ( keys %{ $active_parameter_href->{load_env} } != 1 ) {

        ## Retrieve logger object
        my $log = Log::Log4perl->get_logger($LOG_NAME);

        $log->warn(
            q{Could not use } . $env_name . q{ as there are multiple environments to update} );
        return;
    }

  ENV:
    foreach my $config_env_name ( keys %{ $active_parameter_href->{load_env} } ) {

        ## Set new environment name
        $active_parameter_href->{load_env}{$env_name} =
          delete $active_parameter_href->{load_env}{$config_env_name};
    }
    return;
}

sub set_gender_estimation {

## Function : Set gender estimation for sample_id
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $gender                => Gender to set
##          : $sample_id             => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $gender;
    my $sample_id;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        gender => {
            allow       => [qw{ female male other }],
            defined     => 1,
            required    => 1,
            store       => \$gender,
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

    $active_parameter_href->{gender_estimation}{$sample_id} = $gender;
    return;
}

sub set_gender_sample_ids {

## Function : Set gender sample_ids of the current analysis
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
        1      => q{males},
        2      => q{females},
        female => q{females},
        male   => q{males},
        other  => q{others},
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

            add_gender(
                {
                    active_parameter_href => $active_parameter_href,
                    sample_id             => $sample_id,
                    gender                => $gender_map{$sex},
                }
            );
            next SAMPLE_ID;
        }
        add_gender(
            {
                active_parameter_href => $active_parameter_href,
                sample_id             => $sample_id,
                gender                => q{others},
            }
        );
    }
    set_include_y(
        {
            active_parameter_href => $active_parameter_href,
        }
    );

    return;
}

sub set_include_y {

## Function : Set include_y bolean based on occurenes of male/other
## Returns  :
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

    $active_parameter_href->{include_y} = 0;

  GENDER:
    foreach my $gender ( keys %{ $active_parameter_href->{gender} } ) {

        next GENDER if ( $gender eq q{females} );

        $active_parameter_href->{include_y} =
          @{ $active_parameter_href->{gender}{$gender} } ? 1 : 0;

        last if $active_parameter_href->{include_y} == 1;
    }
    return;
}

sub set_default_container_config_file {

## Function : Set default install config file to active parameters
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

    ## Build default for mip install config
    my $path = catfile( $Bin, qw{ templates mip_container_config.yaml } );

    $active_parameter_href->{$parameter_name} = $path;
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

sub set_programs_for_installation {

## Function : Process the lists of programs that has been selected for installation
##          : and update the environment packages
## Returns  :
## Arguments: $active_parameter_href => The entire active parameter hash {REF}

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

    use Array::Utils qw{ array_minus };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Check that the options supplied are compatible with each other
    if (    @{ $active_parameter_href->{skip_programs} }
        and @{ $active_parameter_href->{select_programs} } )
    {
        $log->fatal(
            q{"--skip_programs" and "--select_programs" are mutually exclusive command line options}
        );
        exit 1;
    }

    ## Set all programs to install for pipelines
    if ( scalar @{ $active_parameter_href->{select_programs} } == 0 ) {

      PIPELINE:
        foreach my $pipeline ( @{ $active_parameter_href->{pipelines} } ) {

            push @{ $active_parameter_href->{select_programs} },
              @{ $active_parameter_href->{$pipeline} };
        }
    }

    ## Remove programs to be skipped from container
    delete @{ $active_parameter_href->{container} }{ @{ $active_parameter_href->{skip_programs} } };

    ## Remove programs to be skipped from select_programs
    @{ $active_parameter_href->{select_programs} } = array_minus(
        @{ $active_parameter_href->{select_programs} },
        @{ $active_parameter_href->{skip_programs} }
    );

    return if ( not @{ $active_parameter_href->{select_programs} } );

    ## Remove all programs that have not been selected for install
    my @container_programs = keys %{ $active_parameter_href->{container} };
    my @not_selected_programs =
      array_minus( @container_programs, @{ $active_parameter_href->{select_programs} } );
    delete @{ $active_parameter_href->{container} }{@not_selected_programs};

    return;
}

sub set_recipe_mode {

## Function : Set recipe mode
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $mode                  => Mode to set
##          : $recipes_ref           => Recipes to set mode for {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $mode;
    my $recipes_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        mode => {
            allow       => [ 0, 1, 2 ],
            defined     => 1,
            required    => 1,
            store       => \$mode,
            strict_type => 1,
        },
        recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$recipes_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $log = Log::Log4perl->get_logger($LOG_NAME);

  RECIPE:
    foreach my $recipe ( @{$recipes_ref} ) {

        $active_parameter_href->{$recipe} = $mode;

        $log->info( qq{Set $recipe to} . $COLON . $SPACE . $mode );
    }
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
        while ( my ( $recipe, $core_number ) = each %{ $active_parameter_href->{$set_hash_key} } ) {

            $active_parameter_href->{$target_hash_key}{$recipe} = $core_number;
        }
    }
    return;
}

sub set_vcfparser_outfile_counter {

## Function : Determine the number of outfile after vcfparser
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis

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

    ## Create link
    my %vcfparser_select_file = (
        sv_vcfparser  => { sv_vcfparser_select_file => q{sv_vcfparser_outfile_count} },
        vcfparser_ar  => { vcfparser_select_file    => q{vcfparser_outfile_count} },
        fusion_report => { fusion_select_file       => q{fusion_outfile_count} },
    );

## Determine if to expect select outfile for vcfparser and sv_vcfparser
  RECIPE:
    foreach my $recipe ( keys %vcfparser_select_file ) {

        next RECIPE if ( not $active_parameter_href->{$recipe} );

      FILES:
        while ( my ( $parameter_name, $parameter_name_counter ) =
            each %{ $vcfparser_select_file{$recipe} } )
        {

            $active_parameter_href->{$parameter_name_counter} = _set_vcfparser_file_counter(
                {
                    parameter_name => $active_parameter_href->{$parameter_name},
                }
            );
        }
    }
    return;
}

sub update_dynamic_config_parameters {

## Function : Updates a scalar reference with dynamic parameter values. Leaves other entries untouched.
## Returns  :
## Arguments: $active_parameter_ref   => Active parameter for this analysis {REF}
##          : $dynamic_parameter_href => Map of dynamic parameters

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_ref;
    my $dynamic_parameter_href;

    my $tmpl = {
        active_parameter_ref => {
            required => 1,
            store    => \$active_parameter_ref,
        },
        dynamic_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$dynamic_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  DYNAMIC_PARAMETER:
    while ( my ( $dynamic_parameter_name, $dynamic_parameter_value ) =
        each %{$dynamic_parameter_href} )
    {

        ${$active_parameter_ref} =~ s/$dynamic_parameter_name!/$dynamic_parameter_value/xsmgi;
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

sub update_recipe_mode_with_dry_run_all {

## Function : Update recipe mode depending on dry_run_all flag
## Returns  :
## Arguments: $active_parameter_href => The active parameters for this analysis hash {REF}
##          : $dry_run_all           => Simulation mode
##          : $recipes_ref           => Recipes in MIP

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $dry_run_all;
    my $recipes_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        dry_run_all => {
            allow       => [ undef, 0, 1, ],
            default     => 0,
            store       => \$dry_run_all,
            strict_type => 1,
        },
        recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$recipes_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( not $dry_run_all );

    ## Activate dry run mode
    my $dry_run_mode = 2;

  RECIPE:
    foreach my $recipe_name ( @{$recipes_ref} ) {

        ## If recipe is activated
        next RECIPE if ( not $active_parameter_href->{$recipe_name} );

        # Change recipe mode to simulation
        $active_parameter_href->{$recipe_name} = $dry_run_mode;
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
            foreach my $parameter_value ( @{ $active_parameter_href->{$parameter_name} } ) {

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

                $parameter_name_href->{$updated_key} = delete $parameter_name_href->{$key};
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

sub update_with_dynamic_config_parameters {

## Function : Updates the config file to particular user/cluster with dynamic config parameters following specifications. Leaves other entries untouched.
## Returns  :
## Arguments: $active_parameter_href  => Active parameters for this analysis hash {REF}
##          : $dynamic_parameter_href => Map of dynamic parameters

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $dynamic_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            required => 1,
            store    => \$active_parameter_href,
        },
        dynamic_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$dynamic_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ update_dynamic_config_parameters };

    ## Return if variable isn't in use
    return if ( not defined $active_parameter_href );

    ## HASH
    if ( ref $active_parameter_href eq q{HASH} ) {

      KEY:
        foreach my $key ( keys %{$active_parameter_href} ) {

            next KEY if ( not defined $active_parameter_href->{$key} );

            if (   ref $active_parameter_href->{$key} eq q{HASH}
                or ref $active_parameter_href->{$key} eq q{ARRAY} )
            {

                update_with_dynamic_config_parameters(
                    {
                        active_parameter_href  => $active_parameter_href->{$key},
                        dynamic_parameter_href => $dynamic_parameter_href,
                    }
                );
            }

            ## Update parameter
            update_dynamic_config_parameters(
                {
                    active_parameter_ref   => \$active_parameter_href->{$key},
                    dynamic_parameter_href => $dynamic_parameter_href,
                }
            );
        }
    }
    ##  ARRAY
    elsif ( ref $active_parameter_href eq q{ARRAY} ) {

      ELEMENT:
        while ( my ( $element_index, $element ) = each @{$active_parameter_href} ) {

            next ELEMENT if ( not defined $element );

            if ( ref $element eq q{HASH} or ref $element eq q{ARRAY} ) {

                update_with_dynamic_config_parameters(
                    {
                        active_parameter_href  => $element,
                        dynamic_parameter_href => $dynamic_parameter_href,
                    }
                );
            }

            ## Update parameter
            update_dynamic_config_parameters(
                {
                    active_parameter_ref   => \$active_parameter_href->[$element_index],
                    dynamic_parameter_href => $dynamic_parameter_href,
                }
            );
        }
    }
    return;
}

sub write_references {

## Function : Write references for this analysis to yaml
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $outfile_path          => Outfile path for reference yaml file
##          : $parameter_href        => Holds all parameters

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $outfile_path;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
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

    use MIP::Io::Write qw{ write_to_file };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };
    use MIP::Parameter qw{ get_parameter_attribute };

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => $LOG_NAME,
        }
    );
    my %reference;

  PARAMETER:
    while ( my ( $parameter_name, $parameter_value ) = each %{$active_parameter_href} ) {

        my $is_reference = get_parameter_attribute(
            {
                attribute      => q{is_reference},
                parameter_href => $parameter_href,
                parameter_name => $parameter_name,
            }
        );
        ## Only defined reference parameters
        if ($is_reference) {

            if ( ref $parameter_value eq q{HASH} ) {
                $reference{$parameter_name} = \%{$parameter_value};
                next PARAMETER;
            }
            if ( ref $parameter_value eq q{ARRAY} ) {
                $reference{$parameter_name} = \@{$parameter_value};
                next PARAMETER;
            }
            $reference{$parameter_name} = $parameter_value;
        }
    }

    # Writes hash to file
    write_to_file(
        {
            data_href => \%reference,
            format    => q{yaml},
            path      => $outfile_path,
        }
    );
    $log->info( q{Wrote reference YAML file to: } . $outfile_path );

    return;
}

sub _check_infile_directory {

## Function : Check if infile directory exists per sample id
## Returns  :
## Arguments: $infile_directory => Infile directory
##          : $sample_id        => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_directory;
    my $sample_id;

    my $tmpl = {
        infile_directory => {
            required    => 1,
            store       => \$infile_directory,
            strict_type => 1,
        },
        sample_id => {
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Log::MIP_log4perl qw{ retrieve_log };

    my $log = retrieve_log( { log_name => $LOG_NAME, } );

    return if ( defined $infile_directory );

    $log->fatal( q{Could not detect any supplied '--infile_dirs' for sample: } . $sample_id );
    exit 1;
}

sub _set_vcfparser_file_counter {

## Function : Return the expected number of outputfile(s) after vcfparser
## Returns  : 1 | 2
## Arguments: $parameter_name => Vcfparser select file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_name;

    my $tmpl =
      { parameter_name => { required => 1, store => \$parameter_name, strict_type => 1, }, };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## To track if vcfparser was used with a vcfparser_select_file (=2) or not (=1)
    # No select file was given
    return 1 if ( not defined $parameter_name );

    ## Select file was given
    return 2;
}

1;
