package MIP::Parameter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw { any uniq };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COMMA $COLON $LOG_NAME $NEWLINE $SINGLE_QUOTE $SPACE $TAB $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_parameter_hash
      check_recipe_vs_binary_name
      get_cache
      get_capture_kit
      get_order_of_parameters
      get_parameter_attribute
      get_program_executables
      get_recipe_attributes
      parse_reference_path
      parse_parameter_files
      parse_parameter_recipe_names
      print_recipe
      set_cache
      set_cache_program_executables
      set_cache_sample_id_parameter
      set_custom_default_to_active_parameter
      set_default
      set_default_to_active_parameter
      set_parameter_build_file_status
    };
}

sub check_parameter_hash {

## Function : Evaluate parameters in parameters hash
## Returns  :
## Arguments: $file_path         => Path to yaml file
##          : $not_required_href => Hash with non required key {REF}
##          : $parameter_href    => Hash with parameters from yaml file {REF}
##          : $required_href     => Hash with required key {REF}

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $file_path;
    my $not_required_href;
    my $parameter_href;
    my $required_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        not_required_href => {
            default     => {},
            required    => 1,
            store       => \$not_required_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        required_href => {
            default     => {},
            required    => 1,
            store       => \$required_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Check that required keys exists for each parameter
    _check_parameter_required_keys_exits(
        {
            file_path      => $file_path,
            parameter_href => $parameter_href,
            required_href  => $required_href,
        }
    );

    ## Test parameter for both required and not_required keys data type and values
    my @keys = ( \%{$required_href}, \%{$not_required_href} );

  KEY_HASH_REF:
    foreach my $key_href (@keys) {

        _check_parameter_keys(
            {
                file_path      => $file_path,
                key_href       => $key_href,
                parameter_href => $parameter_href,
            }
        );
    }
    return 1;
}

sub check_recipe_vs_binary_name {

## Function : Check that recipe name and program binaries are not identical
## Returns  : 1
## Arguments: $parameter_href   => Parameter hash {REF}
##          : $recipe_names_ref => Recipe names {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $recipe_names_ref;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        recipe_names_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$recipe_names_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Get all program executables
    my @program_executables = get_cache(
        {
            parameter_href => $parameter_href,
            parameter_name => q{program_executables},
        }
    );

    ## Create hash map for check
    my %binary_name;
    @binary_name{@program_executables} = undef;

  RECIPE:
    foreach my $recipe ( @{$recipe_names_ref} ) {

        next RECIPE if ( not exists $binary_name{$recipe} );

        my $err_msg =
qq{Identical names for recipe and program: $recipe. Recipes cannot be identical to program binaries };
        croak($err_msg);
    }
    return 1;
}

sub get_cache {

## Function : Get aggregate information from parameter cache
## Returns  : $parameter_value | @{$parameter_value} | %{$parameter_value}
## Arguments: $parameter_href => Parameter hash {REF}
##          : $parameter_name => Parameter name to get chache for

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $parameter_name;

    my $tmpl = {
        parameter_href => {
            default  => {},
            defined  => 1,
            required => 1,
            store    => \$parameter_href,
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
    my $parameter_value = $parameter_href->{cache}{$parameter_name};

    return if ( not defined $parameter_value );

    return @{$parameter_value} if ( ref $parameter_value eq q{ARRAY} );

    return %{$parameter_value} if ( ref $parameter_value eq q{HASH} );

    ## Return scalar parameter cache value
    return $parameter_value;
}

sub get_capture_kit {

## Function : Return a capture kit depending on user info. If $is_set_by_user is
##          : true, go a head and add capture kit no matter what the switch was.
## Returns  : $capture kit, "supported capture kit" or "undef"
## Arguments: $capture_kit                => Capture kit to add
##          : $is_set_by_user             => Has user supplied parameter {OPTIONAL}
##          : $supported_capture_kit_href => Supported capture kits hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $capture_kit;
    my $is_set_by_user;
    my $supported_capture_kit_href;

    my $tmpl = {
        capture_kit                => { store => \$capture_kit,    strict_type => 1, },
        is_set_by_user             => { store => \$is_set_by_user, strict_type => 1, },
        supported_capture_kit_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$supported_capture_kit_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Unpack
    my $supported_capture_kit = $supported_capture_kit_href->{$capture_kit};

    ## Set default or return supplied capture kit
    if ( not defined $is_set_by_user ) {

        return $supported_capture_kit if ( defined $supported_capture_kit );

        ## Return unchanged capture_kit string
        return $capture_kit;
    }
    ## Only add if user supplied no info on parameter
    if ( defined $is_set_by_user
        and not $is_set_by_user )
    {

        return $supported_capture_kit if ( defined $supported_capture_kit );

        ## Return unchanged capture_kit string
        return $capture_kit;
    }
    return;
}

sub get_order_of_parameters {

## Function : Get order of parameters as they appear in definition file(s)
## Returns  : @order_of_parameters
## Arguments: $define_parameters_files_ref => MIPs define parameters file(s)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $define_parameters_files_ref;

    my $tmpl = {
        define_parameters_files_ref => {
            default     => [],
            store       => \$define_parameters_files_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Definition qw{ get_first_level_keys_order_from_definition_file };

    my @order_of_parameters;

  DEFINITION_FILE:
    foreach my $define_parameters_file ( @{$define_parameters_files_ref} ) {

        push @order_of_parameters,
          get_first_level_keys_order_from_definition_file(
            {
                file_path => $define_parameters_file,
            }
          );
    }
    return @order_of_parameters;
}

sub get_parameter_attribute {

## Function : Get parameter attribute
## Returns  : $attribute
## Arguments: $attribute      => Attribute to return
##          : $parameter_href => Parameter hash {REF}
##          : $parameter_name => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $attribute;
    my $parameter_href;
    my $parameter_name;

    my $tmpl = {
        attribute => {
            allow => [
                qw{ analysis_mode
                  associated_recipe
                  build_file
                  chain
                  data_type
                  default
                  exists_check
                  file_tag
                  is_reference
                  mandatory
                  outfile_suffix
                  program_executables
                  recipe_type
                  reference
                  type
                  variant_caller
                  }
            ],
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
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Return entire parameter attribute hash if no specific attribute key supplied
    return %{ $parameter_href->{$parameter_name} } if ( not defined $attribute );

    ## Unpack
    my $parameter_attribute = $parameter_href->{$parameter_name}{$attribute};

    return if ( not defined $parameter_attribute );

    return @{$parameter_attribute} if ( ref $parameter_attribute eq q{ARRAY} );

    return %{$parameter_attribute} if ( ref $parameter_attribute eq q{HASH} );

    ## Return scalar parameter attribute value
    return $parameter_attribute;
}

sub get_program_executables {

## Function : Get the parameter file program executables per recipe
## Returns  : uniq @program_executables
## Arguments: $parameter_href => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Get all program executables
    my @program_executables;

    my $err_msg = q{No keys 'cache' and 'recipe' in parameter hash};
    croak($err_msg) if ( not exists $parameter_href->{cache}{recipe} );

  RECIPE:
    foreach my $recipe ( @{ $parameter_href->{cache}{recipe} } ) {

        next RECIPE if ( not exists $parameter_href->{$recipe}{program_executables} );

        push @program_executables, @{ $parameter_href->{$recipe}{program_executables} };
    }

    ## Make unique and return
    return uniq(@program_executables);
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
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$recipe_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( not exists $parameter_href->{$recipe_name} ) {
        croak(qq{Recipe name: $recipe_name. Does not exists in parameter hash});
    }

    ## Get attribute value
    return $parameter_href->{$recipe_name}{$attribute} if ( defined $attribute and $attribute );

    ## Get recipe attribute hash
    return %{ $parameter_href->{$recipe_name} };
}

sub parse_parameter_files {

## Function : Parse parameter file objects and checks that their paths exist
## Returns  :
## Arguments: $active_parameter_href   => Holds all set parameter for analysis {REF}
##          : $consensus_analysis_type => Consensus analysis type for checking e.g. WGS specific files
##          : $parameter_href          => Holds all parameters {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $consensus_analysis_type;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        consensus_analysis_type => {
            store       => \$consensus_analysis_type,
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

    use MIP::Active_parameter qw{ check_parameter_files };

  PARAMETER:
    foreach my $parameter_name ( keys %{$parameter_href} ) {

        ## Unpack
        my %attribute = get_parameter_attribute(
            {
                parameter_href => $parameter_href,
                parameter_name => $parameter_name,
            }
        );

        next PARAMETER if ( not exists $attribute{exists_check} );

        check_parameter_files(
            {
                active_parameter_href   => $active_parameter_href,
                associated_recipes_ref  => $attribute{associated_recipe},
                build_status            => $attribute{build_file},
                consensus_analysis_type => $consensus_analysis_type,
                parameter_exists_check  => $attribute{exists_check},
                parameter_name          => $parameter_name,
            }
        );
    }
    return 1;
}

sub parse_parameter_recipe_names {

## Function : Parse parameter recipe names
## Returns  : 1
## Arguments: $parameter_href => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Recipes::Check qw{ check_recipe_exists_in_hash };

    ## Parameters with key(s) that have elements as MIP recipe names
    my @parameter_element_to_check = qw{ associated_recipe };

  PARAMETER:
    foreach my $parameter ( keys %{$parameter_href} ) {

      KEY:
        foreach my $parameter_name (@parameter_element_to_check) {

            next KEY if ( not exists $parameter_href->{$parameter}{$parameter_name} );

            ## Test if element from query array exists truth hash
            check_recipe_exists_in_hash(
                {
                    parameter_name => $parameter_name,
                    query_ref      => \@{ $parameter_href->{$parameter}{$parameter_name} },
                    truth_href     => $parameter_href,
                }
            );
        }
    }
    return 1;
}

sub parse_reference_path {

## Function : Parse reference parameteres and updates file name to path for parameters
##          : with key "reference"
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $parameter_href        => Holds all parameters {REF}

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

    use MIP::Active_parameter qw{ update_reference_parameters };

  PARAMETER:
    foreach my $parameter_name ( keys %{$parameter_href} ) {

        ## Expect file to be in reference directory
        next PARAMETER if ( not exists $parameter_href->{$parameter_name}{reference} );

        update_reference_parameters(
            {
                active_parameter_href  => $active_parameter_href,
                associated_recipes_ref =>
                  \@{ $parameter_href->{$parameter_name}{associated_recipe} },
                parameter_name => $parameter_name,
            }
        );
    }
    return;
}

sub print_recipe {

## Function : Print all supported recipes in '-prm' mode if requested and then exit
## Returns  :
## Arguments: $order_parameters_ref => Order of addition to parameter array {REF}
##          : $parameter_href       => Parameter hash {REF}
##          : $print_recipe         => Print recipes switch
##          : $print_recipe_mode    => Mode to run recipes in

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $order_parameters_ref;
    my $parameter_href;
    my $print_recipe;

    ## Default(s)
    my $print_recipe_mode;

    my $tmpl = {
        order_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_parameters_ref,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        print_recipe => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$print_recipe,
            strict_type => 1,
        },
        print_recipe_mode => {
            allow       => [ undef, 0, 1, 2 ],
            default     => $arg_href->{print_recipe_mode} //= 2,
            store       => \$print_recipe_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Definition qw{ get_first_level_keys_order_from_definition_file };
    use MIP::Parameter qw{ set_cache };

    ## Do not print
    return if ( not $print_recipe );

    set_cache(
        {
            aggregates_ref => [q{type:recipe}],
            parameter_href => $parameter_href,
        }
    );

  PARAMETER:
    foreach my $parameter ( @{$order_parameters_ref} ) {

        ## Only process recipes
        if (
            any { $_ eq $parameter }
            @{ $parameter_href->{cache}{recipe} }
          )
        {

            print {*STDOUT} q{--} . $parameter . $SPACE . $print_recipe_mode . $SPACE;

        }
    }
    print {*STDOUT} $NEWLINE;
    exit;
}

sub set_cache {

## Function : Set aggregate information from parameter hash to parameter cache
## Returns  :
## Arguments: $aggregates_ref => Data to aggregate {REF}
##          : $parameter_href => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $aggregates_ref;
    my $parameter_href;

    my $tmpl = {
        aggregates_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$aggregates_ref,
            strict_type => 1,
        },
        parameter_href => {
            default  => {},
            defined  => 1,
            required => 1,
            store    => \$parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Constants
    Readonly my $FIELD_COUNTER => 2;

  PARAMETER:
    foreach my $parameter_name ( keys %{$parameter_href} ) {

      KEY_AND_STRING_TO_MATCH:
        foreach my $aggregate_element ( @{$aggregates_ref} ) {

            ## Split into key and string to match
            my ( $second_key, $string_to_match, $unexpected_data ) =
              split $COLON, $aggregate_element, $FIELD_COUNTER + 1;

            ## Make sure that we get what we expect
            if ( defined $unexpected_data ) {

                carp q{Unexpected trailing garbage at end of aggregate_element '}
                  . $aggregate_element
                  . q{':}, $NEWLINE . $TAB . $unexpected_data . $NEWLINE;
            }

            next KEY_AND_STRING_TO_MATCH
              if ( not defined $parameter_href->{$parameter_name}{$second_key} );

            next KEY_AND_STRING_TO_MATCH
              if ( $parameter_href->{$parameter_name}{$second_key} ne $string_to_match );

            push @{ $parameter_href->{cache}{$string_to_match} }, $parameter_name;
        }
    }
    return;
}

sub set_cache_program_executables {

## Function : Set program executables per recipe to parameter cache
## Returns  :
## Arguments: $parameter_href => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    @{ $parameter_href->{cache}{program_executables} } =
      get_program_executables( { parameter_href => $parameter_href, } );

    return;
}

sub set_cache_sample_id_parameter {

## Function : Set parameter information to parameter cache at sample level
## Returns  :
## Arguments: $parameter_href  => Parameter hash {REF}
##          : $parameter_name  => Parameter to set
##          : $parameter_value => Parmeter value
##          : $sample_id       => Sample id to set parameter cache for

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $parameter_name;
    my $parameter_value;
    my $sample_id;

    my $tmpl = {
        parameter_href => {
            default  => {},
            defined  => 1,
            required => 1,
            store    => \$parameter_href,
        },
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
        parameter_value => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_value,
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

    $parameter_href->{cache}{$sample_id}{$parameter_name} = $parameter_value;
    return;
}

sub set_custom_default_to_active_parameter {

## Function : Checks and sets user input or default values to active_parameters.
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $parameter_href        => Holds all parameters {REF}
##          : $parameter_name        => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;
    my $parameter_name;

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
        parameter_name => { defined => 1, required => 1, store => \$parameter_name, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{
      set_default_analysis_type
      set_default_human_genome
      set_default_infile_dirs
      set_default_install_config_file
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
    };

    ## Set default value only to active_parameter
    my %set_default_parameter = (
        analysis_type => {
            method   => \&set_default_analysis_type,
            arg_href => {
                active_parameter_href => $active_parameter_href,
            },
        },
        bwa_build_reference => {
            method   => \&set_default_human_genome,
            arg_href => {
                active_parameter_href => $active_parameter_href,
                parameter_name        => $parameter_name,
            },
        },
        bwa_mem2_build_reference => {
            method   => \&set_default_human_genome,
            arg_href => {
                active_parameter_href => $active_parameter_href,
                parameter_name        => $parameter_name,
            },
        },
        exome_target_bed => {
            method   => \&_set_default_capture_kit,
            arg_href => {
                active_parameter_href => $active_parameter_href,
                parameter_href        => $parameter_href,
            },
        },
        fusion_select_file => {
            method   => \&set_default_vcfparser_select_file,
            arg_href => {
                active_parameter_href => $active_parameter_href,
                parameter_name        => $parameter_name,
            },
        },
        infile_dirs => {
            method   => \&set_default_infile_dirs,
            arg_href => {
                active_parameter_href => $active_parameter_href,
            },
        },
        install_config_file => {
            method   => \&set_default_install_config_file,
            arg_href => {
                active_parameter_href => $active_parameter_href,
                parameter_name        => $parameter_name,
            },
        },
        pedigree_fam_file => {
            method   => \&set_default_pedigree_fam_file,
            arg_href => {
                active_parameter_href => $active_parameter_href,
            },
        },
        program_test_file => {
            method   => \&set_default_program_test_file,
            arg_href => {
                active_parameter_href => $active_parameter_href,
            },
        },
        reference_dir => {
            method   => \&set_default_reference_dir,
            arg_href => {
                active_parameter_href => $active_parameter_href,
            },
        },
        reference_info_file => {
            method   => \&set_default_reference_info_file,
            arg_href => {
                active_parameter_href => $active_parameter_href,
            },
        },
        rtg_vcfeval_reference_genome => {
            method   => \&set_default_human_genome,
            arg_href => {
                active_parameter_href => $active_parameter_href,
                parameter_name        => $parameter_name,
            },
        },
        salmon_quant_reference_genome => {
            method   => \&set_default_human_genome,
            arg_href => {
                active_parameter_href => $active_parameter_href,
                parameter_name        => $parameter_name,
            },
        },
        sample_info_file => {
            method   => \&_set_default_sample_info_file,
            arg_href => {
                active_parameter_href => $active_parameter_href,
                parameter_href        => $parameter_href,
            },
        },
        select_programs => {
            method   => \&set_default_uninitialized_parameter,
            arg_href => {
                active_parameter_href => $active_parameter_href,
                parameter_name        => $parameter_name,
            },
        },
        skip_programs => {
            method   => \&set_default_uninitialized_parameter,
            arg_href => {
                active_parameter_href => $active_parameter_href,
                parameter_name        => $parameter_name,
            },
        },
        star_aln_reference_genome => {
            method   => \&set_default_human_genome,
            arg_href => {
                active_parameter_href => $active_parameter_href,
                parameter_name        => $parameter_name,
            },
        },
        star_fusion_reference_genome => {
            method   => \&set_default_human_genome,
            arg_href => {
                active_parameter_href => $active_parameter_href,
                parameter_name        => $parameter_name,
            },
        },
        store_file => {
            method   => \&set_default_store_file,
            arg_href => {
                active_parameter_href => $active_parameter_href,
            },
        },
        qccollect_store_metrics_outfile => {
            method   => \&set_default_qccollect_store_metrics_outfile,
            arg_href => {
                active_parameter_href => $active_parameter_href,
            },
        },
        sv_vcfparser_select_file => {
            method   => \&set_default_vcfparser_select_file,
            arg_href => {
                active_parameter_href => $active_parameter_href,
                parameter_name        => $parameter_name,
            },
        },
        temp_directory => {
            method   => \&set_default_temp_directory,
            arg_href => {
                active_parameter_href => $active_parameter_href,
            },
        },
        transcript_annotation_file_endings => {
            method   => \&set_default_transcript_annotation,
            arg_href => {
                active_parameter_href => $active_parameter_href,
                parameter_name        => $parameter_name,
            },
        },
        vcfparser_select_file => {
            method   => \&set_default_vcfparser_select_file,
            arg_href => {
                active_parameter_href => $active_parameter_href,
                parameter_name        => $parameter_name,
            },
        },
    );

    if ( exists $set_default_parameter{$parameter_name} ) {

        $set_default_parameter{$parameter_name}{method}
          ->( { %{ $set_default_parameter{$parameter_name}{arg_href} } } );

    }
    return;
}

sub set_default {

## Function : Set default from parameter hash to active_parameter for uninitilized parameters
## Returns  :
## Arguments: $active_parameter_href         => Active parameters for this analysis hash {REF}
##          : $custom_default_parameters_ref => Custom default parameters that are dependent on other parameters
##          : $parameter_href                => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $custom_default_parameters_ref;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        custom_default_parameters_ref => {
            default     => [],
            store       => \$custom_default_parameters_ref,
            strict_type => 1,
        },
        parameter_href => {
            default  => {},
            defined  => 1,
            required => 1,
            store    => \$parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ### Populate uninitilized active_parameters{parameter_name} with default from parameter
  PARAMETER:
    foreach my $parameter_name ( keys %{$parameter_href} ) {

        ## If hash and set - skip
        next PARAMETER
          if ( ref $active_parameter_href->{$parameter_name} eq qw{HASH}
            && keys %{ $active_parameter_href->{$parameter_name} } );

        ## If array and set - skip
        next PARAMETER
          if ( ref $active_parameter_href->{$parameter_name} eq qw{ARRAY}
            && @{ $active_parameter_href->{$parameter_name} } );

        ## If scalar and set - skip
        next PARAMETER
          if ( defined $active_parameter_href->{$parameter_name}
            and not ref $active_parameter_href->{$parameter_name} );

        if ( any { $_ eq $parameter_name } @{$custom_default_parameters_ref} ) {

            set_custom_default_to_active_parameter(
                {
                    active_parameter_href => $active_parameter_href,
                    parameter_href        => $parameter_href,
                    parameter_name        => $parameter_name,
                }
            );
            next PARAMETER;
        }

        ## Checks and sets user input or default values to active_parameters
        set_default_to_active_parameter(
            {
                active_parameter_href  => $active_parameter_href,
                associated_recipes_ref =>
                  \@{ $parameter_href->{$parameter_name}{associated_recipe} },
                parameter_href => $parameter_href,
                parameter_name => $parameter_name,
            }
        );
    }
    return 1;
}

sub set_default_to_active_parameter {

## Function : Checks and sets user input or default values to active_parameters.
## Returns  :
## Arguments: $active_parameter_href  => Holds all set parameter for analysis {REF}
##          : $associated_recipes_ref => The parameters recipe {REF}
##          : $parameter_href         => Holds all parameters {REF}
##          : $parameter_name         => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $associated_recipes_ref;
    my $parameter_href;
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
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        parameter_name => { defined => 1, required => 1, store => \$parameter_name, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ set_default_parameter };

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %only_wgs = ( gatk_genotypegvcfs_ref_gvcf => 1, );

    ## Alias
    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};

    ## Do nothing since parameter is not required unless exome mode is enabled
    return
      if ( exists $only_wgs{$parameter_name}
        && $consensus_analysis_type =~ / wgs /xsm );

    ## Check all recipes that use parameter
  ASSOCIATED_RECIPE:
    foreach my $associated_recipe ( @{$associated_recipes_ref} ) {

        ## Default exists
        if ( exists $parameter_href->{$parameter_name}{default} ) {

            ## Scalar or array or hash default
            set_default_parameter(
                {
                    active_parameter_href => $active_parameter_href,
                    parameter_name        => $parameter_name,
                    parameter_default     => $parameter_href->{$parameter_name}{default},
                }
            );

            ## Set default - no use in continuing
            return;
        }
        ## No parameter default

        ## Not mandatory paramater - skip
        return
          if ( exists $parameter_href->{$parameter_name}{mandatory}
            && $parameter_href->{$parameter_name}{mandatory} eq q{no} );

        next ASSOCIATED_RECIPE
          if ( not $active_parameter_href->{$associated_recipe} );

        ## Mandatory parameter not supplied
        $log->fatal(
            q{Supply '-} . $parameter_name . q{' if you want to run } . $associated_recipe );
        exit 1;
    }
    return;
}

sub set_parameter_build_file_status {

## Function : Set parameter build file status
## Returns  :
## Arguments: $parameter_href => Parameter hash {REF}
##          : $parameter_name => Parameter name
##          : $status         => Status

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $parameter_name;
    my $status;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
        status => {
            allow       => [ 0, 1 ],
            defined     => 1,
            required    => 1,
            store       => \$status,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    $parameter_href->{$parameter_name}{build_file} = $status;
    return;
}

sub _check_parameter_required_keys_exits {

## Function : Check that required keys exists
## Returns  :
## Arguments: $file_path      => Path to yaml file
##          : $parameter_href => Hash with parameters from yaml file {REF}
##          : $required_href  => Hash with required key {REF}

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $file_path;
    my $parameter_href;
    my $required_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        required_href => {
            default     => {},
            required    => 1,
            store       => \$required_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  PARAMETER:
    foreach my $parameter ( keys %{$parameter_href} ) {

      REQUIRED_KEY:
        foreach my $required_key ( keys %{$required_href} ) {

            next REQUIRED_KEY
              if ( exists $parameter_href->{$parameter}{$required_key} );

            say {*STDERR} q{Missing required key: '}
              . $required_key
              . q{' for parameter: '}
              . $parameter
              . q{' in file: '}
              . $file_path
              . $SINGLE_QUOTE
              . $NEWLINE;
            croak();
        }
    }
    return;
}

sub _check_parameter_keys {

## Function : Evaluate parameter keys in hash
## Returns  :
## Arguments: $file_path      => Path to yaml file
##          : $key_href       => Hash with keys to check for parameter {REF}
##          : $parameter_href => Hash with parameters from yaml file {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $key_href;
    my $parameter_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        key_href => {
            default     => {},
            required    => 1,
            store       => \$key_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  PARAMETER:
    foreach my $parameter ( keys %{$parameter_href} ) {

      KEY:
        foreach my $key ( keys %{$key_href} ) {

            next KEY if ( not exists $parameter_href->{$parameter}{$key} );

            ## Check key data type
            _check_parameter_data_type(
                {
                    file_path      => $file_path,
                    key            => $key,
                    key_href       => $key_href,
                    parameter      => $parameter,
                    parameter_href => $parameter_href,
                }
            );

            ## Check key values
            _check_parameter_values(
                {
                    file_path      => $file_path,
                    key            => $key,
                    key_href       => $key_href,
                    parameter      => $parameter,
                    parameter_href => $parameter_href,
                }
            );
        }
    }
    return;
}

sub _check_parameter_data_type {

## Function : Check key data type
## Returns  :
## Arguments: $file_path      => Path to yaml file
##          : $key            => Hash with non key
##          : $key_href       => Hash with key {REF}
##          : $parameter      => Parameter
##          : $parameter_href => Hash with parameters from yaml file {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $key;
    my $key_href;
    my $parameter;
    my $parameter_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        key      => { defined => 1, required => 1, store => \$key, strict_type => 1, },
        key_href => {
            default     => {},
            required    => 1,
            store       => \$key_href,
            strict_type => 1,
        },
        parameter => {
            defined     => 1,
            required    => 1,
            store       => \$parameter,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Get actual parameter value data type
    my $data_type = ref $parameter_href->{$parameter}{$key};

    ## Array or hash
    if ($data_type) {

        ## Return if parameter value data type equals description
        return if ( $data_type eq $key_href->{$key}{key_data_type} );

        say {*STDERR} q{Found '}
          . $data_type
          . q{' but expected datatype '}
          . $key_href->{$key}{key_data_type}
          . q{' for parameter: '}
          . $parameter
          . q{' in key: '}
          . $key
          . q{' in file: '}
          . $file_path
          . $SINGLE_QUOTE
          . $NEWLINE;
        croak();
    }
    elsif ( $key_href->{$key}{key_data_type} ne q{SCALAR} ) {

        ## Wrong data_type
        say {*STDERR} q{Found 'SCALAR' but expected datatype '}
          . $key_href->{$key}{key_data_type}
          . q{' for parameter: '}
          . $parameter
          . q{' in key: '}
          . $key
          . q{' in file: '}
          . $file_path
          . $SINGLE_QUOTE
          . $NEWLINE;
        croak();
    }
    return;
}

sub _check_parameter_values {

## Function : Evaluate parameter key values
## Returns  :
## Arguments: $file_path      => Path to yaml file
##          : $key            => Hash with non  key
##          : $key_href       => Hash with  key {REF}
##          : $parameter      => Parameter
##          : $parameter_href => Hash with parameters from yaml file {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $key;
    my $key_href;
    my $parameter;
    my $parameter_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        key      => { defined => 1, required => 1, store => \$key, strict_type => 1, },
        key_href => {
            default     => {},
            required    => 1,
            store       => \$key_href,
            strict_type => 1,
        },
        parameter => {
            defined     => 1,
            required    => 1,
            store       => \$parameter,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Check value(s)
    if ( $key_href->{$key}{values} ) {

        my $value = $parameter_href->{$parameter}{$key};

        if ( not( any { $_ eq $value } @{ $key_href->{$key}{values} } ) ) {

            say {*STDERR} q{Found illegal value '}
              . $value
              . q{' for parameter: '}
              . $parameter
              . q{' in key: '}
              . $key
              . q{' in file: '}
              . $file_path
              . $SINGLE_QUOTE
              . $NEWLINE
              . q{Allowed entries: '}
              . join( q{', '}, @{ $key_href->{$key}{values} } )
              . $SINGLE_QUOTE
              . $NEWLINE;
            croak();
        }
    }
    return;
}

sub _set_default_capture_kit {

## Function : Set default capture kit to active parameters
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $parameter_href        => Holds all parameters {REF}

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

    use MIP::Active_parameter qw{ set_exome_target_bed };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ### If capture kit is not set after cmd, config and reading pedigree
    ## Return a default capture kit as user supplied no info
    my $capture_kit = get_capture_kit(
        {
            capture_kit                => q{latest},
            supported_capture_kit_href => $parameter_href->{supported_capture_kit}{default},
        }
    );

    ## Set default
    my $sample_id_string = join $COMMA, @{ $active_parameter_href->{sample_ids} };

    set_exome_target_bed(
        {
            active_parameter_href => $active_parameter_href,
            exome_target_bed_file => $capture_kit,
            sample_id_string      => $sample_id_string,
        }
    );

    $log->warn(
        q{Could not detect a supplied capture kit. Will Try to use 'latest' capture kit: }
          . $capture_kit );
    return;
}

sub _set_default_sample_info_file {

## Function : Set default sample_info_file and qccollect_sampleinfo_file to parameters
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $parameter_href        => Holds all parameters {REF}

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

    ## Set sample info file
    $parameter_href->{sample_info_file}{default} = catfile(
        $active_parameter_href->{outdata_dir},
        $active_parameter_href->{case_id},
        $active_parameter_href->{case_id} . $UNDERSCORE . q{qc_sample_info.yaml}
    );

    ## Set qccollect sampleinfo file input
    $parameter_href->{qccollect_sampleinfo_file}{default} =
      $parameter_href->{sample_info_file}{default};
    return;
}

1;
