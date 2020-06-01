package MIP::Analysis;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## Third party module(s)
use autodie;
use List::MoreUtils qw{ all any };
use Log::Log4perl;
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $COMMA $EMPTY_STR $CLOSE_BRACE $CLOSE_BRACKET $LOG_NAME $NEWLINE $OPEN_BRACE $OPEN_BRACKET $SINGLE_QUOTE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.24;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      broadcast_parameters
      check_analysis_type_to_pipeline
      check_prioritize_variant_callers
      get_overall_analysis_type
      get_vcf_parser_analysis_suffix
      parse_prioritize_variant_callers
      set_parameter_to_broadcast
      update_prioritize_flag
      update_recipe_mode_for_wes
    };
}

sub broadcast_parameters {

## Function : Broadcast parameters to user
## Returns  : 1
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $broadcasts_ref        => Holds the parameters info for broadcasting later {REF}
##          : $order_parameters_ref  => Order of parameters (for structured output) {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $broadcasts_ref;
    my $order_parameters_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        broadcasts_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$broadcasts_ref,
            strict_type => 1,
        },
        order_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_parameters_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    set_parameter_to_broadcast(
        {
            active_parameter_href => $active_parameter_href,
            broadcasts_ref        => $broadcasts_ref,
            order_parameters_ref  => $order_parameters_ref,
        }
    );

    ## Broadcast set parameters info
    foreach my $parameter_info ( @{$broadcasts_ref} ) {

        $log->info($parameter_info);
    }
    return 1;
}

sub check_analysis_type_to_pipeline {

## Function : Check if consensus analysis type is compatible with the pipeline
## Returns  :
## Arguments: $analysis_type => Analysis type
##          : $pipeline      => Pipeline

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_type;
    my $pipeline;

    my $tmpl = {
        analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$analysis_type,
            strict_type => 1,
        },
        pipeline => {
            allow    => [qw{ dragen_rd_dna rd_dna rd_dna_panel rd_dna_vcf_rerun rd_rna }],
            defined  => 1,
            required => 1,
            store    => \$pipeline,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %analysis_pipeline_map = (
        dragen_rd_dna => q{dragen_rd_dna},
        mixed         => q{rd_dna},
        panel         => q{rd_dna_panel},
        vrn           => q{rd_dna_vcf_rerun},
        wes           => q{rd_dna},
        wgs           => q{rd_dna},
        wts           => q{rd_rna},
    );

    if ( $analysis_pipeline_map{$analysis_type} ne $pipeline ) {

        $log->fatal(
qq{Analysis type: $analysis_type is not compatible with MIP pipeline: $pipeline}
        );
        $log->fatal(
qq{Start MIP pipeline: $analysis_pipeline_map{$analysis_type} for this analysis type}
        );
        $log->fatal(q{Aborting run});
        exit 1;
    }
    return 1;
}

sub check_prioritize_variant_callers {

## Function : Check that all active variant callers have a prioritization order and that the prioritization elements
##          : match a supported variant caller
## Returns  :
## Arguments: $active_parameter_href      => Active parameters for this analysis hash {REF}
##          : $parameter_href             => Holds all parameters {REF}
##          : $priority_name_str          => Comma separated priority name str
##          : $variant_caller_recipes_ref => Variant caller recipes to check {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;
    my $priority_name_str;
    my $variant_caller_recipes_ref;

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
        priority_name_str => {
            defined     => 1,
            required    => 1,
            store       => \$priority_name_str,
            strict_type => 1,
        },
        variant_caller_recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$variant_caller_recipes_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Parameter qw{ get_parameter_attribute };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my @priority_order_names = split $COMMA, $priority_name_str;

    my %variant_caller_map;
  RECIPE_NAME:
    foreach my $recipe_name ( @{$variant_caller_recipes_ref} ) {

        $variant_caller_map{$recipe_name} = get_parameter_attribute(
            {
                attribute      => q{variant_caller},
                parameter_href => $parameter_href,
                parameter_name => $recipe_name,
            }
        );
    }

    ## Check that all active variant callers have a priority order
  RECIPE_NAME:
    while ( my ( $recipe_name, $variant_caller ) = each %variant_caller_map ) {

        ## Only active recipes
        if ( $active_parameter_href->{$recipe_name} ) {

            ## If variant caller alias is not part of priority order names
            if ( not any { $_ eq $variant_caller } @priority_order_names ) {

                $log->fatal( $priority_name_str
                      . q{ does not contain active variant caller: '}
                      . $variant_caller
                      . $SINGLE_QUOTE );
                exit 1;
            }
        }
        elsif ( any { $_ eq $variant_caller } @priority_order_names ) {
            ## If variant caller alias is part of priority order names

            $log->fatal( $priority_name_str
                  . q{ contains deactivated variant caller: '}
                  . $variant_caller
                  . $SINGLE_QUOTE );
            exit 1;
        }
    }

    ## Check that prioritize string contains valid variant call names
  PRIO_NAME:
    foreach my $priority_name (@priority_order_names) {

        # If priority order names is part of variant callers
        next PRIO_NAME
          if ( any { $_ eq $priority_name } values %variant_caller_map );

        my $variant_callers_str = join $COMMA, values %variant_caller_map;
        $log->fatal( $priority_name_str . q{: '}
              . $priority_name
              . q{' does not match any supported variant caller: '}
              . $variant_callers_str
              . $SINGLE_QUOTE );
        exit 1;
    }
    return 1;
}

sub get_overall_analysis_type {

## Function : Detect if all samples has the same analysis type and return consensus or mixed
## Returns  : q{consensus} | q{mixed} - analysis_type
## Arguments: $analysis_type_href => Analysis_type hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_type_href;

    my $tmpl = {
        analysis_type_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$analysis_type_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my @analysis_types = qw{ dragen_rd_dna panel vrn wes wgs wts };

  ANALYSIS:
    foreach my $analysis_type (@analysis_types) {

        ## If consensus is reached
        if ( all { $_ eq $analysis_type } values %{$analysis_type_href} ) {

            return $analysis_type;
        }
    }

    ## Check that the user supplied analysis type is supported
    foreach my $user_analysis_type ( values %{$analysis_type_href} ) {

        if ( not any { $_ eq $user_analysis_type } @analysis_types ) {

            $log->fatal(qq{' $user_analysis_type ' is not a supported analysis_type});
            $log->fatal( q{Supported analysis types are }
                  . $SINGLE_QUOTE
                  . join( q{', '}, @analysis_types )
                  . $SINGLE_QUOTE );
            $log->fatal(q{Aborting run});
            exit 1;
        }
    }

    # No consensus, then it must be mixed
    return q{mixed};
}

sub get_vcf_parser_analysis_suffix {

## Function : Get the vcf parser analysis suffix
## Returns  : @analysis_suffixes
## Arguments: $analysis_type           => Analysis type
##          : $vcfparser_outfile_count => Number of user supplied vcf parser outfiles

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_type;
    my $vcfparser_outfile_count;

    my $tmpl = {
        analysis_type => {
            store       => \$analysis_type,
            strict_type => 1,
        },
        vcfparser_outfile_count => {
            defined     => 1,
            required    => 1,
            store       => \$vcfparser_outfile_count,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    Readonly my $VCFPARSER_OUTFILE_COUNT => $vcfparser_outfile_count - 1;

    my $analysis_type_suffix = $EMPTY_STR;
    if ( defined $analysis_type and $analysis_type eq q{panel} ) {

        $analysis_type_suffix = q{all};
    }

    my @analysis_suffixes;

    ## Determined by vcfparser output
    # Set research (="") and selected file suffix
    for my $vcfparser_outfile_counter ( 0 .. $VCFPARSER_OUTFILE_COUNT ) {

        if ( $vcfparser_outfile_counter == 1 ) {

            ## Select file variants
            push @analysis_suffixes, q{selected};
            next;
        }
        push @analysis_suffixes, $analysis_type_suffix;
    }
    return @analysis_suffixes;
}

sub parse_prioritize_variant_callers {

## Function : Parse prioritization string of variant callers merging recipes
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

    use MIP::Analysis qw{ check_prioritize_variant_callers };
    use MIP::Parameter qw{ get_cache };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %priority_call_parameter = (
        variant_callers            => q{gatk_combinevariants_prioritize_caller},
        structural_variant_callers => q{sv_svdb_merge_prioritize},
    );

  PRIO_PARAMETER:
    while ( my ( $variant_caller_type, $prioritize_parameter_name ) =
        each %priority_call_parameter )
    {

        ## Skip if not part of pipeline
        next PRIO_PARAMETER if ( not $parameter_href->{$prioritize_parameter_name} );

        my @variant_caller_recipes = get_cache(
            {
                parameter_href => $parameter_href,
                parameter_name => $variant_caller_type,
            }
        );

        ## Check if we have any active recipes for callers
        my $has_active_recipe =
          grep { defined and $_ >= 1 } @{$active_parameter_href}{@variant_caller_recipes};

        if ($has_active_recipe) {

            check_prioritize_variant_callers(
                {
                    active_parameter_href => $active_parameter_href,
                    parameter_href        => $parameter_href,
                    priority_name_str =>
                      $active_parameter_href->{$prioritize_parameter_name},
                    variant_caller_recipes_ref => \@variant_caller_recipes,
                }
            );
            next PRIO_PARAMETER;
        }

        $log->warn(
qq{Could not find any active variant caller recipes for parameter: $prioritize_parameter_name}
        );
    }
    return 1;
}

sub set_parameter_to_broadcast {

## Function : Set parameters to broadcast message
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $broadcasts_ref        => Holds the parameters info for broadcasting later {REF}
##          : $order_parameters_ref  => Order of parameters (for structured output) {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $broadcasts_ref;
    my $order_parameters_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        broadcasts_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$broadcasts_ref,
            strict_type => 1,
        },
        order_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_parameters_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( not $active_parameter_href->{verbose} );

  PARAMETER:
    foreach my $parameter_name ( @{$order_parameters_ref} ) {

        next PARAMETER
          if ( not defined $active_parameter_href->{$parameter_name} );

        ## Hold parameters info
        my $info = q{Set } . $parameter_name . q{ to: };

        if ( ref $active_parameter_href->{$parameter_name} eq q{ARRAY} ) {

            $info = _parse_parameter_to_broadcast(
                {
                    info  => $info,
                    value => $active_parameter_href->{$parameter_name},
                }
            );

            ## Add info to broadcasts
            push @{$broadcasts_ref}, $info;
        }
        elsif ( ref $active_parameter_href->{$parameter_name} eq q{HASH} ) {

            $info = _parse_parameter_to_broadcast(
                {
                    info  => $info,
                    value => $active_parameter_href->{$parameter_name},
                }
            );

            ## Add info to broadcasts
            push @{$broadcasts_ref}, $info;
        }
        else {

            $info .= $active_parameter_href->{$parameter_name};

            ## Add info to broadcasts
            push @{$broadcasts_ref}, $info;
        }
    }
    return;
}

sub update_prioritize_flag {

## Function : Update prioritize flag depending on analysis run value as some recipes are not applicable for e.g. wes
## Returns  : $prioritize_key
## Arguments: $consensus_analysis_type => Consensus analysis_type
##          : $parameter_href          => Parameter hash {REF}
##          : $prioritize_key          => Prioritize key to update
##          : $recipes_ref             => Recipes to update {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $consensus_analysis_type;
    my $parameter_href;
    my $prioritize_key;
    my $recipes_ref;

    my $tmpl = {
        consensus_analysis_type => {
            defined     => 1,
            required    => 1,
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
        prioritize_key => {
            store       => \$prioritize_key,
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

    use MIP::Parameter qw{ get_parameter_attribute };

    return if ( not $prioritize_key );

    return $prioritize_key if ( $consensus_analysis_type eq q{wgs} );

    ## Split string into array
    my @variant_callers = split $COMMA, $prioritize_key;

  RECIPE:
    foreach my $recipe_name ( @{$recipes_ref} ) {

        my $exclude_variant_caller = get_parameter_attribute(
            {
                attribute      => q{variant_caller},
                parameter_href => $parameter_href,
                parameter_name => $recipe_name,
            }
        );

        ## Keep variant_callers which are not in recipes supplied
        @variant_callers = grep { not $_ eq $exclude_variant_caller } @variant_callers;
    }

    ## Update prioritize parameter
    $prioritize_key = join $COMMA, @variant_callers;
    return $prioritize_key;
}

sub update_recipe_mode_for_wes {

##Function : Update recipe mode depending on analysis run value as some recipes are not applicable for e.g. wes
##Returns  :
##Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##         : $consensus_analysis_type => Consensus analysis_type
##         : $recipes_ref             => Recipes to update {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $consensus_analysis_type;
    my $recipes_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        consensus_analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$consensus_analysis_type,
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

    return if ( $consensus_analysis_type eq q{wgs} );

    my @warning_msgs;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

  RECIPE:
    foreach my $recipe ( @{$recipes_ref} ) {

        ## Update recipe mode
        $active_parameter_href->{$recipe} = 0;

        my $warning_msg =
            q{Turned off: }
          . $recipe
          . q{ as it is not applicable for }
          . $consensus_analysis_type
          . q{ analysis}
          . $NEWLINE;

        push @warning_msgs, $warning_msg;
    }

    ## Broadcast
    if (@warning_msgs) {

      WARNING_MSG:
        foreach my $warning_msg (@warning_msgs) {
            $log->warn($warning_msg);
        }
    }
    return @warning_msgs;
}

sub _parse_parameter_to_broadcast {

## Function : Parse parameter to broadcast
## Returns  : $info
## Arguments: $info  => String to broadcast
##          : $value => Value to parse

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $info;
    my $value;

    my $tmpl = {
        info => {
            defined     => 1,
            required    => 1,
            store       => \$info,
            strict_type => 1,
        },
        value => {
            defined  => 1,
            required => 1,
            store    => \$value,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## HASH
    if ( ref $value eq q{HASH} ) {

        ## Start of hash
        $info .= $OPEN_BRACE;

      KEY:
        foreach my $key ( keys %{$value} ) {

            ## Key-value pairs
            $info .= $key . q{ => };

            if ( ref $value->{$key} eq q{HASH} ) {

                $info = _parse_parameter_to_broadcast(
                    {
                        info  => $info,
                        value => $value->{$key},
                    }
                );
                $info .= $COMMA . $SPACE;
                next KEY;
            }
            if ( ref $value->{$key} eq q{ARRAY} ) {

                $info .= $OPEN_BRACKET;

              ELEMENT:
                foreach my $element ( @{ $value->{$key} } ) {

                    $info = _parse_parameter_to_broadcast(
                        {
                            info  => $info,
                            value => $element,
                        }
                    );
                }
                ## Close array
                $info .= $CLOSE_BRACKET . $COMMA . $SPACE;
                next KEY;
            }
            if ( $value->{$key} ) {

                ## Scalar
                $info .= $value->{$key} . $COMMA . $SPACE;
            }
        }
        ## Close hash
        $info .= $CLOSE_BRACE;
        return $info;
    }
    ## ARRAY
    if ( ref $value eq q{ARRAY} ) {

        ## Open array
        $info .= $OPEN_BRACKET;

      ELEMENT:
        foreach my $element ( @{$value} ) {

            if ( ref $element eq q{HASH} ) {

                $info = _parse_parameter_to_broadcast(
                    {
                        info  => $info,
                        value => $element,
                    }
                );
                $info .= $COMMA . $SPACE;
                next ELEMENT;
            }
            if ( ref $element eq q{ARRAY} ) {

                $info .= $OPEN_BRACKET;

                foreach my $elements_ref ( @{$element} ) {

                    $info = _parse_parameter_to_broadcast(
                        {
                            info  => $info,
                            value => $elements_ref,
                        }
                    );
                }
                ## Close array
                $info .= $CLOSE_BRACKET . $COMMA . $SPACE;
                next ELEMENT;
            }
            if ($element) {

                ## Scalar
                $info .= $element . $COMMA . $SPACE;
            }
        }
        $info .= $CLOSE_BRACKET . $SPACE;
        return $info;
    }

    ## Scalar
    $info .= $value . $COMMA . $SPACE;
    return $info;
}

1;
