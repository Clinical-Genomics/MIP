package MIP::Set::Parameter;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ set_config_to_active_parameters set_custom_default_to_active_parameter set_default_config_dynamic_parameters set_dynamic_parameter set_parameter_reference_dir_path set_parameter_to_broadcast };
}

## Constants
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };
Readonly my $TAB     => qq{\t};

sub set_config_to_active_parameters {

## Function : Add contig parameters to active_parameters if not already initilized from command line
## Returns  :
## Arguments: $active_parameter_href  => Active parameters for this analysis hash {REF}
##          : $config_parameter_href => Config parameters hash

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $config_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        config_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$config_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  PARAMETER:
    foreach my $parmeter_name ( keys %{$config_parameter_href} ) {

        ## Cmd initilized HASH
        next PARAMETER
          if ( ref $active_parameter_href->{$parmeter_name} eq qw{HASH}
            && keys %{ $active_parameter_href->{$parmeter_name} } );

        ## Cmd initilized ARRAY
        next PARAMETER
          if ( ref $active_parameter_href->{$parmeter_name} eq qw{ARRAY}
            && @{ $active_parameter_href->{$parmeter_name} } );

        ## Cmd initilized scalar
        next PARAMETER
          if ( defined $active_parameter_href->{$parmeter_name}
            and not ref $active_parameter_href->{$parmeter_name} );

        ### No input from cmd
        ## Add to active_parameter
        $active_parameter_href->{$parmeter_name} =
          $config_parameter_href->{$parmeter_name};
    }
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
        parameter_name =>
          { defined => 1, required => 1, store => \$parameter_name, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::Parameter qw{ get_capture_kit };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Build default for analysis_type
    if ( $parameter_name eq q{analysis_type} ) {

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            $active_parameter_href->{$parameter_name}{$sample_id} = q{wgs};
        }
        return;
    }

    ## Set bwa build reference
    if ( $parameter_name eq q{bwa_build_reference} ) {

        ## Now we now what human genome reference to build from
        $active_parameter_href->{$parameter_name} =
          $active_parameter_href->{human_genome_reference};

        return;
    }

    ## If capture kit is not set after cmd, config and reading pedigree
    if ( $parameter_name eq q{exome_target_bed} ) {

        ## Return a default capture kit as user supplied no info
        my $capture_kit = get_capture_kit(
            {
                capture_kit => q{latest},
                supported_capture_kit_href =>
                  $parameter_href->{supported_capture_kit},
            }
        );

        ## Set default
        $active_parameter_href->{exome_target_bed}
          {$capture_kit} = join q{,},
          @{ $active_parameter_href->{sample_ids} };

        $log->warn(
q{Could not detect a supplied capture kit. Will Try to use 'latest' capture kit: }
              . $capture_kit );
        return;
    }

    ## Build default for infile_dirs
    if ( $parameter_name eq q{infile_dirs} ) {

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $path = catfile(
                $active_parameter_href->{cluster_constant_path},
                $active_parameter_href->{family_id},
                $active_parameter_href->{analysis_type}{$sample_id},
                $sample_id,
                q{fastq}
            );

            $active_parameter_href->{$parameter_name}{$path} = $sample_id;
        }
        return;
    }
    if ( $parameter_name eq q{sample_info_file} ) {

        $parameter_href->{sample_info_file}{default} = catfile(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
            $active_parameter_href->{family_id}
              . $UNDERSCORE
              . q{qc_sample_info.yaml}
        );

        $parameter_href->{qccollect_sampleinfo_file}{default} =
          $parameter_href->{sample_info_file}{default};
        return;
    }
    return;
}

sub set_default_config_dynamic_parameters {

## Function : Set default for config dynamic parameter using default definitions
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_href        => Parameter hash {REF}
##          : $parameter_names_ref   => MIP activate parameter names {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;
    my $parameter_names_ref;

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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  PARAMETER:
    foreach my $parameter_name ( @{$parameter_names_ref} ) {

        if ( exists $parameter_href->{$parameter_name}{default}
            and not defined $active_parameter_href->{$parameter_name} )
        {

            ## Transfer to active parameter
            $active_parameter_href->{$parameter_name} =
              $parameter_href->{$parameter_name}{default};
        }
    }
    return;
}

sub set_dynamic_parameter {

## Function : Sets dynamic aggregate information from definitions to parameter hash
## Returns  :
## Arguments: $aggregates_ref => The data to aggregate and add to parameter hash{REF}
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
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Constants
    Readonly my $RECORD_SEPARATOR => q{:};
    Readonly my $FIELD_COUNTER    => 2;

  PARAMETER:
    foreach my $parameter_name ( keys %{$parameter_href} ) {

      KEY_AND_STRING_TO_MATCH:
        foreach my $aggregate_element ( @{$aggregates_ref} ) {

            ## Split into key and string to match
            my ( $second_key, $string_to_match, $unexpected_data ) =
              split $RECORD_SEPARATOR, $aggregate_element, $FIELD_COUNTER + 1;

            ## Make sure that we get what we expect
            if ( defined $unexpected_data ) {

                carp
                  q{Unexpected trailing garbage at end of aggregate_element '}
                  . $aggregate_element
                  . q{':}, $NEWLINE . $TAB . $unexpected_data . $NEWLINE;
            }

            if ( defined $parameter_href->{$parameter_name}{$second_key}
                && $parameter_href->{$parameter_name}{$second_key} eq
                $string_to_match )
            {

                push @{ $parameter_href->{dynamic_parameter}{$string_to_match}
                  },
                  $parameter_name;
            }
        }
    }
    return;
}

sub set_parameter_reference_dir_path {

## Function : Set path for supplied reference(s) associated with parameter that should reside in the mip reference directory to full path.
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
            $active_parameter_href->{$parameter_name}
              { catfile( $reference_dir, $file_name ) } =
              delete $active_parameter_href->{$parameter_name}{$file};
        }
        return;
    }
    else {

        ## File
        ## Split to restate
        my ( $volume, $directory, $file_name ) =
          splitpath( $active_parameter_href->{$parameter_name} );

        ## Restate to allow for changing mip reference directory between runs
        $active_parameter_href->{$parameter_name} = $file_name;

        ## Update original value
        my $path =
          catfile( $reference_dir, $active_parameter_href->{$parameter_name} );
        $active_parameter_href->{$parameter_name} = $path;

    }
    return;
}

sub set_parameter_to_broadcast {

## Function : Set parameters to broadcast message
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $broadcasts_ref        => Holds the parameters info for broadcasting later {REF}
##          : $order_parameters_ref  => Order of parameters (for structured output) {REF}
##          : $parameter_href        => Holds all parameters

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $broadcasts_ref;
    my $order_parameters_ref;
    my $parameter_href;

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
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  PARAMETER:
    foreach my $parameter_name ( @{$order_parameters_ref} ) {

        next PARAMETER
          if ( not defined $active_parameter_href->{$parameter_name} );

        ## Hold parameters info
        my $info = q{Set } . $parameter_name . q{ to: };

        if ( ref $active_parameter_href->{$parameter_name} eq q{ARRAY} ) {

            ## Alias
            my $element_separator =
              $parameter_href->{$parameter_name}{element_separator};

            $info .= join $element_separator,
              @{ $active_parameter_href->{$parameter_name} };

            ## Add info to broadcasts
            push @{$broadcasts_ref}, $info;
        }
        elsif ( ref $active_parameter_href->{$parameter_name} eq q{HASH} ) {

            $info .= join q{,},
              map { qq{$_=$active_parameter_href->{$parameter_name}{$_}} }
              ( keys %{ $active_parameter_href->{$parameter_name} } );

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

1;
