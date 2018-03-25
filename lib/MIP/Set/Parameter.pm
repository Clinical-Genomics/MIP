package MIP::Set::Parameter;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ fileparse };
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
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      set_config_to_active_parameters
      set_custom_default_to_active_parameter
      set_default_config_dynamic_parameters
      set_dynamic_parameter
      set_human_genome_reference_features
      set_parameter_reference_dir_path
      set_parameter_to_broadcast
    };
}

## Constants
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $TAB        => qq{\t};
Readonly my $UNDERSCORE => q{_};

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
                  $parameter_href->{supported_capture_kit}{default},
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
    ## Set rtg vcfeval reference genome
    if ( $parameter_name eq q{rtg_vcfeval_reference_genome} ) {

        ## Now we now what human genome reference to build from
        $active_parameter_href->{$parameter_name} =
          $active_parameter_href->{human_genome_reference};

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

sub set_human_genome_reference_features {

## Function : Detect version and source of the human_genome_reference: Source (hg19 or GRCh).
##            Used to change capture kit genome reference version later
## Returns  :
##          : $file_info_href         => File info hash {REF}
##          : $human_genome_reference => The human genome
##          : $log                    => Log

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $human_genome_reference;
    my $log;

    my $tmpl = {
        file_info_href => {
            default     => {},
            strict_type => 1,
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
        },
        human_genome_reference => {
            defined     => 1,
            required    => 1,
            store       => \$human_genome_reference,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Parameter qw{ check_gzipped };

    ## Different regexes for the two sources.
    ## i.e. Don't allow subversion of Refseq genome
    my %genome_source = (
        GRCh => qr/GRCh(\d+\.\d+|\d+)/,
        hg   => qr/hg(\d+)/,
    );

  GENOME_PREFIX:
    foreach my $genome_prefix ( keys %genome_source ) {

        ## Capture version
        my ($genome_version) = $human_genome_reference =~
          m/ $genome_source{$genome_prefix}_homo_sapiens /xms;

        if ($genome_version) {

            $file_info_href->{human_genome_reference_version} = $genome_version;
            $file_info_href->{human_genome_reference_source}  = $genome_prefix;
            last;
        }
    }
    if ( not $file_info_href->{human_genome_reference_version} ) {

        $log->fatal(
q{MIP cannot detect what version of human_genome_reference you have supplied.}
              . $SPACE
              . q{Please supply the reference on this format: [sourceversion]_[species] e.g. 'GRCh37_homo_sapiens' or 'hg19_homo_sapiens'}
              . $NEWLINE );
        exit 1;
    }

    ## Removes ".file_ending" in filename.FILENDING(.gz)
    $file_info_href->{human_genome_reference_name_prefix} =
      fileparse( $human_genome_reference, qr/\.fasta|\.fasta\.gz/ );

    $file_info_href->{human_genome_compressed} =
      check_gzipped( { file_name => $human_genome_reference, } );

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
