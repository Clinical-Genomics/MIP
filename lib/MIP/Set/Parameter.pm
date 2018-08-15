package MIP::Set::Parameter;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ fileparse };
use File::Spec::Functions qw{ catdir catfile splitpath };
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
    our $VERSION = 1.07;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      set_conda_env_names_and_paths
      set_config_to_active_parameters
      set_custom_default_to_active_parameter
      set_default_config_dynamic_parameters
      set_default_to_active_parameter
      set_dynamic_parameter
      set_human_genome_reference_features
      set_no_dry_run_parameters
      set_parameter_reference_dir_path
      set_parameter_to_broadcast
      set_programs_for_installation
    };
}

## Constants
Readonly my $MINUS_ONE  => -1;
Readonly my $MINUS_TWO  => -2;
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
    ## Set sample info file
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

    ## Set default path to expansionhunter repeat specs if needed
    if (    ( $parameter_name eq q{expansionhunter_repeat_specs_dir} )
        and ( not $active_parameter_href->{expansionhunter_repeat_specs_dir} ) )
    {

        $active_parameter_href->{expansionhunter_repeat_specs_dir} =
          _get_default_repeat_specs_dir_path(
            {
                reference_genome_path =>
                  $active_parameter_href->{human_genome_reference},
            }
          );
        return;
    }

    ## Set default dynamic path if needed
    my %dynamic_path = (
        gatk_path => {
            bin_file        => q{gatk3},
            environment_key => q{gatk},
        },
        picardtools_path => {
            bin_file        => q{picard.jar},
            environment_key => q{picardtools},
        },
        snpeff_path => {
            bin_file        => q{snpEff.jar},
            environment_key => q{snpeff},
        },
        vep_directory_path => {
            bin_file        => q{vep},
            environment_key => q{varianteffectpredictor},
        },
    );

    use MIP::Get::Parameter qw{ get_dynamic_conda_path };

    if (    ( defined $dynamic_path{$parameter_name} )
        and ( not $active_parameter_href->{$parameter_name} ) )
    {

        $active_parameter_href->{$parameter_name} = get_dynamic_conda_path(
            {
                active_parameter_href => $active_parameter_href,
                bin_file => $dynamic_path{$parameter_name}{bin_file},
                environment_key =>
                  $dynamic_path{$parameter_name}{environment_key},
            }
        );
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

sub set_default_to_active_parameter {

## Function : Checks and sets user input or default values to active_parameters.
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis
##          : $associated_programs   => The parameters program(s) {array, REF}
##          : $log                   => Log object
##          : $parameter_href        => Holds all parameters
##          : $parameter_name        => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $associated_programs_ref;
    my $log;
    my $parameter_href;
    my $parameter_name;

    ## Default(s)
    my $family_id;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        associated_programs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$associated_programs_ref,
            strict_type => 1,
        },
        log => {
            required => 1,
            defined  => 1,
            store    => \$log
        },
        parameter_name =>
          { defined => 1, required => 1, store => \$parameter_name, },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %only_wgs = ( gatk_genotypegvcfs_ref_gvcf => 1, );

    ## Alias
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};

    ## Do nothing since parameter is not required unless exome mode is enabled
    return
      if ( exists $only_wgs{$parameter_name}
        && $consensus_analysis_type =~ / wgs /xsm );

    ## Check all programs that use parameter
  ASSOCIATED_PROGRAM:
    foreach my $associated_program ( @{$associated_programs_ref} ) {

        ## Only add active programs parameters
        next ASSOCIATED_PROGRAM
          if ( not defined $active_parameter_href->{$associated_program} );

        next ASSOCIATED_PROGRAM
          if ( not $active_parameter_href->{$associated_program} );

        ## Default exists
        if ( exists $parameter_href->{$parameter_name}{default} ) {

            ## Array reference
            if ( $parameter_href->{$parameter_name}{data_type} eq q{ARRAY} ) {

                push
                  @{ $active_parameter_href->{$parameter_name} },
                  @{ $parameter_href->{$parameter_name}{default} };
            }
            elsif ( $parameter_href->{$parameter_name}{data_type} eq q{HASH} ) {
                ## Hash reference

                $active_parameter_href->{$parameter_name} =
                  $parameter_href->{$parameter_name}{default};
            }
            else {
                ## Scalar

                $active_parameter_href->{$parameter_name} =
                  $parameter_href->{$parameter_name}{default};
            }

            ## Set default - no use in continuing
            return;
        }
        else {
            ## No default

            ## Not mandatory - skip
            return
              if ( exists $parameter_href->{$parameter_name}{mandatory}
                && $parameter_href->{$parameter_name}{mandatory} eq q{no} );

            ## Mandatory parameter not supplied
            $log->fatal( q{Supply '-}
                  . $parameter_name
                  . q{' if you want to run }
                  . $associated_program );
            exit 1;
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

## Function : Detect version and source of the human_genome_reference: Source (hg19 or GRCh) as well as compression status.
##            Used to change capture kit genome reference version later
## Returns  :
##          : $file_info_href         => File info hash {REF}
##          : $human_genome_reference => The human genome
##          : $log                    => Log
##          : $parameter_href         => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $human_genome_reference;
    my $log;
    my $parameter_href;

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
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Parameter qw{ check_gzipped };

    ## Different regexes for the two sources.
    ## i.e. Don't allow subversion of Refseq genome
    my %genome_source = (
        GRCh => qr/GRCh(\d+[.]\d+ | \d+)/xsm,
        hg   => qr/hg(\d+)/xsm,
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
      fileparse( $human_genome_reference, qr/[.]fasta | [.]fasta[.]gz/xsm );

    $file_info_href->{human_genome_compressed} =
      check_gzipped( { file_name => $human_genome_reference, } );

    if ( $file_info_href->{human_genome_compressed} ) {

        ## Set build file to one to allow for uncompression before analysis
        $parameter_href->{human_genome_reference_file_endings}{build_file} = 1;
    }
    return;

}

sub set_no_dry_run_parameters {

## Function : Set parameters for true run i.e. not a dry run
## Returns  :
## Arguments: $analysis_date    => Analysis date
##          : $is_dry_run_all   => Dry run boolean
##          : $mip_version      => MIP version
##          : $sample_info_href => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_date;
    my $is_dry_run_all;
    my $mip_version;
    my $sample_info_href;

    my $tmpl = {
        analysis_date => {
            defined     => 1,
            required    => 1,
            store       => \$analysis_date,
            strict_type => 1,
        },
        is_dry_run_all => {
            allow       => [ 0, 1, undef ],
            required    => 1,
            store       => \$is_dry_run_all,
            strict_type => 1,
        },
        mip_version => {
            defined     => 1,
            required    => 1,
            store       => \$mip_version,
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

    return if ($is_dry_run_all);

    my %no_dry_run_info = (
        analysisrunstatus => q{not_finished},
        analysis_date     => $analysis_date,
        mip_version       => $mip_version,
    );

  KEY_VALUE_PAIR:
    while ( my ( $key, $value ) = each %no_dry_run_info ) {

        $sample_info_href->{$key} = $value;
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

          PARAMETER_KEY:
            while ( my ( $key, $value ) =
                each %{ $active_parameter_href->{$parameter_name} } )
            {

                if ( ref $value eq q{ARRAY} ) {

                    $info .= join q{,}, map {
                        qq{$_=} . join $SPACE,
                          @{ $active_parameter_href->{$parameter_name}{$_} }
                    } ( keys %{ $active_parameter_href->{$parameter_name} } );

                    last PARAMETER_KEY;
                }
                else {

                    $info .= join q{,}, map {
                        qq{$_=$active_parameter_href->{$parameter_name}{$_}}
                    } ( keys %{ $active_parameter_href->{$parameter_name} } );

                    last PARAMETER_KEY;
                }
            }

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

sub set_conda_env_names_and_paths {

## Function : Set conda environmnet specific names and paths
## Returns  :
## Arguments: $log            => Log
##          : $parameter_href => The entire parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $log;

    my $tmpl = {
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
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

    use List::Util qw{ first };

    my @environments = @{ $parameter_href->{installations} };

    ## Get array index for the emip environment
    my $emip_idx = first { $environments[$_] eq q{emip} } 0 .. $#environments;

    ## Set up conda prefix path for MIP main environment
    if ( defined $emip_idx ) {

        if ( $parameter_href->{environment_name}{emip} ) {
            $parameter_href->{emip}{conda_prefix_path} =
              catdir( $parameter_href->{conda_dir_path},
                q{envs}, $parameter_href->{environment_name}{emip} );
        }
        else {
            $log->warn(
q{No environment name has been specified for MIP's main environment.}
            );
            $log->warn(q{MIP will be installed in conda's base environment.});
            $parameter_href->{emip}{conda_prefix_path} =
              $parameter_href->{conda_dir_path};
        }

        ## Remove emip from environments array so that the emip conda path is not overwritten later
        splice @environments, $emip_idx, 1;
    }

    ## Set up conda environment names and prefix paths for non mip environmnents
    foreach my $environment (@environments) {

        ## Give the env a default name if not given
        if ( not $parameter_href->{environment_name}{$environment} ) {

            ## Add the env name to mip base name if it is named
            if ( $parameter_href->{environment_name}{emip} ) {

                $parameter_href->{environment_name}{$environment} =
                    $parameter_href->{environment_name}{emip}
                  . $UNDERSCORE
                  . $environment;
            }
            else {
                $parameter_href->{environment_name}{$environment} =
                  $environment;
            }
        }

        ## Add environment specific conda prefix path
        $parameter_href->{$environment}{conda_prefix_path} =
          catdir( $parameter_href->{conda_dir_path},
            q{envs}, $parameter_href->{environment_name}{$environment} );
    }
    return;
}

sub _get_default_repeat_specs_dir_path {

## Function : Return the path to the repeat specs directory in the Expansionhunter directory
## Returns  : $repeat_specs_dir_path
## Arguments: $reference_genome_path => Path to the reference genome used

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $reference_genome_path;

    my $tmpl = {
        reference_genome_path => {
            defined     => 1,
            required    => 1,
            store       => \$reference_genome_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Cwd qw{ abs_path };
    use File::Basename qw{ fileparse };
    use File::Find::Rule;
    use IPC::Cmd qw{ can_run };

    ## Path to set
    my $repeat_specs_dir_path;

    ## Get path to binary
    my $expansionhunter_bin_path = can_run(q{ExpansionHunter});

    ## Return if path not found,
    ## MIP requires a defined variable in order to flag that it can't find the dir
    if ( not $expansionhunter_bin_path ) {
        return q{Failed to find default path};
    }

    ## Follow potential link
    $expansionhunter_bin_path = abs_path($expansionhunter_bin_path);

    ## Get the path to the repeat specs dirs
    my @expansionhunter_dirs = File::Spec->splitdir($expansionhunter_bin_path);
    splice @expansionhunter_dirs, $MINUS_TWO;
    my $parent_repeat_specs_dir_path =
      catdir( @expansionhunter_dirs, qw{ data repeat-specs } );

    ## Get list of genome version directories
    my @repeat_specs_dir_paths =
      File::Find::Rule->directory->in($parent_repeat_specs_dir_path);

    ## Remove top directory
    @repeat_specs_dir_paths =
      grep { !/^$parent_repeat_specs_dir_path$/xms } @repeat_specs_dir_paths;

    ## Find correct repeat spec folder
    my $genome_reference = fileparse($reference_genome_path);
  REPEAT_SPECS_VERSION:
    foreach my $repeat_specs_version (@repeat_specs_dir_paths) {

        ## Get version
        my @genome_version_dirs = File::Spec->splitdir($repeat_specs_version);
        my $genome_version_dir = splice @genome_version_dirs, $MINUS_ONE;

        ## Match version to reference used
        if ( $genome_reference =~ / $genome_version_dir /ixms ) {
            $repeat_specs_dir_path = $repeat_specs_version;
            last;
        }
    }

    ## MIP requires a defined variable in order to flag that it can't find the dir
    if ( not -d $repeat_specs_dir_path ) {
        return q{Failed to find default path};
    }
    return $repeat_specs_dir_path;
}

sub set_programs_for_installation {

## Function : Proccess the lists of programs that has been selected for or omitted from installation
##          : and update the environment packages
## Returns  :
## Arguments: $installation   => Environment to be installed
##          : $log            => Log
##          : $parameter_href => The entire parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $installation;
    my $log;
    my $parameter_href;

    my $tmpl = {
        installation => {
            defined     => 1,
            required    => 1,
            store       => \$installation,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
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

    use Array::Utils qw{ array_minus };
    use Data::Diver qw{ Dive };
    use List::Util qw{ any };
    use MIP::Get::Parameter qw{ get_programs_for_shell_installation };
    use MIP::Check::Installation
      qw{ check_and_add_dependencies check_python_compability };

    ## Check that the options supplied are compatible with each other
    if (    ( scalar @{ $parameter_href->{skip_programs} } > 0 )
        and ( scalar @{ $parameter_href->{select_programs} } > 0 ) )
    {
        $log->fatal(
q{"--skip_programs" and "--select_programs" are mutually exclusive command line options}
        );
        exit 1;
    }

    ## Check that only one environment has been specified for installation if the option select_program is used
    if (    ( scalar @{ $parameter_href->{select_programs} } > 0 )
        and ( scalar @{ $parameter_href->{installations} } > 1 ) )
    {
        $log->fatal(
q{Please select a single installation environment when using the option --select_programs}
        );
        exit 1;
    }

    ## Get programs that are to be installed via shell
    my @shell_programs_to_install = get_programs_for_shell_installation(
        {
            conda_programs_href => $parameter_href->{$installation}{conda},
            log                 => $log,
            prefer_shell        => $parameter_href->{prefer_shell},
            shell_install_programs_ref => $parameter_href->{shell_install},
            shell_programs_href => $parameter_href->{$installation}{shell},
        }
    );

    ## Remove the conda packages that has been selected to be installed via SHELL
    delete @{ $parameter_href->{$installation}{conda} }
      {@shell_programs_to_install};

    ## Special case for snpsift since it is installed together with SnpEff
    ## if shell installation of SnpEff has been requested.
    if ( any { $_ eq q{snpeff} } @shell_programs_to_install ) {
        delete $parameter_href->{$installation}{conda}{snpsift};
    }
    ## Store variable outside of shell hash and use Data::Diver module to avoid autovivification of variable
    $parameter_href->{$installation}{snpeff_genome_versions} = Dive(
        $parameter_href->{$installation},
        qw{ shell snpeff snpeff_genome_versions }
    );

    ## Delete shell programs that are to be installed via conda instead of shell
    my @shell_programs_to_delete =
      keys %{ $parameter_href->{$installation}{shell} };
    @shell_programs_to_delete =
      array_minus( @shell_programs_to_delete, @shell_programs_to_install );
    delete @{ $parameter_href->{$installation}{shell} }
      {@shell_programs_to_delete};

    ## Solve the installation when the skip_program or select_program parameter has been used
  INSTALL_MODE:
    foreach my $install_mode (qw{ conda pip shell }) {

        ## Remove programs that are to be skipped
        delete @{ $parameter_href->{$installation}{$install_mode} }
          { @{ $parameter_href->{skip_programs} } };

        ## Remove all non-selected programs
        if ( scalar @{ $parameter_href->{select_programs} } > 0 ) {
            my @non_selects =
              keys %{ $parameter_href->{$installation}{$install_mode} };
            @non_selects = array_minus( @non_selects,
                @{ $parameter_href->{select_programs} } );
            delete @{ $parameter_href->{$installation}{$install_mode} }
              {@non_selects};
        }
    }

    use Data::Printer;
    ## Check and add dependencies that are needed for shell programs if they are missing from the programs that are to be installed via conda.
  SHELL_PROGRAM:
    foreach
      my $shell_program ( keys %{ $parameter_href->{$installation}{shell} } )
    {
#next SHELL_PROGRAM if (undef %{ $parameter_href->{$installation}{shell}{$shell_program}{conda_dependency} });
        my $dependency_href = Dive( $parameter_href->{$installation},
            q{shell}, $shell_program, q{conda_dependency} );
        next SHELL_PROGRAM if not defined $dependency_href;
        check_and_add_dependencies(
            {
                conda_program_href => $parameter_href->{$installation}{conda},
                dependency_href    => $dependency_href,
                log                => $log,
                shell_program      => $shell_program,
            }
        );
    }

    ## Exit if a python 2 env has ben specified for a python 3 program
    check_python_compability(
        {
            installation_set_href => $parameter_href->{$installation},
            log                   => $log,
            python3_programs_ref  => $parameter_href->{python3_programs},
            python_version => $parameter_href->{$installation}{conda}{python},
            select_programs_ref => $parameter_href->{select_programs},
        }
    );

    return;
}
1;
