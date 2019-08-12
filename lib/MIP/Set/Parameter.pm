package MIP::Set::Parameter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
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
use List::Util qw{ any };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $COLON $COMMA $CLOSE_BRACE $NEWLINE $OPEN_BRACE $SPACE $TAB $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.14;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      set_conda_path
      set_conda_env_names_and_paths
      set_config_to_active_parameters
      set_custom_default_to_active_parameter
      set_default_config_dynamic_parameters
      set_default_to_active_parameter
      set_cache
      set_human_genome_reference_features
      set_nist_file_name_path
      set_no_dry_run_parameters
      set_parameter_reference_dir_path
      set_parameter_to_broadcast
      set_programs_for_installation
      set_recipe_mode
      set_recipe_resource
    };
}

## Constants
Readonly my $TWO         => 2;
Readonly my $ONE_HUNDRED => 100;

sub set_conda_path {

## Function : Set path to conda
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $log                   => Log object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Unix qw{ is_binary_in_path };
    use MIP::Get::Parameter qw{ get_conda_path };

    ## Check if conda is in path
    is_binary_in_path(
        {
            binary => q{conda},
            log    => $log,
        }
    );

    ## Get path to conda
    my $conda_path = get_conda_path( {} );

    ## Set path to conda
    $active_parameter_href->{conda_path} = $conda_path;

    return;
}

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
        parameter_name => { defined => 1, required => 1, store => \$parameter_name, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set default value only to active_parameter
    my %set_to_active_parameter = (
        analysis_type                  => \&_set_analysis_type,
        bwa_build_reference            => \&_set_human_genome,
        fusion_filter_reference_genome => \&_set_human_genome,
        gatk_path                      => \&_set_dynamic_path,
        infile_dirs                    => \&_set_infile_dirs,
        picardtools_path               => \&_set_dynamic_path,
        reference_dir                  => \&_set_reference_dir,
        rtg_vcfeval_reference_genome   => \&_set_human_genome,
        salmon_quant_reference_genome  => \&_set_human_genome,
        select_programs                => \&_set_select_programs,
        shell_install                  => \&_set_shell_install,
        skip_programs                  => \&_set_skip_programs,
        star_aln_reference_genome      => \&_set_human_genome,
        snpeff_path                    => \&_set_dynamic_path,
        temp_directory                 => \&_set_temp_directory,
        vep_directory_path             => \&_set_dynamic_path,
    );

    ## Set default value to parameter and/or active parameter
    my %set_to_parameter = (
        exome_target_bed => \&_set_capture_kit,
        sample_info_file => \&_set_sample_info_file,
    );

    if ( exists $set_to_active_parameter{$parameter_name} ) {

        $set_to_active_parameter{$parameter_name}->(
            {
                active_parameter_href => $active_parameter_href,
                parameter_name        => $parameter_name,
            }
        );
    }

    if ( exists $set_to_parameter{$parameter_name} ) {

        $set_to_parameter{$parameter_name}->(
            {
                active_parameter_href => $active_parameter_href,
                log                   => $log,
                parameter_href        => $parameter_href,
                parameter_name        => $parameter_name,
            }
        );
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
## Arguments: $active_parameter_href  => Holds all set parameter for analysis
##          : $associated_recipes_ref => The parameters recipe {REF}
##          : $log                    => Log object
##          : $parameter_href         => Holds all parameters
##          : $parameter_name         => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $associated_recipes_ref;
    my $log;
    my $parameter_href;
    my $parameter_name;

    ## Default(s)
    my $case_id;

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
        associated_recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$associated_recipes_ref,
            strict_type => 1,
        },
        log => {
            required => 1,
            defined  => 1,
            store    => \$log
        },
        parameter_name => { defined => 1, required => 1, store => \$parameter_name, },
        case_id        => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

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

            next ASSOCIATED_RECIPE
              if ( not $active_parameter_href->{$associated_recipe} );

            ## Mandatory parameter not supplied
            $log->fatal( q{Supply '-}
                  . $parameter_name
                  . q{' if you want to run }
                  . $associated_recipe );
            exit 1;
        }
    }
    return;
}

sub set_cache {

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

                carp q{Unexpected trailing garbage at end of aggregate_element '}
                  . $aggregate_element
                  . q{':}, $NEWLINE . $TAB . $unexpected_data . $NEWLINE;
            }

            if ( defined $parameter_href->{$parameter_name}{$second_key}
                && $parameter_href->{$parameter_name}{$second_key} eq $string_to_match )
            {

                push @{ $parameter_href->{cache}{$string_to_match} }, $parameter_name;
            }
        }
    }
    return;
}

sub set_human_genome_reference_features {

## Function : Detect version and source of the human_genome_reference: Source (hg19 or grch) as well as compression status.
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
        grch => qr/grch(\d+[.]\d+ | \d+)/xsm,
        hg   => qr/hg(\d+)/xsm,
    );

  GENOME_PREFIX:
    foreach my $genome_prefix ( keys %genome_source ) {

        ## Capture version
        my ($genome_version) =
          $human_genome_reference =~ m/ $genome_source{$genome_prefix}_homo_sapiens /xms;

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
              . q{Please supply the reference on this format: [sourceversion]_[species] e.g. 'grch37_homo_sapiens' or 'hg19_homo_sapiens'}
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

sub set_nist_file_name_path {

## Function : Set nist file name path by adding reference directory
## Returns  : 1
## Arguments: $active_parameter_href => Holds all set parameter for analysis
##          : $nist_parameters_ref   => Nist parameters to check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $nist_parameters_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        nist_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$nist_parameters_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Unpack
    my $reference_dir = $active_parameter_href->{reference_dir};

  NIST_PARAMETER:
    foreach my $nist_parameter ( @{$nist_parameters_ref} ) {

        # Alias
        my $nist_href = \%{ $active_parameter_href->{$nist_parameter} };

      NIST_VERSION:
        foreach my $nist_version ( keys %{$nist_href} ) {

          NIST_FILE:
            while ( my ( $nist_id, $file_name ) = each %{ $nist_href->{$nist_version} } )
            {

                ## Add reference directory to path
                $nist_href->{$nist_version}{$nist_id} =
                  catfile( $reference_dir, $file_name );
            }
        }
    }
    return 1;
}

sub set_no_dry_run_parameters {

## Function : Set parameters for true run i.e. not a dry run
## Returns  :
## Arguments: $analysis_date    => Analysis date
##          : $is_dry_run_all   => Dry run boolean
##          : $mip_version      => MIP version
##          : $sample_info_href => Info on samples and case hash {REF}

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
        my $path = catfile( $reference_dir, $active_parameter_href->{$parameter_name} );
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
            my $element_separator = $parameter_href->{$parameter_name}{element_separator};

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

                ## Hash of hash
                if ( ref $value eq q{HASH} ) {

                    $info .= $OPEN_BRACE . $key . q{ => };

                  HASH:
                    foreach my $sec_key ( keys %{$value} ) {

                        $info .= $sec_key . q{=};
                        if ( $value->{$sec_key} ) {

                            $info .= $value->{$sec_key};
                        }
                        $info .= $COMMA;
                    }
                    $info .= $CLOSE_BRACE . $SPACE;
                }
                ## Hash of array
                elsif ( ref $value eq q{ARRAY} ) {

                    $info .= join $COMMA, map {
                        qq{$_=} . join $SPACE,
                          @{ $active_parameter_href->{$parameter_name}{$_} }
                    } ( keys %{ $active_parameter_href->{$parameter_name} } );

                    last PARAMETER_KEY;
                }
                else {

                    $info .= join $COMMA,
                      map { qq{$_=$active_parameter_href->{$parameter_name}{$_}} }
                      ( keys %{ $active_parameter_href->{$parameter_name} } );

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
## Arguments: $parameter_href => The entire parameter hash {REF}

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

    use MIP::Parse::Parameter qw{ parse_conda_env_name };

    ## A default name on which to build environment names if non specified
    my $base_name = $active_parameter_href->{environment_base_name};

    ## Get date and reformat to six digits
    my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime time;
    my $date = sprintf '%02d%02d%02d', $year % $ONE_HUNDRED, $mon + 1, $mday;

    ## Set up conda environment names and prefix paths for environments
  ENVIRONMENT:
    foreach my $environment ( @{ $active_parameter_href->{installations} } ) {

        ## Construct conda environment name
        my $environment_name = parse_conda_env_name(
            {
                base_name      => $base_name,
                date           => $date,
                environment    => $environment,
                parameter_href => $active_parameter_href,
            }
        );

        ## Set names and paths
        $active_parameter_href->{environment_name}{$environment} = $environment_name;
        $active_parameter_href->{$environment}{conda_prefix_path} =
          catdir( $active_parameter_href->{conda_path}, q{envs}, $environment_name );

    }
    return;
}

sub set_programs_for_installation {

## Function : Process the lists of programs that has been selected for or omitted from installation
##          : and update the environment packages
## Returns  :
## Arguments: $installation          => Environment to be installed
##          : $log                   => Log
##          : $active_parameter_href => The entire active parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $installation;
    my $log;
    my $active_parameter_href;

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
    use Data::Diver qw{ Dive };
    use MIP::Get::Parameter qw{ get_programs_for_shell_installation };
    use MIP::Check::Installation
      qw{ check_and_add_dependencies check_python_compability };

    ## Check that the options supplied are compatible with each other
    if (    ( scalar @{ $active_parameter_href->{skip_programs} } > 0 )
        and ( scalar @{ $active_parameter_href->{select_programs} } > 0 ) )
    {
        $log->fatal(
q{"--skip_programs" and "--select_programs" are mutually exclusive command line options}
        );
        exit 1;
    }

    ## Check that only one environment has been specified for installation if the option select_program is used
    if (    ( scalar @{ $active_parameter_href->{select_programs} } > 0 )
        and ( scalar @{ $active_parameter_href->{installations} } > 1 ) )
    {
        $log->fatal(
q{Please select a single installation environment when using the option --select_programs}
        );
        exit 1;
    }

    ## Get programs that are to be installed via shell
    my @shell_programs_to_install = get_programs_for_shell_installation(
        {
            conda_programs_href        => $active_parameter_href->{$installation}{conda},
            log                        => $log,
            prefer_shell               => $active_parameter_href->{prefer_shell},
            shell_install_programs_ref => $active_parameter_href->{shell_install},
            shell_programs_href        => $active_parameter_href->{$installation}{shell},
        }
    );

    ## Remove the conda packages that has been selected to be installed via SHELL
    delete @{ $active_parameter_href->{$installation}{conda} }
      {@shell_programs_to_install};

    ## Special case for snpsift since it is installed together with SnpEff
    ## if shell installation of SnpEff has been requested.
    if ( any { $_ eq q{snpeff} } @shell_programs_to_install ) {
        delete $active_parameter_href->{$installation}{conda}{snpsift};
    }
    ## Store variable outside of shell hash and use Data::Diver module to avoid autovivification of variable
    $active_parameter_href->{$installation}{snpeff_genome_versions} = Dive(
        $active_parameter_href->{$installation},
        qw{ shell snpeff snpeff_genome_versions }
    );

    ## Delete shell programs that are to be installed via conda instead of shell
    my @shell_programs_to_delete =
      keys %{ $active_parameter_href->{$installation}{shell} };
    @shell_programs_to_delete =
      array_minus( @shell_programs_to_delete, @shell_programs_to_install );
    delete @{ $active_parameter_href->{$installation}{shell} }{@shell_programs_to_delete};

    ## Solve the installation when the skip_program or select_program parameter has been used
  INSTALL_MODE:
    foreach my $install_mode (qw{ conda pip shell }) {

        ## Remove programs that are to be skipped
        delete @{ $active_parameter_href->{$installation}{$install_mode} }
          { @{ $active_parameter_href->{skip_programs} } };

        ## Remove all non-selected programs
        if ( scalar @{ $active_parameter_href->{select_programs} } > 0 ) {
            my @non_selects =
              keys %{ $active_parameter_href->{$installation}{$install_mode} };
            @non_selects =
              array_minus( @non_selects, @{ $active_parameter_href->{select_programs} } );
            delete @{ $active_parameter_href->{$installation}{$install_mode} }
              {@non_selects};
        }
    }

    ## Check and add dependencies that are needed for shell programs if they are missing from the programs that are to be installed via conda.
  SHELL_PROGRAM:
    foreach my $shell_program ( keys %{ $active_parameter_href->{$installation}{shell} } )
    {
        my $dependency_href = Dive( $active_parameter_href->{$installation},
            q{shell}, $shell_program, q{conda_dependency} );
        next SHELL_PROGRAM if not defined $dependency_href;
        check_and_add_dependencies(
            {
                conda_program_href => $active_parameter_href->{$installation}{conda},
                dependency_href    => $dependency_href,
                log                => $log,
                shell_program      => $shell_program,
            }
        );
    }

    ## Exit if a python 2 env has ben specified for a python 3 program
    check_python_compability(
        {
            installation_set_href => $active_parameter_href->{$installation},
            log                   => $log,
            python3_programs_ref  => $active_parameter_href->{python3_programs},
            python_version      => $active_parameter_href->{$installation}{conda}{python},
            select_programs_ref => $active_parameter_href->{select_programs},
        }
    );

    return;
}

sub set_recipe_mode {

## Function : Set recipe mode
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $log                   => Log
##          : $mode                  => Mode to set
##          : $recipes_ref           => Recipes to set {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
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
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        mode => {
            allow       => [ 0, 1, $TWO ],
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

    ## Set recipe mode
  RECIPE:
    foreach my $recipe ( @{$recipes_ref} ) {

        $active_parameter_href->{$recipe} = $mode;

        ## Broadcast
        $log->info(
            q{Set} . $SPACE . $recipe . $SPACE . q{to} . $COLON . $SPACE . $mode );
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
        while ( my ( $recipe, $core_number ) =
            each %{ $active_parameter_href->{$set_hash_key} } )
        {

            $active_parameter_href->{$target_hash_key}{$recipe} = $core_number;
        }
    }
    return;
}

sub _set_analysis_type {

## Function : Set default analysis type to active parameters
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

    map { $active_parameter_href->{$parameter_name}{$_} = q{wgs} }
      @{ $active_parameter_href->{sample_ids} };
    return;
}

sub _set_capture_kit {

## Function : Set default capture kit to active parameters
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $log                   => Log object
##          : $parameter_href        => Holds all parameters {REF}
##          : $parameter_name        => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
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
        log => {
            required => 1,
            defined  => 1,
            store    => \$log
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

    use MIP::Get::Parameter qw{ get_capture_kit };

    ### If capture kit is not set after cmd, config and reading pedigree
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
      {$capture_kit} = join $COMMA,
      @{ $active_parameter_href->{sample_ids} };

    $log->warn(
        q{Could not detect a supplied capture kit. Will Try to use 'latest' capture kit: }
          . $capture_kit );
    return;
}

sub _set_dynamic_path {

## Function : Set default dynamic paths to active parameters
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

    use MIP::Get::Parameter qw{ get_dynamic_conda_path };

    ## Already has a path set
    return if ( $active_parameter_href->{$parameter_name} );

    ## Set default dynamic path if needed
    my %dynamic_path = (
        gatk_path => {
            bin_file        => q{gatk3},
            environment_key => q{gatk},
        },
        picardtools_path => {
            bin_file        => q{picard.jar},
            environment_key => q{picard},
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

    ## No defined bin_file or environment key
    return if ( not exists $dynamic_path{$parameter_name} );

    $active_parameter_href->{$parameter_name} = get_dynamic_conda_path(
        {
            active_parameter_href => $active_parameter_href,
            bin_file              => $dynamic_path{$parameter_name}{bin_file},
            environment_key       => $dynamic_path{$parameter_name}{environment_key},
        }
    );
    return;
}

sub _set_human_genome {

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

sub _set_infile_dirs {

## Function : Set default infile dirs to active parameters
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

    ## Build default for infile_dirs
  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        if ( not exists $active_parameter_href->{analysis_type}{$sample_id} ) {

            _set_analysis_type(
                {
                    active_parameter_href => $active_parameter_href,
                    parameter_name        => q{analysis_type},
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

        $active_parameter_href->{$parameter_name}{$path} = $sample_id;
    }
    return;
}

sub _set_reference_dir {

## Function : Set default reference dir to active parameters
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

    $active_parameter_href->{$parameter_name} = cwd();
    return;
}

sub _set_sample_info_file {

## Function : Set default sample_info_file and qccollect_sampleinfo_file to parameters
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $log                   => Log object
##          : $parameter_href        => Holds all parameters {REF}
##          : $parameter_name        => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
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
        log => {
            store => \$log
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

sub _set_temp_directory {

## Function : Set default temp directory to active parameters
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

sub _set_select_programs {

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

    if ( exists $active_parameter_href->{select_programs} ) {

        return;
    }

    $active_parameter_href->{select_programs} = [];

    return;
}

sub _set_shell_install {

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

    if ( exists $active_parameter_href->{shell_install} ) {

        return;
    }

    $active_parameter_href->{shell_install} = [];

    return;
}

sub _set_skip_programs {

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

    if ( exists $active_parameter_href->{skip_programs} ) {

        return;
    }

    $active_parameter_href->{skip_programs} = [];

    return;
}
1;
