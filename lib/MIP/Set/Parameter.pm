package MIP::Set::Parameter;

use 5.026;
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
use List::Util qw{ any };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COLON $COMMA $CLOSE_BRACE $CLOSE_BRACKET $GENOME_VERSION $LOG_NAME
  $NEWLINE $OPEN_BRACE $OPEN_BRACKET $SPACE $TAB };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.28;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      set_conda_path
      set_default_to_active_parameter
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

    use MIP::Check::Unix qw{ is_binary_in_path };
    use MIP::Get::Parameter qw{ get_conda_path };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

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

    ## Set path to conda env
    my $environment_name = $active_parameter_href->{environment_name};
    $active_parameter_href->{conda_prefix_path} =
      catdir( $active_parameter_href->{conda_path}, q{envs}, $environment_name );

    return;
}

sub set_default_to_active_parameter {

## Function : Checks and sets user input or default values to active_parameters.
## Returns  :
## Arguments: $active_parameter_href  => Holds all set parameter for analysis
##          : $associated_recipes_ref => The parameters recipe {REF}
##          : $parameter_href         => Holds all parameters
##          : $parameter_name         => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $associated_recipes_ref;
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
        parameter_name => { defined => 1, required => 1, store => \$parameter_name, },
        case_id        => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

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
    use MIP::Constants qw{ set_genome_build_constants };

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

            ## Only set global constant once
            last GENOME_PREFIX if ($GENOME_VERSION);

            set_genome_build_constants(
                {
                    genome_version => $genome_version,
                    genome_source  => $genome_prefix,
                }
            );
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
    use Data::Diver qw{ Dive };
    use MIP::Get::Parameter qw{ get_programs_for_shell_installation };
    use MIP::Check::Installation qw{ check_and_add_dependencies };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Check that the options supplied are compatible with each other
    if (    ( scalar @{ $active_parameter_href->{skip_programs} } > 0 )
        and ( scalar @{ $active_parameter_href->{select_programs} } > 0 ) )
    {
        $log->fatal(
q{"--skip_programs" and "--select_programs" are mutually exclusive command line options}
        );
        exit 1;
    }

    ## Get programs that are to be installed via shell
    my @shell_programs_to_install = get_programs_for_shell_installation(
        {
            conda_programs_href        => $active_parameter_href->{conda},
            log                        => $log,
            prefer_shell               => $active_parameter_href->{prefer_shell},
            shell_install_programs_ref => $active_parameter_href->{shell_install},
            shell_programs_href        => $active_parameter_href->{shell},
        }
    );

    ## Remove the conda packages that has been selected to be installed via SHELL
    delete @{ $active_parameter_href->{conda} }{@shell_programs_to_install};

    ## Delete shell programs that are to be installed via conda instead of shell
    my @shell_programs_to_delete = keys %{ $active_parameter_href->{shell} };
    @shell_programs_to_delete =
      array_minus( @shell_programs_to_delete, @shell_programs_to_install );
    delete @{ $active_parameter_href->{shell} }{@shell_programs_to_delete};

    ## Solve the installation when the skip_program or select_program parameter has been used
  INSTALL_MODE:
    foreach my $install_mode (qw{ conda pip shell singularity }) {

        ## Remove programs that are to be skipped
        delete @{ $active_parameter_href->{$install_mode} }
          { @{ $active_parameter_href->{skip_programs} } };

        ## Remove all non-selected programs
        if ( scalar @{ $active_parameter_href->{select_programs} } > 0 ) {
            my @non_selects = keys %{ $active_parameter_href->{$install_mode} };
            @non_selects =
              array_minus( @non_selects, @{ $active_parameter_href->{select_programs} } );
            delete @{ $active_parameter_href->{$install_mode} }{@non_selects};
        }
    }

    ## Check and add dependencies that are needed for shell programs if they are missing from the programs that are to be installed via conda.
  SHELL_PROGRAM:
    foreach my $shell_program ( keys %{ $active_parameter_href->{shell} } ) {
        my $dependency_href =
          Dive( $active_parameter_href->{shell}, $shell_program, q{conda_dependency} );

        next SHELL_PROGRAM if ( not defined $dependency_href );
        check_and_add_dependencies(
            {
                conda_program_href => $active_parameter_href->{conda},
                dependency_href    => $dependency_href,
                log                => $log,
                shell_program      => $shell_program,
            }
        );
    }
    return;
}

sub set_recipe_mode {

## Function : Set recipe mode
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $mode                  => Mode to set
##          : $recipes_ref           => Recipes to set {REF}

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

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

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
