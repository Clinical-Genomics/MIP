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
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COMMA $DOT $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      get_user_supplied_pedigree_parameter
      set_default_analysis_type
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
      set_parameter_reference_dir_path
      set_pedigree_sample_id_parameter
      update_reference_parameters
      update_to_absolute_path
    };
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

sub set_default_analysis_type {

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

            set_default_analysis_type(
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

sub set_default_pedigree_fam_file {

## Function : Set default pedigree_fam_file to active parameters
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

    ## Set pedigree fam file
    $active_parameter_href->{$parameter_name} = catfile(
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

    return if ( $active_parameter_href->{$parameter_name} );

    $active_parameter_href->{$parameter_name} =
      catfile( $Bin, qw{templates program_test_cmds.yaml } );

    return;
}

sub set_default_reference_dir {

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

    ## Set reference dir to current working dir
    $active_parameter_href->{$parameter_name} = cwd();
    return;
}

sub set_default_reference_info_file {

## Function : Set default reference_info_file
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

    ## Set reference info file
    $active_parameter_href->{reference_info_file} =
      catfile( $active_parameter_href->{outdata_dir}, q{reference_info.yaml} );
    return;
}

sub set_default_store_file {

## Function : Set default store_file to active parameters
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
        parameter_name =>
          { defined => 1, required => 1, store => \$parameter_name, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set store file
    $active_parameter_href->{$parameter_name} =
      catfile( $active_parameter_href->{outdata_dir}, q{store_info.yaml} );
    return;
}

sub set_default_temp_directory {

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

    ## Build default for vcfparser select file
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

sub update_reference_parameters {

## Function : Update reference parameters with mip_reference directory path
## Returns  :
## Arguments: $active_parameter_href   => Holds all set parameter for analysis
##          : $associated_recipes_ref  => The parameters recipe(s) {REF}
##          : $parameter_name          => Parameter name

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

1;
