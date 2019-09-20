package MIP::Update::Parameters;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
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
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ update_dynamic_config_parameters
      update_exome_target_bed
      update_reference_parameters
      update_vcfparser_outfile_counter };
}

sub update_dynamic_config_parameters {

## Function : Updates the config file to particular user/cluster for dynamic config parameters following specifications. Leaves other entries untouched.
## Returns  :
## Arguments: $active_parameter_href  => Active parameters for this analysis hash {REF}
##          : $dynamic_parameter_href => Map of dynamic parameters
##          : $parameter_name         => MIP Parameter to update

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $dynamic_parameter_href;
    my $parameter_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        dynamic_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$dynamic_parameter_href,
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

    return if ( not defined $active_parameter_href->{$parameter_name} );

    if ( ref $active_parameter_href->{$parameter_name} eq q{HASH} ) {

      KEY:
        foreach my $key ( keys %{ $active_parameter_href->{$parameter_name} } ) {

          DYNAMIC_PARAMETER:
            while ( my ( $dynamic_parameter_name, $dynamic_parameter_value ) =
                each %{$dynamic_parameter_href} )
            {

                ## Replace dynamic config parameters with actual value that is now set from cmd or config
	      next KEY if(not $active_parameter_href->{$parameter_name}{$key});

                $active_parameter_href->{$parameter_name}{$key} =~
                  s/$dynamic_parameter_value!/$dynamic_parameter_value/smgi;
            }

            update_dynamic_config_parameters(
                {
                    active_parameter_href  => $active_parameter_href->{$parameter_name},
                    dynamic_parameter_href => $dynamic_parameter_href,
                    parameter_name         => $key,
                }
            );
        }
    }

  DYNAMIC_PARAMETER:
    while ( my ( $dynamic_parameter_name, $dynamic_parameter_value ) =
        each %{$dynamic_parameter_href} )
    {

        ## Replace dynamic config parameters with actual value that is now set from cmd or config
        $active_parameter_href->{$parameter_name} =~
          s/$dynamic_parameter_name!/$dynamic_parameter_value/smgi;
    }
    return;
}

sub update_exome_target_bed {

## Function : Update exome_target_bed files with human genome reference source and version
## Returns  :
## Arguments: $exome_target_bed_file_href     => Exome target bed
##          : $human_genome_reference_source  => Human genome reference source
##          : $human_genome_reference_version => Human genome reference version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $exome_target_bed_file_href;
    my $human_genome_reference_source;
    my $human_genome_reference_version;

    my $tmpl = {
        exome_target_bed_file_href =>
          { required => 1, store => \$exome_target_bed_file_href, },
        human_genome_reference_source => {
            defined     => 1,
            required    => 1,
            store       => \$human_genome_reference_source,
            strict_type => 1,
        },
        human_genome_reference_version => {
            defined     => 1,
            required    => 1,
            store       => \$human_genome_reference_version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  EXOME_FILE:
    foreach my $exome_target_bed_file ( keys %{$exome_target_bed_file_href} ) {

        my $original_file_name = $exome_target_bed_file;

        ## Replace with actual version
        if ( $exome_target_bed_file =~
            s/genome_reference_source/$human_genome_reference_source/xsm
            && $exome_target_bed_file =~ s/_version/$human_genome_reference_version/xsm )
        {

            ## The delete operator returns the value being deleted i.e. updating hash key while preserving original info
            $exome_target_bed_file_href->{$exome_target_bed_file} =
              delete $exome_target_bed_file_href->{$original_file_name};
        }
    }
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

    use MIP::Set::Parameter qw{ set_parameter_reference_dir_path };

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

sub update_vcfparser_outfile_counter {

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
        sv_vcfparser => { sv_vcfparser_select_file => q{sv_vcfparser_outfile_count} },
        vcfparser_ar => { vcfparser_select_file    => q{vcfparser_outfile_count} },
    );

## Determine if to expect select outfile for vcfparser and sv_vcfparser
  RECIPE:
    foreach my $recipe ( keys %vcfparser_select_file ) {

        next RECIPE if ( not $active_parameter_href->{$recipe} );

      FILES:
        while ( my ( $parameter_name, $parameter_name_counter ) =
            each %{ $vcfparser_select_file{$recipe} } )
        {

            $active_parameter_href->{$parameter_name_counter} =
              _set_vcfparser_file_counter(
                {
                    parameter_name => $active_parameter_href->{$parameter_name},
                }
              );
        }
    }
    return;
}

sub _set_vcfparser_file_counter {

## Function : Return the expected number of outputfile(s) after vcfparser
## Returns  : 1 | 2
## Arguments: $parameter_name => Vcfparser select file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_name;

    my $tmpl =
      { parameter_name => { required => 1, store => \$parameter_name, strict_type => 1, },
      };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## To track if vcfparser was used with a vcfparser_select_file (=2) or not (=1)
    # No select file was given
    return 1 if ( not defined $parameter_name );

    ## Select file was given
    return 2;
}

1;
