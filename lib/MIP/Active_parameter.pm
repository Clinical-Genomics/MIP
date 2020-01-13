package MIP::Active_parameter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
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
use MIP::Constants qw{ $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_user_supplied_pedigree_parameter
      set_pedigree_sample_id_parameter
      update_to_absolute_path };
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
