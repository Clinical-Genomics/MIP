package MIP::Set::Pedigree;

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
use List::MoreUtils qw { any };
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ set_active_parameter_pedigree_keys set_pedigree_capture_kit_info set_pedigree_case_info set_pedigree_phenotype_info set_pedigree_sample_info set_pedigree_sex_info };
}

## Constants
Readonly my $SPACE => q{ };

sub set_active_parameter_pedigree_keys {

## Function : Get the pedigree case keys and values
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $pedigree_href           => YAML pedigree info hash {REF}
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $user_supply_switch_href => The user supplied info switch {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $pedigree_href;
    my $sample_info_href;
    my $user_supply_switch_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        user_supply_switch_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$user_supply_switch_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @pedigree_keys = qw{ analysis_type expected_coverage time_point };
    my @sample_ids    = @{ $pedigree_href->{samples} };

  SAMPLE_HREF:
    foreach my $pedigree_sample_href ( @{ $pedigree_href->{samples} } ) {

        ## Alias
        my $sample_id = $pedigree_sample_href->{sample_id};

      PEDIGREE_KEY:
        foreach my $pedigree_key (@pedigree_keys) {

            ## Alias
            my $pedigree_value =
              $sample_info_href->{sample}{$sample_id}{$pedigree_key};

            ## Key was not set in pedigree (this run or previous run)
            ## Hence not stored in sample_info
            next PEDIGREE_KEY if ( not defined $pedigree_value );

            ## User supplied info on cmd or config
            next PEDIGREE_KEY if ( $user_supply_switch_href->{$pedigree_key} );

            ## Add value for sample_id using pedigree info
            $active_parameter_href->{$pedigree_key}{$sample_id} = $pedigree_value;
        }
    }
    return;
}

sub set_pedigree_capture_kit_info {

## Function : Add and set capture kit for each individual
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $pedigree_href           => YAML pedigree info hash {REF}
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $user_supply_switch_href => The user supplied info switch {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;
    my $pedigree_href;
    my $sample_info_href;
    my $user_supply_switch_href;

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
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        user_supply_switch_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$user_supply_switch_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::Parameter qw{ get_capture_kit };

    my %exom_target_bed_file_tracker;

  SAMPLE_HREF:
    foreach my $pedigree_sample_href ( @{ $pedigree_href->{samples} } ) {

        ## Alias
        my $sample_id = $pedigree_sample_href->{sample_id};
        my $capture_kit =
          $sample_info_href->{sample}{$sample_id}{capture_kit};

        ## No recorded capture kit from pedigree or previous run
        next SAMPLE_HREF if ( not $capture_kit );

        ## Return a capture kit depending on user info
        my $exome_target_bed_file = get_capture_kit(
            {
                capture_kit => $capture_kit,
                supported_capture_kit_href =>
                  $parameter_href->{supported_capture_kit}{default},
                user_supplied_parameter_switch =>
                  $user_supply_switch_href->{exome_target_bed},
            }
        );

        ## No capture kit returned
        next SAMPLE_HREF if ( not $exome_target_bed_file );

        push @{ $exom_target_bed_file_tracker{$exome_target_bed_file} }, $sample_id;
    }

  BED_FILE:
    foreach my $exome_target_bed_file ( keys %exom_target_bed_file_tracker ) {

        ## We have read capture kits from pedigree and
        ## need to transfer to active_parameters
        $active_parameter_href->{exome_target_bed}{$exome_target_bed_file}
          = join q{,},
          @{ $exom_target_bed_file_tracker{$exome_target_bed_file} };
    }
    return;
}

sub set_pedigree_case_info {

## Function : Get the pedigree case keys and values
## Returns  :
## Arguments: $pedigree_href    => YAML pedigree info hash {REF}
##          : $sample_info_href => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $pedigree_href;
    my $sample_info_href;

    my $tmpl = {
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
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

    ### Add key and values for case level info
  KEY:
    foreach my $key ( keys %{$pedigree_href} ) {

        ## Do not add sample level info
        next KEY if ( $key eq q{samples} );

        $sample_info_href->{$key} = $pedigree_href->{$key};
    }
    return;
}

sub set_pedigree_phenotype_info {

    ## Function : Store phenotype and plink phenotype in dynamic parameters. Reformats pedigree phenotype to plink format before adding
## Returns  :
## Arguments: $parameter_href => Parameter hash {REF}
##          : $pedigree_href  => YAML pedigree info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $pedigree_href;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %plink_phenotype = (
        unaffected => 1,
        affected   => 2,
        other      => q{other},
        unknown    => 0,
    );
  SAMPLE_HREF:
    foreach my $pedigree_sample_href ( @{ $pedigree_href->{samples} } ) {

        next SAMPLE_HREF if ( not exists $pedigree_sample_href->{phenotype} );

        # Alias
        my $sample_id = $pedigree_sample_href->{sample_id};
        my $phenotype = $pedigree_sample_href->{phenotype};

        push
          @{ $parameter_href->{cache}{ $pedigree_sample_href->{phenotype} } },
          $sample_id;

        if ( exists $plink_phenotype{$phenotype} ) {

            $parameter_href->{cache}{$sample_id}{plink_phenotype} =
              $plink_phenotype{$phenotype};
        }
    }
    return;
}

sub set_pedigree_sample_info {

## Function : Get the pedigree sample keys and values
## Returns  :
## Arguments: $active_parameter_href     => Active parameters for this analysis hash {REF}
##          : $pedigree_href             => YAML pedigree info hash {REF}
##          : $sample_info_href          => Info on samples and case hash {REF}
##          : $user_supply_switch_href   => User supplied info switch {REF}
##          : $user_input_sample_ids_ref => User supplied sample_ids via cmd or config

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $pedigree_href;
    my $sample_info_href;
    my $user_supply_switch_href;
    my $user_input_sample_ids_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        user_supply_switch_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$user_supply_switch_href,
            strict_type => 1,
        },
        user_input_sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$user_input_sample_ids_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Store pedigree sample_ids
    my @pedigree_sample_ids;

    ### Add values sample level info
  SAMPLE_HREF:
    foreach my $pedigree_sample_href ( @{ $pedigree_href->{samples} } ) {

        # Alias
        my $sample_id = $pedigree_sample_href->{sample_id};

        ## Save pedigree sample_id info for later check
        push @pedigree_sample_ids, $sample_id;

        ## No sample_ids supplied by user
        if ( not $user_supply_switch_href->{sample_ids} ) {

            ## Save sample_id info for analysis
            push @{ $active_parameter_href->{sample_ids} }, $sample_id;

            ## Add input to sample_info hash for at sample level
            _add_pedigree_sample_info(
                {
                    pedigree_sample_href => $pedigree_sample_href,
                    sample_id            => $sample_id,
                    sample_info_href     => $sample_info_href,
                }
            );
            next SAMPLE_HREF;
        }

        ### Sample_ids supplied by user
        ## Update sample_id with pedigree info for user supplied sample_ids
        if ( any { $_ eq $sample_id } @{$user_input_sample_ids_ref} ) {

            ## Set pedigree sample info
            _add_pedigree_sample_info(
                {
                    pedigree_sample_href => $pedigree_sample_href,
                    sample_id            => $sample_id,
                    sample_info_href     => $sample_info_href,
                }
            );
        }
    }
    return @pedigree_sample_ids;
}

sub set_pedigree_sex_info {

    ## Function : Store sex and plink sex in dynamic parameters. Reformats pedigree sex to plink format before adding
## Returns  :
## Arguments: $parameter_href => Parameter hash {REF}
##          : $pedigree_href  => YAML pedigree info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $pedigree_href;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %plink_sex = (
        male    => 1,
        female  => 2,
        other   => q{other},
        unknown => q{other},
    );
  SAMPLE_HREF:
    foreach my $pedigree_sample_href ( @{ $pedigree_href->{samples} } ) {

        next SAMPLE_HREF if ( not exists $pedigree_sample_href->{sex} );

        # Alias
        my $sample_id = $pedigree_sample_href->{sample_id};
        my $sex       = $pedigree_sample_href->{sex};

        push @{ $parameter_href->{cache}{ $pedigree_sample_href->{sex} } }, $sample_id;

        if ( exists $plink_sex{$sex} ) {

            $parameter_href->{cache}{$sample_id}{plink_sex} =
              $plink_sex{$sex};
        }
    }
    return;
}

sub _add_pedigree_sample_info {

## Function : Add pedigree sample level keys and values to sample info
## Returns  :
## Arguments: $pedigree_sample_href => YAML sample info hash {REF}
##          : $sample_id            => Sample ID
##          : $sample_info_href     => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $pedigree_sample_href;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        pedigree_sample_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_sample_href,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
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

    ## Add input to sample_info hash for at sample level
  PEDIGREE_SAMPLE_KEY:
    foreach my $key ( keys %{$pedigree_sample_href} ) {

        $sample_info_href->{sample}{$sample_id}{$key} =
          $pedigree_sample_href->{$key};
    }
    return;
}

1;
