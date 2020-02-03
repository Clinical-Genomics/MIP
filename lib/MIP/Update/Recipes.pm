package MIP::Update::Recipes;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use List::MoreUtils qw { any };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## MIPs lib/
use MIP::Constants qw{ $COMMA $LOG_NAME $NEWLINE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      update_prioritize_flag
      update_recipe_mode_for_analysis_type
      update_recipe_mode_for_pedigree
      update_recipe_mode_with_dry_run_all
      update_recipe_mode_with_start_with
    };
}

sub update_prioritize_flag {

##Function : Update prioritize flag depending on analysis run value as some recipes are not applicable for e.g. wes
##Returns  :
##Arguments: $consensus_analysis_type => Consensus analysis_type
##         : $prioritize_key          => Prioritize key to update
##         : $recipes_ref             => Recipes to update {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $consensus_analysis_type;
    my $prioritize_key;
    my $recipes_ref;

    my $tmpl = {
        consensus_analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$consensus_analysis_type,
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

    if (   $consensus_analysis_type ne q{wgs}
        && $prioritize_key )
    {

        ## Split string into array
        my @callers = split $COMMA, $prioritize_key;

      CALLER:
        foreach my $caller ( @{$recipes_ref} ) {

            ## Remove all wgs specific callers
            @callers = grep { $caller !~ /^$_/sxm } @callers;
        }

        ## Update sv_svdb_merge_prioritize flag
        $prioritize_key = join q{,}, @callers;
    }
    return $prioritize_key;
}

sub update_recipe_mode_for_analysis_type {

##Function : Update recipe mode depending on analysis run value as some recipes are not applicable for e.g. wes
##Returns  :
##Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##         : $consensus_analysis_type => Consensus analysis_type
##         : $log                     => Log
##         : $recipes_ref             => Recipes to update {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $consensus_analysis_type;
    my $log;
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
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
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

    if ( $consensus_analysis_type ne q{wgs} ) {

        my @warning_msgs;

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

            return @warning_msgs;
        }
    }
    return;
}

sub update_recipe_mode_for_pedigree {

## Function : Update recipe mode depending on analysis run value as some recipes are not applicable for e.g. wts
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $recipes_ref           => Recipes to update {REF}
##          : $sample_info_href      => Sample info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $recipes_ref;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$recipes_ref,
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

    use MIP::Sample_info qw{ get_pedigree_sample_id_attributes };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %phenotype_count;

  SAMPLE_ID:
    foreach my $sample_id ( keys %{ $sample_info_href->{sample} } ) {

        ## Get phenotype
        my $phenotype = get_pedigree_sample_id_attributes(
            {
                attribute        => q{phenotype},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );
        $phenotype_count{$phenotype}++;
    }

    ## Return if the phenotypes are affected and unaffected
    if (    not $phenotype_count{unknown}
        and $phenotype_count{affected}
        and $phenotype_count{unaffected} )
    {

        return;
    }

    ## Otherwise turn off and warn
  RECIPE:
    foreach my $recipe ( @{$recipes_ref} ) {

        ## Don't warn when recipe is already turned off
        next RECIPE if $active_parameter_href->{$recipe} == 0;

        ## Update recipe mode
        $active_parameter_href->{$recipe} = 0;

        $log->warn(
            q{Turned off: } . $recipe . q{ as it is not compatible with this pedigree.} );
    }
    return;
}

sub update_recipe_mode_with_dry_run_all {

## Function : Update recipe mode depending on dry_run_all flag
## Returns  :
## Arguments: $active_parameter_href => The active parameters for this analysis hash {REF}
##          : $dry_run_all           => Simulation mode
##          : $recipes_ref           => Recipes in MIP

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $dry_run_all;
    my $recipes_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        dry_run_all => {
            allow       => [ undef, 0, 1, 2 ],
            default     => 0,
            store       => \$dry_run_all,
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

    ## Activated simulation mode
    my $simulation_mode = 2;

    if ($dry_run_all) {

      RECIPE:
        foreach my $recipe_name ( @{$recipes_ref} ) {

            ## If recipe is activated
            if ( $active_parameter_href->{$recipe_name} ) {

                # Change recipe mode to simulation
                $active_parameter_href->{$recipe_name} = $simulation_mode;
            }
        }
    }
    return;
}

sub update_recipe_mode_with_start_with {

## Function : Update recipe mode depending on start with flag
## Returns  :
## Arguments: $active_parameter_href   => The active parameters for this analysis hash {REF}
##          : $recipes_ref             => Recipes in MIP
##          : $start_with_recipes_ref  => Recipes to run

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $recipes_ref;
    my $start_with_recipes_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$recipes_ref,
            strict_type => 1,
        },
        start_with_recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$start_with_recipes_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Active run mode
    my $run_mode = 1;

    ## Activated simulation mode
    my $simulation_mode = 2;

  RECIPE:
    foreach my $recipe_name ( @{$recipes_ref} ) {

        next RECIPE if ( not $active_parameter_href->{$recipe_name} );

        ## If recipe is uppstream of start recipe
        if ( not any { $_ eq $recipe_name } @{$start_with_recipes_ref} ) {

            # Change recipe mode to simulation
            $active_parameter_href->{$recipe_name} = $simulation_mode;
        }
        else {
            #Recipe or downstream dependency recipe

            # Change recipe mode to active
            $active_parameter_href->{$recipe_name} = $run_mode;
        }
    }
    return;
}

1;
