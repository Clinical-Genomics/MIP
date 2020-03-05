package MIP::Update::Recipes;

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

## MIPs lib/
use MIP::Constants qw{ $COMMA $LOG_NAME $NEWLINE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      update_prioritize_flag
      update_recipe_mode_for_analysis_type
      update_recipe_mode_for_pedigree
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

1;
