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
use MIP::Constants qw{ $COMMA $LOG_NAME $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.08;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      update_recipe_mode_for_pedigree
    };
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
