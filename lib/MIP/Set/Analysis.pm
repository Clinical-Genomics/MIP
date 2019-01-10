package MIP::Set::Analysis;

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

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ set_recipe_on_analysis_type set_rankvariants_ar };
}

## Constants
Readonly my $SPACE => q{ };

sub set_recipe_on_analysis_type {

## Function : Set which recipe to use dependeing on consensus analysis type
## Returns  :
## Arguments: $analysis_recipe_href    => Analysis recipe hash {REF}
##          : $consensus_analysis_type => Consensus analysis type

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_recipe_href;
    my $consensus_analysis_type;

    my $tmpl = {
        analysis_recipe_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$analysis_recipe_href,
            strict_type => 1,
        },
        consensus_analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$consensus_analysis_type,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Recipes::Analysis::Gatk_variantrecalibration
      qw{ analysis_gatk_variantrecalibration_wes analysis_gatk_variantrecalibration_wgs };
    use MIP::Recipes::Analysis::Mip_vcfparser
      qw{ analysis_vcfparser_sv_wes analysis_vcfparser_sv_wgs };
    use MIP::Recipes::Analysis::Vep qw{ analysis_vep_sv_wes analysis_vep_sv_wgs };

    my %analysis_type_recipe = (
        vrn => {
            sv_varianteffectpredictor => \&analysis_vep_sv_wgs,
            sv_vcfparser              => \&analysis_vcfparser_sv_wgs,
        },
        wes => {
            gatk_variantrecalibration => \&analysis_gatk_variantrecalibration_wes,
            sv_varianteffectpredictor => \&analysis_vep_sv_wes,
            sv_vcfparser              => \&analysis_vcfparser_sv_wes,
        },
        wgs => {
            gatk_variantrecalibration => \&analysis_gatk_variantrecalibration_wgs,
            sv_varianteffectpredictor => \&analysis_vep_sv_wgs,
            sv_vcfparser              => \&analysis_vcfparser_sv_wgs,
        },
    );

    ## If not a defined consensus analysis type e.g. "mixed"
    if ( not exists $analysis_type_recipe{$consensus_analysis_type} ) {

        ## Use wgs as fallback
        $consensus_analysis_type = q{wgs};
    }

  ANALYSIS_RECIPE:
    foreach my $recipe_name ( keys %{$analysis_recipe_href} ) {

        next ANALYSIS_RECIPE
          if ( not exists $analysis_type_recipe{$consensus_analysis_type}{$recipe_name} );

        $analysis_recipe_href->{$recipe_name} =
          $analysis_type_recipe{$consensus_analysis_type}{$recipe_name};
    }
    return;
}

sub set_rankvariants_ar {

## Function : Update which rankvariants recipe to use
## Returns  :
## Arguments: $sample_ids_ref   => Active parameters for this analysis hash {REF}
##          : $analysis_recipe_href    => Analysis recipe hash {REF}
##          : $log                     => Log object to write to
##          : $parameter_href          => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_ids_ref;
    my $analysis_recipe_href;
    my $log;
    my $parameter_href;

    my $tmpl = {
        analysis_recipe_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$analysis_recipe_href,
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
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Recipes::Analysis::Rankvariant
      qw{ analysis_rankvariant analysis_rankvariant_unaffected analysis_rankvariant_sv analysis_rankvariant_sv_unaffected };

    if ( defined $parameter_href->{cache}{unaffected}
        && @{ $parameter_href->{cache}{unaffected} } eq @{$sample_ids_ref} )
    {

        $log->warn(
q{Only unaffected sample(s) in pedigree - skipping genmod 'models', 'score' and 'compound'}
        );

        $analysis_recipe_href->{sv_rankvariant} = \&analysis_rankvariant_sv_unaffected;
        $analysis_recipe_href->{rankvariant}    = \&analysis_rankvariant_unaffected;
    }
    else {

        $analysis_recipe_href->{sv_rankvariant} = \&analysis_rankvariant_sv;
        $analysis_recipe_href->{rankvariant}    = \&analysis_rankvariant;
    }
    return;
}

1;
