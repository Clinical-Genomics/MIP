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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ set_recipe_on_analysis_type };
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
    use MIP::Recipes::Analysis::Vep
      qw{ analysis_vep_sv_wes analysis_vep_sv_wgs };

    my %analysis_type_recipe = (
        wes => {
            gatk_variantrecalibration =>
              \&analysis_gatk_variantrecalibration_wes,
            sv_vcfparser              => \&analysis_vcfparser_sv_wes,
            sv_varianteffectpredictor => \&analysis_vep_sv_wes,
        },
        wgs => {
            gatk_variantrecalibration =>
              \&analysis_gatk_variantrecalibration_wgs,
            sv_vcfparser              => \&analysis_vcfparser_sv_wgs,
            sv_varianteffectpredictor => \&analysis_vep_sv_wgs,
        },
    );

  PROGRAM:
    while ( my ( $program_name, $recipe_cref ) =
        each %{ $analysis_type_recipe{$consensus_analysis_type} } )
    {

        next PROGRAM if ( not exists $analysis_recipe_href->{$program_name} );

        if ( exists $analysis_type_recipe{$consensus_analysis_type} ) {
            $analysis_recipe_href->{$program_name} = $recipe_cref;
        }
        else {

            ## Use wgs as default
            my $recipe_wgs_cref = $analysis_type_recipe{wgs}{$program_name};

            # Set recipe
            $analysis_recipe_href->{$program_name} = $recipe_wgs_cref;
        }
    }
    return;
}

1;
