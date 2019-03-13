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

## MIPs lib/
use MIP::Constants qw{ $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ set_recipe_on_analysis_type set_recipe_bwa_mem set_recipe_gatk_variantrecalibration set_rankvariants_ar };
}

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

sub set_recipe_bwa_mem {

## Function : Set correct bwa_mem recipe depending on version and source of the human_genome_reference: Source (hg19 or GRCh)
## Returns  :
## Arguments: $analysis_recipe_href           => Analysis recipe hash {REF}
##          : $human_genome_reference_source  => Human genome reference source
##          : $human_genome_reference_version => Human genome reference version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_recipe_href;
    my $human_genome_reference_source;
    my $human_genome_reference_version;

    my $tmpl = {
        analysis_recipe_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$analysis_recipe_href,
            strict_type => 1,
        },
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

    use MIP::Recipes::Analysis::Bwa_mem qw{ analysis_bwa_mem analysis_run_bwa_mem };

    Readonly my $GENOME_BUILD_VERSION_GRCH_PRIOR_ALTS => 37;
    Readonly my $GENOME_BUILD_VERSION_HG_PRIOR_ALTS   => 19;

    if ( $human_genome_reference_source eq q{GRCh} ) {

        # Human genome version > GRCh37
        if ( $human_genome_reference_version > $GENOME_BUILD_VERSION_GRCH_PRIOR_ALTS ) {

            # Use bwa run mem recipe
            $analysis_recipe_href->{bwa_mem} = \&analysis_run_bwa_mem;
            return;
        }

        # Human genome version <= GRCh37
        # Use bwa mem recipe
        $analysis_recipe_href->{bwa_mem} = \&analysis_bwa_mem;
        return;
    }

    # hgXX build

    # Human genome version > hg19
    if ( $human_genome_reference_version > $GENOME_BUILD_VERSION_HG_PRIOR_ALTS ) {

        # Use bwa run mem recipe
        $analysis_recipe_href->{bwa_mem} = \&analysis_run_bwa_mem;
        return;
    }

    ## Human genome version <= hg19
    # Use bwa mem recipe
    $analysis_recipe_href->{bwa_mem} = \&analysis_bwa_mem;

    return;
}

sub set_recipe_gatk_variantrecalibration {

## Function : Update which gatk variant recalibration to use depending on number of samples
## Returns  :
## Arguments: $analysis_recipe_href => Analysis recipe hash {REF}
##          : $log                  => Log object to write to
##          : $sample_ids_ref       => Sample ids {REF}
##          : $use_cnnscorevariants => Use gatk cnnscorevariants recipe

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_recipe_href;
    my $log;
    my $sample_ids_ref;
    my $use_cnnscorevariants;

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
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
        use_cnnscorevariants => {
            store       => \$use_cnnscorevariants,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Recipes::Analysis::Gatk_cnnscorevariants
      qw{ analysis_gatk_cnnscorevariants };

    ## Use already set gatk_variantrecalibration recipe
    return if ( @{$sample_ids_ref} != 1 );

    return if ( not $use_cnnscorevariants );

    $log->warn(
q{Switched from VariantRecalibration to CNNScoreVariants for single sample analysis}
    );

    ## Use new CNN recipe for single samples
    $analysis_recipe_href->{gatk_variantrecalibration} = \&analysis_gatk_cnnscorevariants;

    return;
}

sub set_rankvariants_ar {

## Function : Update which rankvariants recipe to use
## Returns  :
## Arguments: $analysis_recipe_href => Analysis recipe hash {REF}
##          : $log                  => Log object to write to
##          : $parameter_href       => Parameter hash {REF}
##          : $sample_ids_ref       => Sample ids {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_recipe_href;
    my $log;
    my $parameter_href;
    my $sample_ids_ref;

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
