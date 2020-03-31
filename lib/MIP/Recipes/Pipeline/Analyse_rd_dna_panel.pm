package MIP::Recipes::Pipeline::Analyse_rd_dna_panel;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use List::MoreUtils qw { any };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $SPACE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.09;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ pipeline_analyse_rd_dna_panel };
}

sub pipeline_analyse_rd_dna_panel {

## Function : Pipeline recipe for dna panel data analysis.
## Returns  :
## Arguments: $active_parameter_href           => Active parameters for this analysis hash {REF}
##          : $broadcasts_ref                  => Holds the parameters info for broadcasting later {REF}
##          : $file_info_href                  => File info hash {REF}
##          : $infile_both_strands_prefix_href => The infile(s) without the ".ending" and strand info {REF}
##          : $infile_lane_prefix_href         => Infile(s) without the ".ending" {REF}
##          : $job_id_href                     => Job id hash {REF}
##          : $log                             => Log object to write to
##          : $order_parameters_ref            => Order of parameters (for structured output) {REF}
##          : $order_recipes_ref               => Order of recipes
##          : $parameter_href                  => Parameter hash {REF}
##          : $sample_info_href                => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $broadcasts_ref;
    my $file_info_href;
    my $infile_both_strands_prefix_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $log;
    my $order_parameters_ref;
    my $order_recipes_ref;
    my $parameter_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        broadcasts_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$broadcasts_ref,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        infile_both_strands_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_both_strands_prefix_href,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        order_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_parameters_ref,
            strict_type => 1,
        },
        order_recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_recipes_ref,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
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

    use MIP::Check::Pipeline qw{ check_rd_dna_panel };
    use MIP::Constants qw{ set_singularity_constants };
    use MIP::Log::MIP_log4perl qw{ log_display_recipe_for_user };
    use MIP::Parse::Reference qw{ parse_reference_for_vt };
    use MIP::Set::Analysis qw{ set_recipe_bwa_mem };

    ## Recipes
    use MIP::Recipes::Analysis::Analysisrunstatus qw{ analysis_analysisrunstatus };
    use MIP::Recipes::Analysis::Cadd qw{ analysis_cadd_panel };
    use MIP::Recipes::Analysis::Fastqc qw{ analysis_fastqc };
    use MIP::Recipes::Analysis::Frequency_annotation
      qw{ analysis_frequency_annotation_panel };
    use MIP::Recipes::Analysis::Frequency_filter qw{ analysis_frequency_filter_panel };
    use MIP::Recipes::Analysis::Gatk_baserecalibration
      qw{ analysis_gatk_baserecalibration_panel };
    use MIP::Recipes::Analysis::Gatk_combinevariantcallsets
      qw{ analysis_gatk_combinevariantcallsets };
    use MIP::Recipes::Analysis::Gatk_haplotypecaller
      qw{ analysis_gatk_haplotypecaller_panel};
    use MIP::Recipes::Analysis::Gzip_fastq qw{ analysis_gzip_fastq };
    use MIP::Recipes::Analysis::Gatk_gathervcfs qw{ analysis_gatk_gathervcfs };
    use MIP::Recipes::Analysis::Gatk_genotypegvcfs qw{ analysis_gatk_genotypegvcfs };
    use MIP::Recipes::Analysis::Gatk_variantevalall qw{ analysis_gatk_variantevalall };
    use MIP::Recipes::Analysis::Gatk_variantrecalibration
      qw{ analysis_gatk_variantrecalibration_wes };
    use MIP::Recipes::Analysis::Markduplicates qw{ analysis_markduplicates_panel };
    use MIP::Recipes::Analysis::Mip_vercollect qw{ analysis_mip_vercollect };
    use MIP::Recipes::Analysis::Multiqc qw{ analysis_multiqc };
    use MIP::Recipes::Analysis::Picardtools_collecthsmetrics
      qw{ analysis_picardtools_collecthsmetrics };
    use MIP::Recipes::Analysis::Picardtools_collectmultiplemetrics
      qw{ analysis_picardtools_collectmultiplemetrics };
    use MIP::Recipes::Analysis::Rtg_vcfeval qw{ analysis_rtg_vcfeval  };
    use MIP::Recipes::Analysis::Sambamba_depth qw{ analysis_sambamba_depth };
    use MIP::Recipes::Analysis::Samtools_merge qw{ analysis_samtools_merge_panel };
    use MIP::Recipes::Analysis::Variant_integrity qw{ analysis_variant_integrity };
    use MIP::Recipes::Analysis::Vep qw{ analysis_vep };
    use MIP::Recipes::Analysis::Vt qw{ analysis_vt_panel };
    use MIP::Recipes::Build::Rd_dna qw{ build_rd_dna_meta_files };

    ### Pipeline specific checks
    check_rd_dna_panel(
        {
            active_parameter_href           => $active_parameter_href,
            broadcasts_ref                  => $broadcasts_ref,
            file_info_href                  => $file_info_href,
            infile_both_strands_prefix_href => $infile_both_strands_prefix_href,
            infile_lane_prefix_href         => $infile_lane_prefix_href,
            order_parameters_ref            => $order_parameters_ref,
            parameter_href                  => $parameter_href,
            sample_info_href                => $sample_info_href,
        }
    );

    ## Set analysis constants
    set_singularity_constants( { active_parameter_href => $active_parameter_href, } );

    ### Build recipes
    $log->info(q{[Reference check - Reference prerequisites]});

    build_rd_dna_meta_files(
        {
            active_parameter_href   => $active_parameter_href,
            file_info_href          => $file_info_href,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            log                     => $log,
            parameter_href          => $parameter_href,
            sample_info_href        => $sample_info_href,
        }
    );

    ## Check if vt has processed references
    ## If not try to reprocesses them before launching recipes
    $log->info(q{[Reference check - Reference processed by VT]});
    parse_reference_for_vt(
        {
            active_parameter_href   => $active_parameter_href,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            log                     => $log,
            parameter_href          => $parameter_href,
        }
    );

    ### Analysis recipes
    ## Create code reference table for pipeline analysis recipes
    my %analysis_recipe = (
        analysisrunstatus            => \&analysis_analysisrunstatus,
        bwa_mem                      => undef,
        cadd_ar                      => \&analysis_cadd_panel,
        fastqc_ar                    => \&analysis_fastqc,
        frequency_annotation         => \&analysis_frequency_annotation_panel,
        frequency_filter             => \&analysis_frequency_filter_panel,
        gatk_baserecalibration       => \&analysis_gatk_baserecalibration_panel,
        gatk_combinevariantcallsets  => \&analysis_gatk_combinevariantcallsets,
        gatk_haplotypecaller         => \&analysis_gatk_haplotypecaller_panel,
        gatk_gathervcfs              => \&analysis_gatk_gathervcfs,
        gatk_genotypegvcfs           => \&analysis_gatk_genotypegvcfs,
        gatk_variantevalall          => \&analysis_gatk_variantevalall,
        gatk_variantrecalibration    => \&analysis_gatk_variantrecalibration_wes,
        gzip_fastq                   => \&analysis_gzip_fastq,
        markduplicates               => \&analysis_markduplicates_panel,
        multiqc_ar                   => \&analysis_multiqc,
        picardtools_collecthsmetrics => \&analysis_picardtools_collecthsmetrics,
        picardtools_collectmultiplemetrics =>
          \&analysis_picardtools_collectmultiplemetrics,
        rtg_vcfeval            => \&analysis_rtg_vcfeval,
        sambamba_depth         => \&analysis_sambamba_depth,
        samtools_merge         => \&analysis_samtools_merge_panel,
        variant_integrity_ar   => \&analysis_variant_integrity,
        varianteffectpredictor => \&analysis_vep,
        version_collect_ar     => \&analysis_mip_vercollect,
        vt_ar                  => \&analysis_vt_panel,
    );

    ## Set correct bwa_mem recipe depending on version and source of the human_genome_reference: Source (hg19 or grch)
    set_recipe_bwa_mem(
        {
            analysis_recipe_href => \%analysis_recipe,
            human_genome_reference_source =>
              $file_info_href->{human_genome_reference_source},
            human_genome_reference_version =>
              $file_info_href->{human_genome_reference_version},
        }
    );

  RECIPE:
    foreach my $recipe ( @{$order_recipes_ref} ) {

        ## Skip not active recipes
        next RECIPE if ( not $active_parameter_href->{$recipe} );

        ## Skip recipe if not part of dispatch table (such as gzip_fastq)
        next RECIPE if ( not $analysis_recipe{$recipe} );

        ### Analysis recipes
        ## For displaying
        log_display_recipe_for_user(
            {
                log    => $log,
                recipe => $recipe,
            }
        );
        ## Sample mode
        if ( $parameter_href->{$recipe}{analysis_mode} eq q{sample} ) {

          SAMPLE_ID:
            foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

                $analysis_recipe{$recipe}->(
                    {
                        active_parameter_href   => $active_parameter_href,
                        file_info_href          => $file_info_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        job_id_href             => $job_id_href,
                        parameter_href          => $parameter_href,
                        recipe_name             => $recipe,
                        sample_id               => $sample_id,
                        sample_info_href        => $sample_info_href,
                    }
                );
            }
        }

        ## Family mode
        elsif ( $parameter_href->{$recipe}{analysis_mode} eq q{case} ) {

            $analysis_recipe{$recipe}->(
                {
                    active_parameter_href   => $active_parameter_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    parameter_href          => $parameter_href,
                    recipe_name             => $recipe,
                    sample_info_href        => $sample_info_href,
                }
            );
        }

    }
    return;
}

1;
