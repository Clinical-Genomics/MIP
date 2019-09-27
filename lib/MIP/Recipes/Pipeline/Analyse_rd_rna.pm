package MIP::Recipes::Pipeline::Analyse_rd_rna;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $CLOSE_BRACKET $OPEN_BRACKET $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.23;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ pipeline_analyse_rd_rna };
}

sub pipeline_analyse_rd_rna {

## Function : Pipeline recipe for rare disease rna data analysis.
## Returns  :
## Arguments: $active_parameter_href           => Active parameters for this analysis hash {REF}
##          : $broadcasts_ref                  => Holds the parameters info for broadcasting later {REF}
##          : $file_info_href                  => File info hash {REF}
##          : $infile_both_strands_prefix_href => The infile(s) without the ".ending" and strand info {REF}
##          : $infile_lane_prefix_href         => Infile(s) without the ".ending" {REF}
##          : $job_id_href                     => Job id hash {REF}
##          : $log                             => Log object to write to
##          : $order_parameters_ref            => Order of parameters (for structured output) {REF}
##          : $order_recipes_ref               => Execution order of recipes {REF}
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

    use MIP::Check::Pipeline qw{ check_rd_rna };

    ## Recipes
    use MIP::Log::MIP_log4perl qw{ log_display_recipe_for_user };
    use MIP::Recipes::Analysis::Analysisrunstatus qw{ analysis_analysisrunstatus };
    use MIP::Recipes::Analysis::Bcftools_merge qw{ analysis_bcftools_merge };
    use MIP::Recipes::Analysis::Blobfish qw{ analysis_blobfish };
    use MIP::Recipes::Analysis::BootstrapAnn qw{ analysis_bootstrapann };
    use MIP::Recipes::Analysis::Fastqc qw{ analysis_fastqc };
    use MIP::Recipes::Analysis::Gatk_asereadcounter qw{ analysis_gatk_asereadcounter };
    use MIP::Recipes::Analysis::Gatk_baserecalibration
      qw{ analysis_gatk_baserecalibration };
    use MIP::Recipes::Analysis::Gatk_haplotypecaller qw{ analysis_gatk_haplotypecaller };
    use MIP::Recipes::Analysis::Gatk_splitncigarreads
      qw{ analysis_gatk_splitncigarreads };
    use MIP::Recipes::Analysis::Gatk_variantfiltration
      qw{ analysis_gatk_variantfiltration };
    use MIP::Recipes::Analysis::Genebody_coverage qw{ analysis_genebody_coverage };
    use MIP::Recipes::Analysis::Gffcompare qw{ analysis_gffcompare };
    use MIP::Recipes::Analysis::Gzip_fastq qw{ analysis_gzip_fastq };
    use MIP::Recipes::Analysis::Markduplicates qw{ analysis_markduplicates };
    use MIP::Recipes::Analysis::Mip_vercollect qw{ analysis_mip_vercollect };
    use MIP::Recipes::Analysis::Multiqc qw{ analysis_multiqc };
    use MIP::Recipes::Analysis::Picardtools_mergesamfiles
      qw{ analysis_picardtools_mergesamfiles };
    use MIP::Recipes::Analysis::Preseq qw{ analysis_preseq };
    use MIP::Recipes::Analysis::Rseqc qw{ analysis_rseqc };
    use MIP::Recipes::Analysis::Sacct qw{ analysis_sacct };
    use MIP::Recipes::Analysis::Salmon_quant qw{ analysis_salmon_quant };
    use MIP::Recipes::Analysis::Star_aln qw{ analysis_star_aln };
    use MIP::Recipes::Analysis::Star_fusion qw{ analysis_star_fusion };
    use MIP::Recipes::Analysis::Stringtie qw{ analysis_stringtie };
    use MIP::Recipes::Analysis::Trim_galore qw{ analysis_trim_galore };
    use MIP::Recipes::Analysis::Vep qw{ analysis_vep_rna };
    use MIP::Recipes::Build::Rd_rna qw{build_rd_rna_meta_files};

    ### Pipeline specific checks
    check_rd_rna(
        {
            active_parameter_href           => $active_parameter_href,
            broadcasts_ref                  => $broadcasts_ref,
            file_info_href                  => $file_info_href,
            infile_both_strands_prefix_href => $infile_both_strands_prefix_href,
            infile_lane_prefix_href         => $infile_lane_prefix_href,
            log                             => $log,
            order_parameters_ref            => $order_parameters_ref,
            parameter_href                  => $parameter_href,
            sample_info_href                => $sample_info_href,
        }
    );

    ### Build recipes
    $log->info(q{[Reference check - Reference prerequisites]});

    build_rd_rna_meta_files(
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

    ## Dispatch table
    my %analysis_recipe = (
        analysisrunstatus         => \&analysis_analysisrunstatus,
        bcftools_merge            => \&analysis_bcftools_merge,
        blobfish                  => \&analysis_blobfish,
        bootstrapann              => \&analysis_bootstrapann,
        fastqc_ar                 => \&analysis_fastqc,
        gatk_asereadcounter       => \&analysis_gatk_asereadcounter,
        gatk_baserecalibration    => \&analysis_gatk_baserecalibration,
        gatk_haplotypecaller      => \&analysis_gatk_haplotypecaller,
        gatk_splitncigarreads     => \&analysis_gatk_splitncigarreads,
        gatk_variantfiltration    => \&analysis_gatk_variantfiltration,
        genebody_coverage         => \&analysis_genebody_coverage,
        gffcompare_ar             => \&analysis_gffcompare,
        markduplicates            => \&analysis_markduplicates,
        multiqc_ar                => \&analysis_multiqc,
        picardtools_mergesamfiles => \&analysis_picardtools_mergesamfiles,
        preseq_ar                 => \&analysis_preseq,
        rseqc                     => \&analysis_rseqc,
        sacct                     => \&analysis_sacct,
        salmon_quant              => \&analysis_salmon_quant,
        star_aln                  => \&analysis_star_aln,
        star_fusion               => \&analysis_star_fusion,
        stringtie_ar              => \&analysis_stringtie,
        trim_galore_ar            => \&analysis_trim_galore,
        varianteffectpredictor    => \&analysis_vep_rna,
        version_collect_ar        => \&analysis_mip_vercollect,
    );

  RECIPE:
    foreach my $recipe ( @{$order_recipes_ref} ) {

        ## Skip not active recipes
        next RECIPE if ( not $active_parameter_href->{$recipe} );

        ## Skip recipe if not part of dispatch table (such as gzip fastq)
        next RECIPE if ( not $analysis_recipe{$recipe} );

        ## For displaying
        log_display_recipe_for_user(
            {
                log    => $log,
                recipe => $recipe,
            }
        );

        ## Sample mode
        if ( $parameter_href->{$recipe}{analysis_mode} eq q{sample} ) {

          SAMPLE:
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
