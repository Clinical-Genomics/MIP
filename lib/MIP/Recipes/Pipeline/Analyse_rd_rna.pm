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
    our $VERSION = 1.37;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ parse_rd_rna pipeline_analyse_rd_rna };
}

sub parse_rd_rna {

## Function : Rare disease RNA pipeline specific checks and parsing
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $broadcasts_ref        => Holds the parameters info for broadcasting later {REF}
##          : $file_info_href        => File info hash {REF}
##          : $order_parameters_ref  => Order of parameters (for structured output) {REF}
##          : $parameter_href        => Parameter hash {REF}
##          : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $broadcasts_ref;
    my $file_info_href;
    my $order_parameters_ref;
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
        order_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_parameters_ref,
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

    use MIP::Active_parameter qw{
      check_sample_id_in_hash_parameter
      check_sample_id_in_hash_parameter_path
      parse_infiles
      write_references
    };
    use MIP::Analysis qw{ broadcast_parameters check_ids_in_dna_vcf
      update_recipe_mode_for_fastq_compatibility
      update_recipe_mode_for_pedigree };
    use MIP::Config qw{ write_mip_config };
    use MIP::Contigs qw{ update_contigs_for_run };
    use MIP::Fastq qw{ parse_fastq_infiles };
    use MIP::File_info qw{ check_parameter_metafiles };
    use MIP::Parameter qw{ get_cache };
    use MIP::Sample_info qw{ set_parameter_in_sample_info };
    use MIP::Analysis qw{ set_ase_chain_recipes };
    use MIP::Star qw{ check_interleaved_files_for_star };

    ## Constants
    Readonly my @REMOVE_CONFIG_KEYS => qw{ associated_recipe };

    my $consensus_analysis_type = get_cache(
        {
            parameter_href => $parameter_href,
            parameter_name => q{consensus_analysis_type},
        }
    );

    ## Checks parameter metafile exists and set build_file parameter
    check_parameter_metafiles(
        {
            active_parameter_href => $active_parameter_href,
            file_info_href        => $file_info_href,
            parameter_href        => $parameter_href,
        }
    );

    ## Check sample_id provided in hash parameter is included in the analysis
    check_sample_id_in_hash_parameter(
        {
            active_parameter_href => $active_parameter_href,
            parameter_names_ref   => [qw{ analysis_type }],
            parameter_href        => $parameter_href,
            sample_ids_ref        => \@{ $active_parameter_href->{sample_ids} },
        }
    );

    ## Check sample_id provided in hash path parameter is included in the analysis and only represented once
    check_sample_id_in_hash_parameter_path(
        {
            active_parameter_href => $active_parameter_href,
            parameter_names_ref   => [qw{ infile_dirs }],
            sample_ids_ref        => \@{ $active_parameter_href->{sample_ids} },
        }
    );

    ## Check dna vcf
    check_ids_in_dna_vcf(
        {
            active_parameter_href => $active_parameter_href,
            dna_vcf_file          => $active_parameter_href->{dna_vcf_file},
            sample_info_href      => $sample_info_href,
        }
    );

    ## Set ASE recipes depending on previous check
    set_ase_chain_recipes(
        {
            active_parameter_href => $active_parameter_href,
        }
    );

    broadcast_parameters(
        {
            active_parameter_href => $active_parameter_href,
            broadcasts_ref        => $broadcasts_ref,
            order_parameters_ref  => $order_parameters_ref,
        }
    );

    ## Write references for this analysis to yaml
    write_references(
        {
            active_parameter_href => $active_parameter_href,
            outfile_path          => $active_parameter_href->{reference_info_file},
            parameter_href        => $parameter_href,
        }
    );

    ## Update recipes depending on pedigree
    update_recipe_mode_for_pedigree(
        {
            active_parameter_href => $active_parameter_href,
            recipes_ref           => [qw{ blobfish }],
            sample_info_href      => $sample_info_href,
        }
    );

    ## Write config file for case
    write_mip_config(
        {
            active_parameter_href => $active_parameter_href,
            remove_keys_ref       => \@REMOVE_CONFIG_KEYS,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Update contigs depending on settings in run (wes or if only male samples)
    update_contigs_for_run(
        {
            consensus_analysis_type => $consensus_analysis_type,
            exclude_contigs_ref     => \@{ $active_parameter_href->{exclude_contigs} },
            file_info_href          => $file_info_href,
            include_y               => $active_parameter_href->{include_y},
        }
    );

    ## Get the ".fastq(.gz)" files from the supplied infiles directory. Checks if the files exist
    parse_infiles(
        {
            active_parameter_href => $active_parameter_href,
            file_info_href        => $file_info_href,
        }
    );

    ## Reformat file names to MIP format, get file name info and add info to sample_info
    parse_fastq_infiles(
        {
            active_parameter_href => $active_parameter_href,
            file_info_href        => $file_info_href,
            sample_info_href      => $sample_info_href,
        }
    );

    check_interleaved_files_for_star(
        {
            file_info_href => $file_info_href,
            sample_ids_ref => $active_parameter_href->{sample_ids},
        }
    );

    ## Add to sample info
    set_parameter_in_sample_info(
        {
            active_parameter_href => $active_parameter_href,
            file_info_href        => $file_info_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Check recipe compability with fastq files
    my @recipes_to_check = qw{ arriba_ar salmon_quant };

  RECIPE:
    foreach my $recipe (@recipes_to_check) {

        update_recipe_mode_for_fastq_compatibility(
            {
                active_parameter_href => $active_parameter_href,
                file_info_href        => $file_info_href,
                parameter_href        => $parameter_href,
                recipe_name           => $recipe,
            }
        );
    }

    return;
}

sub pipeline_analyse_rd_rna {

## Function : Pipeline recipe for rare disease rna data analysis.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $broadcasts_ref        => Holds the parameters info for broadcasting later {REF}
##          : $file_info_href        => File info hash {REF}
##          : $job_id_href           => Job id hash {REF}
##          : $log                   => Log object to write to
##          : $order_parameters_ref  => Order of parameters (for structured output) {REF}
##          : $order_recipes_ref     => Execution order of recipes {REF}
##          : $parameter_href        => Parameter hash {REF}
##          : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $broadcasts_ref;
    my $file_info_href;
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

    use MIP::Constants qw{ set_singularity_constants };

    ## Recipes
    use MIP::Log::MIP_log4perl qw{ log_display_recipe_for_user };
    use MIP::Recipes::Analysis::Analysisrunstatus qw{ analysis_analysisrunstatus };
    use MIP::Recipes::Analysis::Arriba qw{ analysis_arriba };
    use MIP::Recipes::Analysis::Bcftools_merge qw{ analysis_bcftools_merge };
    use MIP::Recipes::Analysis::Blobfish qw{ analysis_blobfish };
    use MIP::Recipes::Analysis::BootstrapAnn qw{ analysis_bootstrapann };
    use MIP::Recipes::Analysis::Fastqc qw{ analysis_fastqc };
    use MIP::Recipes::Analysis::Gatk_asereadcounter qw{ analysis_gatk_asereadcounter };
    use MIP::Recipes::Analysis::Gatk_baserecalibration
      qw{ analysis_gatk_baserecalibration_rna };
    use MIP::Recipes::Analysis::Gatk_haplotypecaller qw{ analysis_gatk_haplotypecaller };
    use MIP::Recipes::Analysis::Gatk_splitncigarreads
      qw{ analysis_gatk_splitncigarreads };
    use MIP::Recipes::Analysis::Gatk_variantfiltration
      qw{ analysis_gatk_variantfiltration };
    use MIP::Recipes::Analysis::Genebody_coverage qw{ analysis_genebody_coverage };
    use MIP::Recipes::Analysis::Gffcompare qw{ analysis_gffcompare };
    use MIP::Recipes::Analysis::Gzip_fastq qw{ analysis_gzip_fastq };
    use MIP::Recipes::Analysis::Markduplicates qw{ analysis_markduplicates_rna };
    use MIP::Recipes::Analysis::Mip_qccollect qw{ analysis_mip_qccollect };
    use MIP::Recipes::Analysis::Mip_vercollect qw{ analysis_mip_vercollect };
    use MIP::Recipes::Analysis::Multiqc qw{ analysis_multiqc };
    use MIP::Recipes::Analysis::Picardtools_collectrnaseqmetrics
      qw{ analysis_picardtools_collectrnaseqmetrics };
    use MIP::Recipes::Analysis::Picardtools_mergesamfiles
      qw{ analysis_picardtools_mergesamfiles };
    use MIP::Recipes::Analysis::Preseq qw{ analysis_preseq };
    use MIP::Recipes::Analysis::Rseqc qw{ analysis_rseqc };
    use MIP::Recipes::Analysis::Sacct qw{ analysis_sacct };
    use MIP::Recipes::Analysis::Salmon_quant qw{ analysis_salmon_quant };
    use MIP::Recipes::Analysis::Star_fusion qw{ analysis_star_fusion };
    use MIP::Recipes::Analysis::Stringtie qw{ analysis_stringtie };
    use MIP::Recipes::Analysis::Trim_galore qw{ analysis_trim_galore };
    use MIP::Recipes::Analysis::Vcf_ase_reformat qw{ analysis_vcf_ase_reformat};
    use MIP::Recipes::Analysis::Vep qw{ analysis_vep };
    use MIP::Recipes::Build::Rd_rna qw{ build_rd_rna_meta_files };
    use MIP::Set::Analysis qw{ set_recipe_star_aln };

    ### Pipeline specific checks
    parse_rd_rna(
        {
            active_parameter_href => $active_parameter_href,
            broadcasts_ref        => $broadcasts_ref,
            file_info_href        => $file_info_href,
            order_parameters_ref  => $order_parameters_ref,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Set analysis constants
    set_singularity_constants( { active_parameter_href => $active_parameter_href, } );

    ### Build recipes
    $log->info(q{[Reference check - Reference prerequisites]});

    build_rd_rna_meta_files(
        {
            active_parameter_href => $active_parameter_href,
            file_info_href        => $file_info_href,
            job_id_href           => $job_id_href,
            log                   => $log,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Dispatch table
    my %analysis_recipe = (
        analysisrunstatus                => \&analysis_analysisrunstatus,
        arriba_ar                        => \&analysis_arriba,
        bcftools_merge                   => \&analysis_bcftools_merge,
        blobfish                         => \&analysis_blobfish,
        bootstrapann                     => \&analysis_bootstrapann,
        dna_vcf_reformat                 => \&analysis_vcf_ase_reformat,
        fastqc_ar                        => \&analysis_fastqc,
        gatk_asereadcounter              => \&analysis_gatk_asereadcounter,
        gatk_baserecalibration           => \&analysis_gatk_baserecalibration_rna,
        gatk_haplotypecaller             => \&analysis_gatk_haplotypecaller,
        gatk_splitncigarreads            => \&analysis_gatk_splitncigarreads,
        gatk_variantfiltration           => \&analysis_gatk_variantfiltration,
        genebody_coverage                => \&analysis_genebody_coverage,
        gffcompare_ar                    => \&analysis_gffcompare,
        markduplicates                   => \&analysis_markduplicates_rna,
        multiqc_ar                       => \&analysis_multiqc,
        picardtools_collectrnaseqmetrics => \&analysis_picardtools_collectrnaseqmetrics,
        picardtools_mergesamfiles        => \&analysis_picardtools_mergesamfiles,
        preseq_ar                        => \&analysis_preseq,
        qccollect_ar                     => \&analysis_mip_qccollect,
        rseqc                            => \&analysis_rseqc,
        sacct                            => \&analysis_sacct,
        salmon_quant                     => \&analysis_salmon_quant,
        star_aln                         => undef,
        star_fusion                      => \&analysis_star_fusion,
        stringtie_ar                     => \&analysis_stringtie,
        trim_galore_ar                   => \&analysis_trim_galore,
        varianteffectpredictor           => \&analysis_vep,
        version_collect_ar               => \&analysis_mip_vercollect,
    );

    ## Update which star recipe to use depending on fastq infile mix
    set_recipe_star_aln(
        {
            analysis_recipe_href => \%analysis_recipe,
            file_info_href       => $file_info_href,
            sample_ids_ref       => $active_parameter_href->{sample_ids},
        }
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
                        active_parameter_href => $active_parameter_href,
                        file_info_href        => $file_info_href,
                        job_id_href           => $job_id_href,
                        parameter_href        => $parameter_href,
                        recipe_name           => $recipe,
                        sample_id             => $sample_id,
                        sample_info_href      => $sample_info_href,
                    }
                );
            }
        }
        ## Family mode
        elsif ( $parameter_href->{$recipe}{analysis_mode} eq q{case} ) {

            $analysis_recipe{$recipe}->(
                {
                    active_parameter_href => $active_parameter_href,
                    file_info_href        => $file_info_href,
                    job_id_href           => $job_id_href,
                    parameter_href        => $parameter_href,
                    recipe_name           => $recipe,
                    sample_info_href      => $sample_info_href,
                }
            );
        }
    }

    return;
}

1;
