package MIP::Cli::Mip::Analyse::Rd_rna;

use 5.026;
use Carp;
use open qw{ :encoding(UTF-8) :std };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use MooseX::App::Command;
use MooseX::Types::Moose qw{ Str Int HashRef Num Bool ArrayRef };
use Moose::Util::TypeConstraints;

## MIPs lib
use MIP::Main::Analyse qw{ mip_analyse };

extends(qw{ MIP::Cli::Mip::Analyse });

command_short_description(q{Rare disease RNA analysis});

command_long_description(q{Rare diseae RNA analysis on wts sequence data});

command_usage(q{mip <analyse> <rd_rna> <case_id> --config <config_file> });

## Define, check and get Cli supplied parameters
_build_usage();

sub run {
    my ($arg_href) = @_;

    ## Remove Moose::App extra variable
    delete $arg_href->{extra_argv};

    ## Input from Cli
    my %active_parameter = %{$arg_href};

    use MIP::Definition qw{ get_dependency_tree_from_definition_file
      get_first_level_keys_order_from_definition_file
      get_parameter_definition_file_paths
      get_parameter_from_definition_files };
    use MIP::Dependency_tree qw{ get_dependency_tree_chain set_dependency_tree_order };
    use MIP::Parameter qw{ get_cache get_order_of_parameters print_recipe };

    ## %parameter holds all defined parameters for MIP analyse rd_rna
    ## CLI commands inheritance
    my $level     = q{rd_rna};
    my %parameter = get_parameter_from_definition_files( { level => $level, } );

    my @rd_rna_definition_file_paths = get_parameter_definition_file_paths( { level => $level, } );

    ### To write parameters and their values to log in logical order
    ## Adds the order of first level keys from definition files to array
    my @order_parameters =
      get_order_of_parameters( { define_parameters_files_ref => \@rd_rna_definition_file_paths, } );

    ## Print recipes if requested and exit
    print_recipe(
        {
            order_parameters_ref => \@order_parameters,
            parameter_href       => \%parameter,
            print_recipe         => $active_parameter{print_recipe},
            print_recipe_mode    => $active_parameter{print_recipe_mode},
        }
    );

    ## Get dependency tree and store in parameter hash
    %{ $parameter{dependency_tree_href} } =
      get_dependency_tree_from_definition_file( { level => $level, } );

## Sets chain id to parameters hash from the dependency tree
    get_dependency_tree_chain(
        {
            dependency_tree_href => $parameter{dependency_tree_href},
            parameter_href       => \%parameter,
        }
    );

    ## Set order of recipes according to dependency tree
    set_dependency_tree_order(
        {
            dependency_tree_href => $parameter{dependency_tree_href},
            recipes_ref          => \@{ $parameter{cache}{order_recipes_ref} },
        }
    );

    ## File info hash
    my %file_info = (
        human_genome_reference_file_endings => [qw{ .dict .fai }],
        salmon_quant_reference_genome       => [qw{ _salmon_quant_genome_dir }],
        star_aln_reference_genome           => [qw{ _star_genome_dir }],
        star_fusion_reference_genome        => [qw{ _star_fusion_genome_dir }],
        transcript_annotation_file_endings  => [qw{ .refFlat .rrna.interval_list }],
    );

    mip_analyse(
        {
            active_parameter_href => \%active_parameter,
            file_info_href        => \%file_info,
            order_parameters_ref  => \@order_parameters,
            parameter_href        => \%parameter,
        }
    );

    return;
}

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

    option(
        q{arriba_ar} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Detect and visualize fusions using Arriba},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{arriba_blacklist_path} => (
            cmd_tags      => [q{Recipe argument}],
            documentation => q{Path to arriba blacklist file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{arriba_cytoband_path} => (
            cmd_tags      => [q{Recipe argument}],
            documentation => q{Path to arriba cytoband file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{arriba_proteindomain_path} => (
            cmd_tags      => [q{Recipe argument}],
            documentation => q{Path to arriba protein domain file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{arriba_use_sample_id_as_display_name} => (
            cmd_tags      => [q{Default: 0}],
            documentation => q{Use sample id as display name for arriba fusion report},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{bcftools_merge} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Merge vcfs before annotation},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{blobfish} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Run BlobFish on salmon quant files},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{bootstrapann} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Run BootstrapAnn on ASE file},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{dna_vcf_file} => (
            cmd_flag      => q{dna_vcf_file},
            cmd_tags      => [q{Format: vcf | bcf}],
            documentation => q{Variantcalls made on wgs or wes data },
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{force_dna_ase} => (
            documentation => q{Force ASE analysis on partially matching dna-rna samples},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_disable_auto_index_and_file_lock} => (
            cmd_flag      => q{gatk_dis_auto_ind_fl},
            documentation => q{Disable auto index creation and locking when reading rods},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{genebody_coverage} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Run geneBody_coverage2.py on bam},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{infile_dirs} => (
            cmd_aliases   => [qw{ ifd }],
            cmd_tags      => [q{infile_dirs=sample_id}],
            documentation => q{Infile directory(s)},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{picardtools_path} => (
            documentation => q{Path to Picardtools},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{time_point} => (
            cmd_tags      => [q{sample_id=time_point}],
            documentation => q{Time point of replicate},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{split_fastq_file} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Split fastq files in batches of X reads and exits},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{split_fastq_file_read_batch} => (
            cmd_flag      => q{spt_fsq_rd_bt},
            cmd_tags      => [q{Default: 25,000,000}],
            documentation => q{Number of sequence reads to place in each batch},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{gzip_fastq} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Gzip fastq files},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{fastqc_ar} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Sequence quality analysis using FastQC},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{picardtools_mergesamfiles} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Merge bam files using Picardtools},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{picardtools_collectrnaseqmetrics} => (
            cmd_flag      => q{ppt_col_rna_met},
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Collect RNA seq metrics},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_baserecalibration} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Recalibration of bases using GATK BaseReCalibrator/PrintReads},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_baserecalibration_covariates} => (
            cmd_flag => q{gatk_baserecal_covariates},
            cmd_tags => [
q{Default: ReadGroupCovariate, ContextCovariate, CycleCovariate, QualityScoreCovariate}
            ],
            documentation => q{GATK BaseReCalibration covariates},
            is            => q{rw},
            isa           => ArrayRef [
                enum(
                    [
                        qw{ ContextCovariate CycleCovariate QualityScoreCovariate ReadGroupCovariate RepeatLengthCovariate RepeatUnitCovariate RepeatUnitAndLengthCovariate }
                    ]
                )
            ],
        )
    );

    option(
        q{gatk_baserecalibration_disable_indel_qual} => (
            cmd_flag      => q{gatk_baserecal_dis_indel_q},
            documentation => q{Disable indel quality scores},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_baserecalibration_known_sites} => (
            cmd_flag => q{gatk_baserecal_ks},
            cmd_tags => [
q{Default: grch37_dbsnp_-138-.vcf, grch37_1000g_indels_-phase1-.vcf, grch37_mills_and_1000g_indels_-gold_standard-.vcf}
            ],
            documentation => q{GATK BaseReCalibration known SNV and INDEL sites},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{gatk_baserecalibration_read_filters} => (
            cmd_flag      => q{gatk_baserecal_read_filts},
            cmd_tags      => [q{Default: OverclippedRead}],
            documentation => q{Filter out reads according to set filter},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{gatk_baserecalibration_static_quantized_quals} => (
            cmd_flag      => q{gatk_baserecal_sta_qua_qua},
            cmd_tags      => [q{Default: 10,20,30,40}],
            documentation => q{Static binning of base quality scores},
            is            => q{rw},
            isa           => ArrayRef [Int],
        )
    );

    option(
        q{salmon_quant} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Quantify transcripts using salmon},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{star_aln} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Align reads using Star aln},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{star_ulimit_n} => (
            documentation => q{Set ulimit -n for star recipe},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{align_intron_max} => (
            cmd_tags      => [q{Default: 100,000}],
            documentation => q{Maximum intron size},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{align_mates_gap_max} => (
            cmd_tags      => [q{Default: 100,000}],
            documentation => q{Maximum gap between two mates},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{align_sjdb_overhang_min} => (
            cmd_tags      => [q{Default: 10}],
            documentation => q{Minimum overhang (i.e. block size) for spliced alignments},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{chim_junction_overhang_min} => (
            cmd_tags      => [q{Default: 12}],
            documentation => q{Minimum overhang for a chimeric junction},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{chim_out_type} => (
            cmd_tags      => [q{Default: WithinBam}],
            documentation => q{Type of chimeric output},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{chim_segment_min} => (
            cmd_tags      => [q{Default: 12}],
            documentation => q{Minimum length of chimeric segment},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{pe_overlap_nbases_min} => (
            cmd_tags      => [q{Default: 10}],
            documentation => q{Min bases overlap to merge reads whne alingning },
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{two_pass_mode} => (
            cmd_tags      => [q{Default: Basic}],
            documentation => q{Two pass mode setting},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{star_fusion} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Detect fusion transcripts with star fusion},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{star_fusion_min_junction_reads} => (
            documentation => q{STAR-Fusion: Minimum junction spanning reads},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{star_fusion_pfam_db} => (
            documentation => q{Pfam db, used when building Star-Fusion references},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{rseq} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Qc using rseqc},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{rseqc_transcripts_file} => (
            cmd_tags      => [q{Rseqc transcripts file: Format: bed}],
            documentation => q{Input for rseqc to build transcript bed format file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{gatk_haplotypecaller} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Variant discovery using GATK HaplotypeCaller},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_haplotypecaller_annotation} => (
            cmd_flag => q{gatk_haplotype_ann},
            cmd_tags => [
q{Default: BaseQualityRankSumTest, ChromosomeCounts, Coverage, DepthPerAlleleBySample, FisherStrand, MappingQualityRankSumTest, QualByDepth, RMSMappingQuality, ReadPosRankSumTest, StrandOddsRatio}
            ],
            documentation => q{GATK HaploTypeCaller annotations},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{gatk_haplotypecaller_emit_ref_confidence} => (
            cmd_flag      => q{gatk_haplotype_emit_ref_conf},
            cmd_tags      => [q{Default: NONE}],
            documentation => q{VCF to produce},
            is            => q{rw},
            isa           => ArrayRef [ enum( [qw{ NONE BP_RESOLUTION GVCF }] ), ],
        )
    );

    option(
        q{gatk_haplotypecaller_pcr_indel_model} => (
            cmd_flag      => q{gatk_haplotype_pcr_ind_mod},
            cmd_tags      => [q{Default: CONSERVATIVE; Set to "0" to disable}],
            documentation => q{PCR indel model to use},
            is            => q{rw},
            isa => ArrayRef [ enum( [ 0, qw{ AGGRESSIVE CONSERVATIVE HOSTILE NONE } ] ), ],
        )
    );

    option(
        q{gatk_haplotypecaller_snp_known_set} => (
            cmd_flag      => q{gatk_haplotype_snp_ks},
            cmd_tags      => [q{Default: grch37_dbsnp_-138-.vcf}],
            documentation => q{GATK HaplotypeCaller dbSNP set for annotating ID columns},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{markduplicates} => (
            cmd_flag      => q{markduplicates},
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Markduplicate reads},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{markduplicates_picardtools_opt_dup_dist} => (
            cmd_flag      => q{picard_mdup_odd},
            cmd_tags      => [q{Default: 2500}],
            documentation => q{Picardtools markduplicates optical duplicate distance},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{gatk_splitncigarrreads} => (
            cmd_flag      => q{gatk_splitncigarreads},
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Split reads that contain Ns in their cigar},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_haplotypecaller} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Variant discovery using GATK HaplotypeCaller},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_haplotypecaller_annotation} => (
            cmd_flag => q{gatk_haplotype_ann},
            cmd_tags => [
q{Default: BaseQualityRankSumTest, ChromosomeCounts, Coverage, DepthPerAlleleBySample, FisherStrand, MappingQualityRankSumTest, QualByDepth, RMSMappingQuality, ReadPosRankSumTest, StrandOddsRatio}
            ],
            documentation => q{GATK HaploTypeCaller annotations},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{gatk_haplotypecaller_no_soft_clipped_bases} => (
            cmd_flag      => q{gatk_haplotype_no_soft_cb},
            documentation => q{Do not include soft clipped bases in the variant calling},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_haplotypecaller_pcr_indel_model} => (
            cmd_flag      => q{gatk_haplotype_pcr_ind_mod},
            cmd_tags      => [q{Default: None; Set to "0" to disable}],
            documentation => q{PCR indel model to use},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_haplotypecaller_snp_known_set} => (
            cmd_flag      => q{gatk_haplotype_snp_ks},
            cmd_tags      => [q{Default: grch37_dbsnp_-138-.vcf}],
            documentation => q{GATK HaplotypeCaller dbSNP set for annotating ID columns},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{gatk_asereadcounter} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Allel specific expression},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_variantfiltration} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Hard filterering of variants},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_variantfiltration_cluster_size} => (
            cmd_flag      => q{gatk_variantfiltration_cluster_size},
            cmd_tags      => [q{Default: 3}],
            documentation => q{GATK VariantFiltration, the number of SNPs which make up a cluster},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{gatk_variantfiltration_filter} => (
            cmd_flag      => q{gatk_variantfiltration_filter},
            cmd_tags      => [q{filter_name=filter_expression}],
            documentation => q{GATK VariantFiltration, the filter to apply},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{gatk_variantfiltration_cluster_window_size} => (
            cmd_flag      => q{gatk_variantfiltration_cluster_window_size},
            cmd_tags      => [q{Default: 35}],
            documentation =>
              q{GATK VariantFiltration, window size (in bases) in which to evaluate clustered SNPs},
            is  => q{rw},
            isa => Int,
        )
    );

    option(
        q{gffcompare_ar} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Compare RNA transcripts to reference using GffCompare},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{qccollect_ar} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Collect QC metrics from recipes output},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{qccollect_eval_metric_file} => (
            cmd_tags      => [q{Default: qc_eval_metric_-v1.3-.yaml}],
            documentation => q{File containing Qccollect evaluation metrics},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{qccollect_regexp_file} => (
            cmd_tags      => [q{Default: qc_regexp_-v1.25-.yaml}],
            documentation =>
q{Regular expression file containing the regular expression to be used for each program},
            is  => q{rw},
            isa => Str,
        )
    );

    option(
        q{qccollect_sampleinfo_file} => (
            cmd_tags      => [q{Default: {outdata_dir}/{case_id}/{case_id}_qc_sample_info.yaml}],
            documentation =>
              q{Sample info file containing info on what to parse from this analysis run},
            is  => q{rw},
            isa => Str,
        )
    );

    option(
        q{qccollect_skip_evaluation} => (
            documentation => q{Skip evaluation step in qccollect},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{qccollect_store_metrics_outfile} => (
            cmd_tags      => [q{Default: {outdata_dir}/{case_id}_metrics_deliverables.yaml}],
            documentation => q{File containing metrics from this analysis run},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{multiqc_ar} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Create aggregate bioinformatics analysis report across many samples},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{multiqc_per_sample} => (
            documentation => q{Generate sample specific reports},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{preseq_ar} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Estimate library complexity},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{recipe_core_number} => (
            cmd_tags      => [q{recipe_name=X(cores)}],
            documentation => q{Set the number of cores for each recipe},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{recipe_memory} => (
            cmd_tags      => [q{recipe_name=X(G)}],
            documentation => q{Set the memory for each recipe},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{recipe_time} => (
            cmd_tags      => [q{recipe_name=time(hours)}],
            documentation => q{Set the time allocation for each recipe},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{stringtie_ar} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Assemble alignments using StringTie},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{stringtie_junction_reads} => (
            documentation => q{StringTie min junction spanning reads},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{stringtie_minimum_coverage} => (
            documentation => q{StringTie min transcript coverage},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{library_type} => (
            cmd_flag      => q{library_type},
            cmd_tags      => [q{Default: reverse_stranded}],
            documentation => q{Strandedness of library},
            is            => q{rw},
            isa => ArrayRef [ enum( [qw{ unstranded forward_stranded reverse_stranded }] ), ],
        )
    );

    option(
        q{transcript_annotation} => (
            cmd_tags      => [q{Transcripts file: Format: GTF}],
            documentation => q{Transcript file for the rd_rna pipeline},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{trim_galore_ar} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Trim fastq files using Trim galore},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{varianteffectpredictor} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Annotate variants using VEP},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{vep_custom_annotation} => (
            documentation => q{VEP custom annotation},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{vep_directory_cache} => (
            documentation => q{Specify the cache directory to use},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{vep_features} => (
            cmd_tags => [
q{Default: hgvs, symbol, numbers, sift, polyphen, humdiv, domains, protein, ccds, uniprot, biotype, regulatory, tsl, canonical, per_gene, appris}
            ],
            documentation => q{VEP features},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{vep_plugins_dir_path} => (
            documentation => q{Path to directory with VEP plugins},
            is            => q{rw},
            isa           => Str,
        )
    );

    return;
}

1;
