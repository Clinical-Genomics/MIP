package MIP::Cli::Mip::Analyse::Rd_rna;

use 5.026;
use Carp;
use File::Spec::Functions qw{ catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use strict;
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

our $VERSION = 1.07;

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

    use MIP::File::Format::Parameter qw{ parse_definition_file  };
    use MIP::File::Format::Yaml qw{ load_yaml order_parameter_names };
    use MIP::Get::Analysis
      qw{ get_dependency_tree_chain get_dependency_tree_order print_recipe };

    ## Mip analyse rd_rna parameters
    ## CLI commands inheritance
    my @definition_files = (
        catfile( $Bin, qw{ definitions mip_parameters.yaml } ),
        catfile( $Bin, qw{ definitions analyse_parameters.yaml } ),
        catfile( $Bin, qw{ definitions rd_rna_parameters.yaml } ),
    );

    ## Non mandatory parameter definition keys to check
    my $non_mandatory_parameter_keys_path =
      catfile( $Bin, qw{ definitions non_mandatory_parameter_keys.yaml } );

    ## Mandatory parameter definition keys to check
    my $mandatory_parameter_keys_path =
      catfile( $Bin, qw{ definitions mandatory_parameter_keys.yaml } );

    ### %parameter holds all defined parameters for MIP
    ### mip analyse rd_rna
    my %parameter;
    foreach my $definition_file (@definition_files) {

        %parameter = (
            %parameter,
            parse_definition_file(
                {
                    define_parameters_path => $definition_file,
                    non_mandatory_parameter_keys_path =>
                      $non_mandatory_parameter_keys_path,
                    mandatory_parameter_keys_path => $mandatory_parameter_keys_path,
                }
            ),
        );
    }

    ## Print recipes if requested and exit
    print_recipe(
        {
            define_parameters_files_ref => \@definition_files,
            parameter_href              => \%parameter,
            print_recipe                => $active_parameter{print_recipe},
            print_recipe_mode           => $active_parameter{print_recipe_mode},
        }
    );

    my %dependency_tree = load_yaml(
        {
            yaml_file => catfile( $Bin, qw{ definitions rd_rna_initiation_map.yaml } ),
        }
    );

## Sets chain id to parameters hash from the dependency tree
    get_dependency_tree_chain(
        {
            dependency_tree_href => \%dependency_tree,
            parameter_href       => \%parameter,
        }
    );

## Order recipes - Parsed from initiation file
    get_dependency_tree_order(
        {
            dependency_tree_href => \%dependency_tree,
            recipes_ref          => \@{ $parameter{cache}{order_recipes_ref} },
        }
    );

### To write parameters and their values to log in logical order
### Actual order of parameters in definition parameters file(s) does not matter
## Adds the order of first level keys from yaml files to array
    my @order_parameters;
    foreach my $define_parameters_file (@definition_files) {

        push @order_parameters,
          order_parameter_names(
            {
                file_path => $define_parameters_file,
            }
          );
    }

## File info hash
    my %file_info = (

        fusion_filter_reference_genome      => [qw{ _fusion_filter_genome_dir }],
        human_genome_reference_file_endings => [qw{ .dict .fai }],
        salmon_quant_reference_genome       => [qw{ _salmon_quant_genome_dir }],
        star_aln_reference_genome           => [qw{ _star_genome_dir }],
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
        q{blobfish} => (
            cmd_aliases   => [qw{ blb }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Run BlobFish on salmon quant files},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{bootstrapann} => (
            cmd_aliases   => [qw{ ba }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Run BootstrapAnn on ASE file},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_bundle_download_version} => (
            cmd_aliases   => [qw{ gbdv }],
            cmd_tags      => [q{Default: 2.8}],
            documentation => q{GATK FTP bundle download version},
            is            => q{rw},
            isa           => Num,
        )
    );

    option(
        q{gatk_disable_auto_index_and_file_lock} => (
            cmd_aliases   => [qw{ gdai }],
            cmd_flag      => q{gatk_dis_auto_ind_fl},
            documentation => q{Disable auto index creation and locking when reading rods},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_downsample_to_coverage} => (
            cmd_aliases   => [qw{ gdco }],
            cmd_tags      => [q{Default: 1000}],
            documentation => q{Coverage to downsample to at any given locus},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{gatk_use_new_qual_calculator} => (
            cmd_aliases   => [qw{ gatknq }],
            cmd_flag      => q{gatk_new_qual},
            documentation => q{Use new qual calculator},
            is            => q{rw},
            isa           => Bool,
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
        q{java_use_large_pages} => (
            cmd_aliases   => [qw{ jul }],
            documentation => q{Use large page memory},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{recipe_core_number} => (
            cmd_aliases   => [qw{ rcn }],
            cmd_tags      => [q{recipe_name=X(cores)}],
            documentation => q{Set the number of cores for each recipe},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{recipe_time} => (
            cmd_aliases   => [qw{ rot }],
            cmd_tags      => [q{recipe_name=time(hours)}],
            documentation => q{Set the time allocation for each recipe},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{picardtools_path} => (
            cmd_aliases   => [qw{ ptp }],
            documentation => q{Path to Picardtools},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sample_origin} => (
            cmd_aliases   => [qw{ samo }],
            cmd_tags      => [q{sample_id=sample_origin}],
            documentation => q{Sample origin of replicate},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{time_point} => (
            cmd_aliases   => [qw{ timp }],
            cmd_tags      => [q{sample_id=time_point}],
            documentation => q{Time point of replicate},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{split_fastq_file} => (
            cmd_aliases   => [qw{ sfq }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Split fastq files in batches of X reads and exits},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{split_fastq_file_read_batch} => (
            cmd_aliases   => [qw{ sfqrdb }],
            cmd_flag      => q{spt_fsq_rd_bt},
            cmd_tags      => [q{Default: 25,000,000}],
            documentation => q{Number of sequence reads to place in each batch},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{gzip_fastq} => (
            cmd_aliases   => [qw{ gz }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Gzip fastq files},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{fastqc_ar} => (
            cmd_aliases   => [qw{ fqc }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Sequence quality analysis using FastQC},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{picardtools_mergesamfiles} => (
            cmd_aliases   => [qw{ pms }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Merge bam files using Picardtools},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_baserecalibration} => (
            cmd_aliases => [qw{ gbr }],
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
              q{Recalibration of bases using GATK BaseReCalibrator/PrintReads},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_baserecalibration_covariates} => (
            cmd_aliases => [qw{ gbrcov }],
            cmd_flag    => q{gatk_baserecal_covariates},
            cmd_tags    => [
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
            cmd_aliases   => [qw{ gbrdiq }],
            cmd_flag      => q{gatk_baserecal_dis_indel_q},
            documentation => q{Disable indel quality scores},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_baserecalibration_known_sites} => (
            cmd_aliases => [qw{ gbrkst }],
            cmd_flag    => q{gatk_baserecal_ks},
            cmd_tags    => [
q{Default: GRCh37_dbsnp_-138-.vcf, GRCh37_1000g_indels_-phase1-.vcf, GRCh37_mills_and_1000g_indels_-gold_standard-.vcf}
            ],
            documentation => q{GATK BaseReCalibration known SNV and INDEL sites},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{gatk_baserecalibration_read_filters} => (
            cmd_aliases   => [qw{ gbrrf }],
            cmd_flag      => q{gatk_baserecal_read_filts},
            cmd_tags      => [q{Default: OverclippedRead}],
            documentation => q{Filter out reads according to set filter},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{gatk_baserecalibration_static_quantized_quals} => (
            cmd_aliases   => [qw{ gbrsqq }],
            cmd_flag      => q{gatk_baserecal_sta_qua_qua},
            cmd_tags      => [q{Default: 10,20,30,40}],
            documentation => q{Static binning of base quality scores},
            is            => q{rw},
            isa           => ArrayRef [Int],
        )
    );

    option(
        q{salmon_quant} => (
            cmd_aliases   => [qw{ sqt }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Quantify transcripts using salmon},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{star_aln} => (
            cmd_aliases   => [qw{ stn }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Align reads using Star aln},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{align_intron_max} => (
            cmd_aliases   => [qw{ stn_aim }],
            cmd_tags      => [q{Default: 100,000}],
            documentation => q{Maximum intron size},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{align_mates_gap_max} => (
            cmd_aliases   => [qw{ stn_amg }],
            cmd_tags      => [q{Default: 100,000}],
            documentation => q{Maximum gap between two mates},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{align_sjdb_overhang_min} => (
            cmd_aliases   => [qw{ stn_asom }],
            cmd_tags      => [q{Default: 10}],
            documentation => q{Minimum overhang (i.e. block size) for spliced alignments},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{chim_junction_overhang_min} => (
            cmd_aliases   => [qw{ stn_cjom }],
            cmd_tags      => [q{Default: 12}],
            documentation => q{Minimum overhang for a chimeric junction},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{chim_segment_min} => (
            cmd_aliases   => [qw{ stn_csm }],
            cmd_tags      => [q{Default: 12}],
            documentation => q{Minimum length of chimaeric segment},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{two_pass_mode} => (
            cmd_aliases   => [qw{ stn_tpm }],
            cmd_tags      => [q{Default: Basic}],
            documentation => q{Two pass mode setting},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{star_fusion} => (
            cmd_aliases   => [qw{ stf }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Detect fusion transcripts with star fusion},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{rseq} => (
            cmd_aliases   => [qw{ rseq }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Qc using rseqc},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{rseqc_transcripts_file} => (
            cmd_aliases   => [qw{ rseqctf }],
            cmd_tags      => [q{Rseqc transcripts file: Format: GTF}],
            documentation => q{Input for rseqc to build transcript bed format file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{gatk_haplotypecaller} => (
            cmd_aliases   => [qw{ ghc }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Variant discovery using GATK HaplotypeCaller},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_haplotypecaller_annotation} => (
            cmd_aliases => [qw{ ghcann }],
            cmd_flag    => q{gatk_haplotype_ann},
            cmd_tags    => [
q{Default: BaseQualityRankSumTest, ChromosomeCounts, Coverage, DepthPerAlleleBySample, FisherStrand, MappingQualityRankSumTest, QualByDepth, RMSMappingQuality, ReadPosRankSumTest, StrandOddsRatio}
            ],
            documentation => q{GATK HaploTypeCaller annotations},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{gatk_haplotypecaller_emit_ref_confidence} => (
            cmd_aliases   => [qw{ ghcerc }],
            cmd_flag      => q{gatk_haplotype_emit_ref_conf},
            cmd_tags      => [q{Default: NONE}],
            documentation => q{VCF to produce},
            is            => q{rw},
            isa           => ArrayRef [ enum( [qw{ NONE BP_RESOLUTION GVCF }] ), ],
        )
    );

    option(
        q{gatk_haplotypecaller_pcr_indel_model} => (
            cmd_aliases   => [qw{ ghcpim }],
            cmd_flag      => q{gatk_haplotype_pcr_ind_mod},
            cmd_tags      => [q{Default: CONSERVATIVE; Set to "0" to disable}],
            documentation => q{PCR indel model to use},
            is            => q{rw},
            isa =>
              ArrayRef [ enum( [ 0, qw{ AGGRESSIVE CONSERVATIVE HOSTILE NONE } ] ), ],
        )
    );

    option(
        q{gatk_haplotypecaller_snp_known_set} => (
            cmd_aliases   => [qw{ ghckse }],
            cmd_flag      => q{gatk_haplotype_snp_ks},
            cmd_tags      => [q{Default: GRCh37_dbsnp_-138-.vcf}],
            documentation => q{GATK HaplotypeCaller dbSNP set for annotating ID columns},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{markduplicates} => (
            cmd_aliases   => [qw{ pmd }],
            cmd_flag      => q{markduplicates},
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Markduplicate reads},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{markduplicates_picardtools_markduplicates} => (
            cmd_aliases   => [qw{ mdpmd }],
            cmd_flag      => q{picard_markduplicates},
            documentation => q{Markduplicates using Picardtools markduplicates},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{markduplicates_sambamba_markdup} => (
            cmd_aliases   => [qw{ mdsmd }],
            cmd_flag      => q{sambamba_markdup},
            documentation => q{Markduplicates using Sambamba markduplicates},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{markduplicates_sambamba_markdup_hash_table_size} => (
            cmd_aliases   => [qw{ mdshts }],
            cmd_flag      => q{sba_mdup_hts},
            cmd_tags      => [q{Default: 262144}],
            documentation => q{Sambamba size of hash table for finding read pairs},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{markduplicates_sambamba_markdup_io_buffer_size} => (
            cmd_aliases => [qw{ mdsibs }],
            cmd_flag    => q{sba_mdup_ibs},
            cmd_tags    => [q{Default: 2048}],
            documentation =>
q{Sambamba size of the io buffer for reading and writing BAM during the second pass},
            is  => q{rw},
            isa => Int,
        )
    );

    option(
        q{markduplicates_sambamba_markdup_overflow_list_size} => (
            cmd_aliases   => [qw{ mdsols }],
            cmd_flag      => q{sba_mdup_ols},
            cmd_tags      => [q{Default: 200000}],
            documentation => q{Sambamba size of the overflow list},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{gatk_splitncigarrreads} => (
            cmd_aliases   => [qw{ gs }],
            cmd_flag      => q{gatk_splitncigarreads},
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Split reads that contain Ns in their cigar},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_haplotypecaller} => (
            cmd_aliases   => [qw{ ghc }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Variant discovery using GATK HaplotypeCaller},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_haplotypecaller_annotation} => (
            cmd_aliases => [qw{ ghcann }],
            cmd_flag    => q{gatk_haplotype_ann},
            cmd_tags    => [
q{Default: BaseQualityRankSumTest, ChromosomeCounts, Coverage, DepthPerAlleleBySample, FisherStrand, MappingQualityRankSumTest, QualByDepth, RMSMappingQuality, ReadPosRankSumTest, StrandOddsRatio}
            ],
            documentation => q{GATK HaploTypeCaller annotations},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{gatk_haplotypecaller_no_soft_clipped_bases} => (
            cmd_aliases   => [qw{ ghcscb }],
            cmd_flag      => q{gatk_haplotype_no_soft_cb},
            documentation => q{Do not include soft clipped bases in the variant calling},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_haplotypecaller_pcr_indel_model} => (
            cmd_aliases   => [qw{ ghcpim }],
            cmd_flag      => q{gatk_haplotype_pcr_ind_mod},
            cmd_tags      => [q{Default: None; Set to "0" to disable}],
            documentation => q{PCR indel model to use},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_haplotypecaller_snp_known_set} => (
            cmd_aliases   => [qw{ ghckse }],
            cmd_flag      => q{gatk_haplotype_snp_ks},
            cmd_tags      => [q{Default: GRCh37_dbsnp_-138-.vcf}],
            documentation => q{GATK HaplotypeCaller dbSNP set for annotating ID columns},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{gatk_asereadcounter} => (
            cmd_aliases   => [qw{ gae }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Allel specific expression},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_variantfiltration} => (
            cmd_aliases   => [qw{ gvf }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Hard filterering of variants},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_variantfiltration_cluster_size} => (
            cmd_aliases => [qw{ gvfc }],
            cmd_flag    => q{gatk_variantfiltration_cluster_size},
            cmd_tags    => [q{Default: 3}],
            documentation =>
              q{GATK VariantFiltration, the number of SNPs which make up a cluster},
            is  => q{rw},
            isa => Int,
        )
    );

    option(
        q{gatk_variantfiltration_filter} => (
            cmd_aliases   => [qw{ gvff }],
            cmd_flag      => q{gatk_variantfiltration_filter},
            cmd_tags      => [q{filter_name=filter_expression}],
            documentation => q{GATK VariantFiltration, the filter to apply},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{gatk_variantfiltration_cluster_window_size} => (
            cmd_aliases => [qw{ gvfw }],
            cmd_flag    => q{gatk_variantfiltration_cluster_window_size},
            cmd_tags    => [q{Default: 35}],
            documentation =>
q{GATK VariantFiltration, window size (in bases) in which to evaluate clustered SNPs},
            is  => q{rw},
            isa => Int,
        )
    );

    option(
        q{gffcompare} => (
            cmd_aliases   => [qw{ gffcmp }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Compare RNA transcripts to reference using GffCompare},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{multiqc_ar} => (
            cmd_aliases => [qw{ mqc }],
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
              q{Create aggregate bioinformatics analysis report across many samples},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{multiqc_per_sample} => (
            cmd_aliases   => [qw{ mqcps }],
            documentation => q{Generate sample specific reports},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{sacct} => (
            cmd_aliases => [qw{ sac }],
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
              q{Generating sbatch script for SLURM info on each submitted job},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{sacct_format_fields} => (
            cmd_aliases => [qw{ sacfrf }],
            cmd_tags    => [
q{Default: jobid, jobname%50, account, partition, alloccpus, TotalCPU, elapsed, start, end, state, exitcode}
            ],
            documentation => q{Format and fields of sacct output},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{stringtie} => (
            cmd_aliases   => [qw{ strg }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Assemble alignments using StringTie},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{library_type} => (
            cmd_aliases   => [qw{ lt }],
            cmd_flag      => q{library_type},
            cmd_tags      => [q{Default: reverse_stranded}],
            documentation => q{Strandedness of library},
            is            => q{rw},
            isa =>
              ArrayRef [ enum( [qw{ unstranded forward_stranded reverse_stranded }] ), ],
        )
    );

    option(
        q{transcript_annotation} => (
            cmd_aliases   => [qw{ ttf }],
            cmd_tags      => [q{Transcripts file: Format: GTF}],
            documentation => q{Transcript file for the rd_rna pipeline},
            is            => q{rw},
            isa           => Str,
        )
    );

    return;
}

1;
