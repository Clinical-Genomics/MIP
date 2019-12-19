package MIP::Cli::Mip::Analyse::Rd_dna;

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
use List::MoreUtils qw { any };
use MooseX::App::Command;
use MooseX::Types::Moose qw{ ArrayRef Bool HashRef Int Num Str };
use Moose::Util::TypeConstraints;

## MIPs lib
use MIP::Main::Analyse qw{ mip_analyse };

our $VERSION = 1.37;

extends(qw{ MIP::Cli::Mip::Analyse });

command_short_description(q{Rare disease DNA analysis});

command_long_description(q{Rare disease DNA analysis on wes, wgs or mixed sequence data});

command_usage(q{mip <analyse> <rd_dna> <case_id> --config <config_file> });

## Define, check and get Cli supplied parameters
_build_usage();

sub run {
    my ($arg_href) = @_;

    ## Remove Moose::App extra variable
    delete $arg_href->{extra_argv};

    ## Input from Cli
    my %active_parameter = %{$arg_href};

    use MIP::Definition
      qw{ get_definition_file_paths get_parameter_hash_from_definition_files };
    use MIP::File::Format::Yaml qw{ load_yaml order_parameter_names };
    use MIP::Get::Analysis
      qw{ get_dependency_tree_chain get_dependency_tree_order print_recipe };

    ## %parameter holds all defined parameters for MIP analyse rd_dna
    ## CLI commands inheritance level
    my %parameter = get_parameter_hash_from_definition_files( { level => q{rd_dna}, } );

    my @rd_dna_definition_file_paths =
      get_definition_file_paths( { level => q{rd_dna}, } );

    ## Print recipes if requested and exit
    print_recipe(
        {
            define_parameters_files_ref => \@rd_dna_definition_file_paths,
            parameter_href              => \%parameter,
            print_recipe                => $active_parameter{print_recipe},
            print_recipe_mode           => $active_parameter{print_recipe_mode},
        }
    );

    ## Get dependency tree and store in parameter hash
    my %dependency_tree = load_yaml(
        {
            yaml_file => catfile( $Bin, qw{ definitions rd_dna_initiation_map.yaml } ),
        }
    );
    $parameter{dependency_tree} = \%dependency_tree;

    ## Sets chain id to parameters hash from the dependency tree
    get_dependency_tree_chain(
        {
            dependency_tree_href => $parameter{dependency_tree},
            parameter_href       => \%parameter,
        }
    );

    ## Order recipes - Parsed from initiation file
    get_dependency_tree_order(
        {
            dependency_tree_href => $parameter{dependency_tree},
            recipes_ref          => \@{ $parameter{cache}{order_recipes_ref} },
        }
    );

    ### To write parameters and their values to log in logical order
    ### Actual order of parameters in definition parameters file(s) does not matter
    ## Adds the order of first level keys from yaml files to array
    my @order_parameters;

  DEFINITION_FILE:
    foreach my $define_parameters_file (@rd_dna_definition_file_paths) {

        push @order_parameters,
          order_parameter_names(
            {
                file_path => $define_parameters_file,
            }
          );
    }

    ## File info hash
    my %file_info = (

        # BWA human genome reference file endings
        bwa_build_reference => [qw{ .bwt .ann .amb .pac .sa }],

        exome_target_bed => [qw{ .interval_list .pad100.interval_list }],

        # Human genome meta files
        human_genome_reference_file_endings => [qw{ .dict .fai }],

        # RTG human genome reference file endings
        rtg_vcfeval_reference_genome => [qw{ _sdf_dir }],
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
        q{decompose_normalize_references} => (
            cmd_aliases => [qw{ dnr }],
            cmd_flag    => q{dec_norm_ref},
            cmd_tags    => [
q{gatk_baserecalibration_known_sites, gatk_haplotypecaller_snp_known_set, gatk_variantrecalibration_resource_snv, gatk_variantrecalibration_resource_indel, frequency_genmod_filter_1000g, gatk_varianteval_gold, gatk_varianteval_dbsnp}
            ],
            documentation => q{Set the references to be decomposed and normalized},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{exome_target_bed} => (
            cmd_aliases => [qw{ extb }],
            cmd_tags => [q{file.bed=Sample_id; Default: latest_supported_capturekit.bed}],
            documentation => q{Exome target bed file per sample id},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{expected_coverage} => (
            cmd_aliases   => [qw{ ec }],
            cmd_tags      => [q{sample_id=expected_coverage}],
            documentation => q{Expected mean target coverage for analysis},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{frequency_annotation} => (
            cmd_aliases   => [qw{ fqa }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Annotate vcf with allele frequencies},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{fqa_vcfanno_config} => (
            cmd_aliases   => [qw{ fqavac }],
            documentation => q{Frequency vcfanno toml config},
            is            => q{rw},
            isa           => Str,
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
        q{gatk_use_new_qual_calculator} => (
            cmd_aliases   => [qw{ gatknq }],
            cmd_flag      => q{gatk_new_qual},
            documentation => q{Use new qual calculator},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{genomic_set} => (
            cmd_aliases   => [qw{ ges }],
            cmd_tags      => [q{sorted BED}],
            documentation => q{Selection of relevant regions post alignment},
            is            => q{ro},
            isa           => Str,
        )
    );

    option(
        q{human_genome_reference} => (
            cmd_aliases   => [qw{ hgr }],
            cmd_tags      => [q{Default: grch37_homo_sapiens_-d5-.fasta}],
            documentation => q{Human genome reference},
            is            => q{rw},
            isa           => Str,
        )
    );

    has(
        q{recipe_core_number} => (
            cmd_aliases   => [qw{ rcn }],
            cmd_tags      => [q{recipe_name=X(cores)}],
            documentation => q{Set the number of cores for each recipe},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{set_recipe_core_number} => (
            cmd_aliases   => [qw{ srcn }],
            cmd_tags      => [q{recipe_name=X(cores)}],
            documentation => q{Set the number of cores for specific recipe(s)},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    has(
        q{recipe_memory} => (
            cmd_aliases   => [qw{ rm }],
            cmd_tags      => [q{recipe_name=X(G)}],
            documentation => q{Set the memory for each recipe},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{set_recipe_memory} => (
            cmd_aliases   => [qw{ srm }],
            cmd_tags      => [q{recipe_name=X(G)}],
            documentation => q{Set the memory for specific recipe(s)},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    has(
        q{recipe_time} => (
            cmd_aliases   => [qw{ rot }],
            cmd_tags      => [q{recipe_name=time(hours)}],
            documentation => q{Set the time allocation for each recipe},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{set_recipe_time} => (
            cmd_aliases   => [qw{ srot }],
            cmd_tags      => [q{recipe_name=time(hours)}],
            documentation => q{Set the time allocation for specific recipe(s)},
            is            => q{rw},
            isa           => HashRef,
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
            cmd_aliases   => [qw{ ptp }],
            documentation => q{Path to Picardtools},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{replace_iupac} => (
            cmd_aliases   => [qw{ riu }],
            documentation => q{Replace IUPAC code in alternative alleles with N},
            is            => q{rw},
            isa           => Bool,
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
        q{bwa_mem} => (
            cmd_aliases   => [qw{ mem }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Align reads using Bwa Mem},
            is            => q{rw},
            isa           => enum( [ 1, 2 ] ),
        )
    );

    option(
        q{bwa_mem_bamstats} => (
            cmd_aliases   => [qw{ memsts }],
            documentation => q{Collect statistics from BAM files},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{bwa_mem_cram} => (
            cmd_aliases   => [qw{ memcrm }],
            documentation => q{Use CRAM-format for additional output file},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{bwa_mem_hla} => (
            cmd_aliases   => [qw{ memhla }],
            documentation => q{Apply HLA typing},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{bwa_sambamba_sort_memory_limit} => (
            cmd_aliases => [qw{ memssm }],
            cmd_flag    => q{bwa_sbm_srt_ml},
            cmd_tags    => [q{Default: 32G}],
            documentation =>
              q{Set the memory limit for Sambamba sort after bwa alignment},
            is  => q{rw},
            isa => Str,
        )
    );

    option(
        q{bwa_soft_clip_sup_align} => (
            cmd_aliases   => [qw{ memscsa }],
            documentation => q{Use soft clipping for supplementary alignments},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{picardtools_mergesamfiles} => (
            cmd_aliases => [qw{ ptm }],
            cmd_flag    => q{picardtools_mergesamfiles},
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
              q{Merge (BAM file(s) ) or rename single samples for downstream processing},
            is  => q{rw},
            isa => enum( [ 1, 2 ] ),
        )
    );

    option(
        q{markduplicates} => (
            cmd_aliases   => [qw{ md }],
            cmd_flag      => q{markduplicates},
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Markduplicate reads},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{markduplicates_no_bam_to_cram} => (
            cmd_aliases   => [qw{ mdnbtc }],
            cmd_flag      => q{markduplicates_nbtc},
            documentation => q{Generate CRAM from BAM},
            is            => q{rw},
            isa           => Bool,
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
        q{gatk_baserecalibration_no_bam_to_cram} => (
            cmd_aliases   => [qw{ gbrnbtc }],
            cmd_flag      => q{gatk_baserecal_nbtc},
            documentation => q{Generate CRAM from BAM},
            is            => q{rw},
            isa           => Bool,
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
q{Default: grch37_dbsnp_-138-.vcf, grch37_1000g_indels_-phase1-.vcf, grch37_mills_and_1000g_indels_-gold_standard-.vcf}
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
        q{chanjo_sexcheck} => (
            cmd_aliases   => [qw{ phs }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Predicts gender from sex chromosome coverage},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{chanjo_sexcheck_log_level} => (
            cmd_aliases   => [qw{ chslle }],
            cmd_flag      => q{chanjo_sexcheck_ll},
            documentation => q{Set chanjo sex log level},
            is            => q{rw},
            isa           => enum( [qw{ DEBUG INFO WARNING ERROR CRITICAL }] ),
        )
    );

    option(
        q{sambamba_depth} => (
            cmd_aliases   => [qw{ sdt }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Sambamba depth coverage analysis},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{sambamba_depth_bed} => (
            cmd_aliases   => [qw{ sdtbed }],
            documentation => q{Reference bed file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sambamba_depth_base_quality} => (
            cmd_aliases   => [qw{ sdtbaq }],
            cmd_flag      => q{sba_depth_bq},
            cmd_tags      => [q{Default: 10}],
            documentation => q{Do not count bases with lower base quality},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{sambamba_depth_cutoffs} => (
            cmd_aliases   => [qw{ sdtcut }],
            cmd_flag      => q{sba_depth_co},
            documentation => q{Read depth cutoff},
            is            => q{rw},
            isa           => ArrayRef [Int],
        )
    );

    option(
        q{sambamba_depth_mapping_quality} => (
            cmd_aliases   => [qw{ sdtmaq }],
            cmd_flag      => q{sba_depth_mq},
            cmd_tags      => [q{Default: 10}],
            documentation => q{Do not count reads with lower mapping quality},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{sambamba_depth_mode} => (
            cmd_aliases   => [qw{ sdtmod }],
            documentation => q{Mode unit to print the statistics on},
            is            => q{rw},
            isa           => enum( [qw{ base region window }] ),
        )
    );

    option(
        q{sambamba_depth_noduplicates} => (
            cmd_aliases   => [qw{ sdtndu }],
            cmd_flag      => q{sba_depth_nod},
            documentation => q{Do not include duplicates in coverage calculation},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{sambamba_depth_quality_control} => (
            cmd_aliases   => [qw{ sdtfqc }],
            cmd_flag      => q{sba_depth_qc},
            documentation => q{Do not include reads with failed quality control},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{tiddit_coverage} => (
            cmd_aliases   => [qw{ tcv }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Generate coverage data from alignment},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{picardtools_collectmultiplemetrics} => (
            cmd_aliases   => [qw{ ptcmm }],
            cmd_flag      => q{ppt_col_mul_met},
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Qc metrics calculation},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{picardtools_collecthsmetrics} => (
            cmd_aliases   => [qw{ ptchs }],
            cmd_flag      => q{ppt_col_hs_met},
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Qc metrics calculation for capture},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{cnvnator_ar} => (
            cmd_aliases   => [qw{ cnv }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Structural variant calling using CNVnator},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{cnv_bin_size} => (
            cmd_aliases   => [qw{ cnvhbs }],
            cmd_tags      => [q{Default: 1000}],
            documentation => q{CNVnator bin size},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{delly_call} => (
            cmd_aliases   => [qw{ delc }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Structural variant calling using Delly},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{delly_reformat} => (
            cmd_aliases   => [qw{ del }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Merge, regenotype and filter using Delly},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{delly_exclude_file} => (
            cmd_aliases   => [qw{ delexc }],
            cmd_tags      => [q{Default: hg19_human_excl_-0.7.6-.tsv}],
            documentation => q{Exclude centomere and telemore regions in delly calling},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{delly_types} => (
            cmd_aliases   => [qw{ deltyp }],
            cmd_tags      => [q{Default: DEL,DUP,INV,INS}],
            documentation => q{Type of SV to call},
            is            => q{rw},
            isa           => ArrayRef [ enum( [qw{ DEL DUP INV INS TRA }] ), ],
        )
    );

    option(
        q{expansionhunter} => (
            cmd_aliases   => [qw{ exp }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Anaylse expansions of Short Tandem Repeats},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{expansionhunter_variant_catalog_file_path} => (
            cmd_aliases   => [qw{ exphun_vcfp }],
            cmd_flag      => q{exphun_var_cat_fp},
            documentation => q{Path to variant catalog json file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{manta} => (
            cmd_aliases   => [qw{ mna }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Structural variant calling using Manta},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{manta_call_regions_file_path} => (
            cmd_aliases   => [qw{ mna_cr }],
            documentation => q{Path to manta call regions file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{tiddit} => (
            cmd_aliases   => [qw{ tid }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Structural variant calling using Tiddit},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{tiddit_coverage_bin_size} => (
            cmd_aliases   => [qw{ tidbin }],
            cmd_tags      => [q{Default: 500}],
            documentation => q{Size of coverage bins in calculation},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{tiddit_minimum_number_supporting_pairs} => (
            cmd_aliases   => [qw{ tidmsp }],
            cmd_flag      => q{tid_min_num_sp},
            cmd_tags      => [q{Default: 6}],
            documentation => q{Minimum number of supporting reads},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{sv_combinevariantcallsets} => (
            cmd_aliases   => [qw{ svc }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Combine structural variant call sets},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{sv_combinevariantcallsets_bcf_file} => (
            cmd_aliases => [qw{ svcbcf }],
            cmd_flag    => q{sv_comb_vcs_bf},
            documentation =>
              q{Produce a bcf from the CombineStructuralVariantCallSet vcf},
            is  => q{rw},
            isa => Bool,
        )
    );

    option(
        q{sv_svdb_merge_prioritize} => (
            cmd_aliases   => [qw{ svsvdbmp }],
            documentation => q{Prioritization order of structural variant callers},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sv_vt_decompose} => (
            cmd_aliases   => [qw{ svcvtd }],
            documentation => q{Split multi allelic records into single records},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{sv_annotate} => (
            cmd_aliases   => [qw{ svan }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Annotate and filter structural variant calls},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{sv_fqa_vcfanno_config} => (
            cmd_aliases   => [qw{ svfqav }],
            documentation => q{Frequency vcfanno toml config},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sv_fqa_annotations} => (
            cmd_aliases   => [qw{ svfqaa }],
            documentation => q{Frequency annotations to use when filtering },
            is            => q{rw},
            isa           => ArrayRef,
        )
    );

    option(
        q{sv_frequency_filter} => (
            cmd_aliases   => [qw{ svcgmf }],
            documentation => q{Remove common structural variants from vcf},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{sv_svdb_query} => (
            cmd_aliases   => [qw{ svcdbq }],
            documentation => q{Annotate structural variants using svdb query},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{sv_svdb_query_db_files} => (
            cmd_aliases   => [qw{ svcdbqd }],
            cmd_tags      => [q{file.vcf=vcf_info_key}],
            documentation => q{Database file(s) for annotation},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{vcf2cytosure_ar} => (
            cmd_aliases => [qw{ v2cs }],
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
q{Convert a VCF with structural variants to the “.CGH” format used by the commercial Cytosure software},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{vcf2cytosure_exclude_filter} => (
            cmd_aliases   => [qw{ vc2csef }],
            cmd_flag      => q{vcf2cytosure_ex_fi},
            documentation => q{Filter vcf using bcftools exclude filter string},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{vcf2cytosure_freq} => (
            cmd_aliases   => [qw{ v2csfq }],
            cmd_tags      => [q{Default: 0.01}],
            documentation => q{Specify maximum frequency},
            is            => q{rw},
            isa           => Num,
        )
    );

    option(
        q{vcf2cytosure_freq_tag} => (
            cmd_aliases   => [qw{ v2csfqt }],
            cmd_tags      => [q{Default: FRQ}],
            documentation => q{Specify frequency tag},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{vcf2cytosure_maxbnd} => (
            cmd_aliases   => [qw{ v2csmb }],
            cmd_tags      => [q{Default: 5000}],
            documentation => q{Specify maximum BND},
            is            => q{rw},
            isa           => Num,
        )
    );

    option(
        q{vcf2cytosure_no_filter} => (
            cmd_aliases   => [qw{ v2csnf }],
            documentation => q{Do not use any filtering},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{vcf2cytosure_var_size} => (
            cmd_aliases   => [qw{ v2csvs }],
            cmd_tags      => [q{Default: 5000}],
            documentation => q{Specify minimum variant size},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{sv_varianteffectpredictor} => (
            cmd_aliases   => [qw{ svv }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Annotate SV variants using VEP},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{sv_vep_features} => (
            cmd_aliases => [qw{ svvepf }],
            cmd_tags    => [
q{Default: hgvs, symbol, numbers, sift, polyphen, humdiv, domains, protein, ccds, uniprot, biotype, regulatory, tsl, canonical, per_gene, appris}
            ],
            documentation => q{VEP features},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{sv_vcfparser} => (
            cmd_aliases   => [qw{ svvcp }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Parse structural variants using vcfParser.pl},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{sv_vcfparser_add_all_mt_var} => (
            cmd_aliases   => [qw{ svvcpamt }],
            cmd_flag      => q{sv_vcfparser_all_mt},
            documentation => q{Add all MT variants in select vcf},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{sv_vcfparser_per_gene} => (
            cmd_aliases   => [qw{ svvcppg }],
            documentation => q{Keep only most severe consequence per gene},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{sv_vcfparser_range_feature_annotation_columns} => (
            cmd_aliases   => [qw{ svvcprfa }],
            cmd_flag      => q{sv_vcfparser_fac},
            documentation => q{Range annotations feature columns},
            is            => q{rw},
            isa           => ArrayRef [Int],
        )
    );

    option(
        q{sv_vcfparser_range_feature_file} => (
            cmd_aliases   => [qw{ svvcprff }],
            cmd_flag      => q{sv_vcfparser_rff},
            cmd_tags      => [q{Format: tsv}],
            documentation => q{Range annotations file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sv_vcfparser_select_feature_annotation_columns} => (
            cmd_aliases   => [qw{ svvcpsfa }],
            cmd_flag      => q{sv_vcfparser_slt_fac},
            documentation => q{Feature columns to use in annotation},
            is            => q{rw},
            isa           => ArrayRef [Int],
        )
    );

    option(
        q{sv_vcfparser_select_file} => (
            cmd_aliases   => [qw{ svvcpsf }],
            cmd_flag      => q{sv_vcfparser_slt_fl},
            cmd_tags      => [q{Format: tsv; HGNC Symbol required in file}],
            documentation => q{Select file with list of genes to analyse separately},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sv_vcfparser_select_file_matching_column} => (
            cmd_aliases   => [qw{ svvcpsfm }],
            cmd_flag      => q{sv_vcfparser_slt_fmc},
            documentation => q{Position of HGNC Symbol column in select file},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{sv_vcfparser_vep_transcripts} => (
            cmd_aliases   => [qw{ svvcvt }],
            cmd_flag      => q{sv_vcfparser_vtr},
            documentation => q{Parse VEP transcript specific entries},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{sv_rankvariant} => (
            cmd_aliases   => [qw{ svr }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Ranking of annotated SV variants},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{sv_genmod_annotate_regions} => (
            cmd_aliases => [qw{ svravanr }],
            cmd_flag    => q{sv_genmod_ann_reg},
            documentation =>
              q{Use predefined gene annotation supplied with genmod for defining genes},
            is  => q{rw},
            isa => Bool,
        )
    );

    option(
        q{sv_genmod_models_case_type} => (
            cmd_aliases   => [qw{ svravgft }],
            cmd_flag      => q{sv_genmod_mod_fam_typ},
            cmd_tags      => [q{Default: mip}],
            documentation => q{Use one of the known setups},
            is            => q{rw},
            isa           => enum( [qw{ped alt cmms mip}] ),
        )
    );

    option(
        q{sv_genmod_models_reduced_penetrance_file} => (
            cmd_aliases   => [qw{ svravrpf }],
            cmd_flag      => q{sv_genmod_mod_red_pen_f},
            documentation => q{File containing genes with reduced penetrance},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sv_genmod_models_whole_gene} => (
            cmd_aliases   => [qw{ svravwg }],
            cmd_flag      => q{sv_genmod_mod_whl_gene},
            documentation => q{Allow compound pairs in intronic regions},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{sv_rank_model_file} => (
            cmd_aliases   => [qw{ svravrm }],
            documentation => q{Rank model config file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sv_reformat} => (
            cmd_aliases   => [qw{ svre }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Concatenating files},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{sv_reformat_remove_genes_file} => (
            cmd_aliases   => [qw{ svrergf }],
            cmd_flag      => q{sv_reformat_rem_gen_f},
            documentation => q{Remove variants with hgnc_ids from file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sv_rankvariant_binary_file} => (
            cmd_aliases => [qw{ svrevbf }],
            documentation =>
              q{Produce binary file from the rank variant chromosome sorted vcfs},
            is  => q{rw},
            isa => Bool,
        )
    );

    option(
        q{bcftools_mpileup} => (
            cmd_aliases   => [qw{ bmp }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Variant calling using bcftools mpileup},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{bcftools_mpileup_constrain} => (
            cmd_aliases   => [qw{ bmpcon }],
            cmd_flag      => q{bcftools_mpileup_constrain},
            documentation => q{Use contrain in trio calling},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{bcftools_mpileup_filter_variant} => (
            cmd_aliases   => [qw{ bmpfv }],
            cmd_flag      => q{bcftools_mpileup_fil_var},
            documentation => q{Use standard bcftools filters},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{bcftools_mpileup_keep_unnormalised} => (
            cmd_aliases   => [qw{ bmpkn }],
            cmd_flag      => q{bcftools_mpileup_keep_unn},
            documentation => q{Do not normalise variants},
            is            => q{rw},
            isa           => Bool,
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
            cmd_tags      => [q{Default: GVCF}],
            documentation => q{VCF to produce},
            is            => q{rw},
            isa           => ArrayRef [ enum( [qw{ NONE BP_RESOLUTION GVCF }] ), ],
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
            cmd_tags      => [q{Default: NONE; Set to "0" to disable}],
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
            cmd_tags      => [q{Default: grch37_dbsnp_-138-.vcf}],
            documentation => q{GATK HaplotypeCaller dbSNP set for annotating ID columns},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{gatk_genotypegvcfs} => (
            cmd_aliases   => [qw{ ggt }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Merge gVCF records using GATK GenotypeGVCFs},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_genotypegvcfs_ref_gvcf} => (
            cmd_aliases => [qw{ ggtgrl }],
            cmd_flag    => q{gatk_genotype_ref_gvcf},
            documentation =>
              q{GATK GenoTypeGVCFs gVCF reference infile list for joint genotyping},
            is  => q{rw},
            isa => Str,
        )
    );

    option(
        q{gatk_genotypegvcfs_all_sites} => (
            cmd_aliases   => [qw{ ggtas }],
            cmd_flag      => q{gatk_genotype_all_sit},
            documentation => q{Include loci found to be non-variant after genotyping},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_gathervcfs} => (
            cmd_aliases   => [qw{ gcgt }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Concatenate gVCF records using GATK Concatenate variants},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_gathervcfs_bcf_file} => (
            cmd_aliases   => [qw{ gcgbcf }],
            cmd_flag      => q{gatk_genotype_bcf_f},
            documentation => q{Produce a bcf from the GATK ConcatenateGenoTypeGVCFs vcf},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_variantrecalibration} => (
            cmd_aliases => [qw{ gvr }],
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
              q{Variant recalibration using GATK VariantRecalibrator/ApplyRecalibration},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_variantrecalibration_annotations} => (
            cmd_aliases   => [qw{ gvrann }],
            cmd_flag      => q{gatk_varrecal_ann},
            cmd_tags      => [q{Default: QD, MQRankSum, ReadPosRankSum, FS, SOR, DP}],
            documentation => q{Annotations to use with GATK VariantRecalibrator},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{gatk_calculategenotypeposteriors} => (
            cmd_aliases   => [qw{ gcgp }],
            cmd_flag      => q{gatk_calculategenotypeposteriors},
            documentation => q{Perform gatk calculate genotype posterior},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_cnnscorevariants} => (
            cmd_aliases => [qw{ gcnn }],
            cmd_flag    => q{gatk_cnnscorevariants},
            documentation =>
              q{Perform gatk cnnscorevariants instead of gatk variantscore recalibration},
            is  => q{rw},
            isa => Bool,
        )
    );

    option(
        q{gatk_variantrecalibration_dp_annotation} => (
            cmd_aliases   => [qw{ gvrdpa }],
            cmd_flag      => q{gatk_varrecal_dp_ann},
            documentation => q{Use the DP annotation in variant recalibration},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_variantrecalibration_indel_max_gaussians} => (
            cmd_aliases   => [qw{ gvrimg }],
            cmd_flag      => q{gatk_varrecal_indel_max_gau},
            documentation => q{Use hard filtering for indels},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_variantrecalibration_indel_tsfilter_level} => (
            cmd_aliases => [qw{ gvritf }],
            cmd_flag    => q{gatk_varrecal_indel_ts_fl},
            cmd_tags    => [q{Defaults: 99.9}],
            documentation =>
              q{Truth sensitivity level for indels at which to start filtering},
            is  => q{rw},
            isa => Num,
        )
    );

    option(
        q{gatk_variantrecalibration_keep_unnormalised} => (
            cmd_aliases   => [qw{ gvrkn }],
            cmd_flag      => q{gatk_variantrecalibration_keep_unn},
            documentation => q{Do not normalise variants},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_variantrecalibration_resource_indel} => (
            cmd_aliases => [qw{ gvrrei }],
            cmd_flag    => q{gatk_varrecal_res_indel},
            cmd_tags    => [
q{file.vcf=settings; Default: grch37_dbsnp_-138-.vcf="dbsnp,known=true,training=false,truth=false,prior=2.0", grch37_mills_and_1000g_indels_-gold_standard-.vcf="mills,VCF,known=true,training=true,truth=true,prior=12.0"}
            ],
            documentation =>
              q{Resource to use with GATK VariantRecalibrator in INDEL|BOTH},
            is  => q{rw},
            isa => HashRef,
        )
    );

    option(
        q{gatk_variantrecalibration_resource_snv} => (
            cmd_aliases => [qw{ gvrres }],
            cmd_flag    => q{gatk_varrecal_res_snv},
            cmd_tags    => [
q{file.vcf=settings; Default: grch37_dbsnp_-138-.vcf="dbsnp,known=true,training=false,truth=false,prior=2.0", grch37_hapmap_-3.3-.vcf="hapmap,VCF,known=false,training=true,truth=true,prior=15.0", grch37_1000g_omni_-2.5-.vcf="omni,VCF,known=false,training=true,truth=false,prior=12.0", grch37_1000g_snps_high_confidence_-phase1-.vcf="1000G,known=false,training=true,truth=false,prior=10.0"}
            ],
            documentation =>
              q{Resource to use with GATK VariantRecalibrator in SNV|BOTH mode},
            is  => q{rw},
            isa => HashRef,
        )
    );

    option(
        q{gatk_variantrecalibration_snv_max_gaussians} => (
            cmd_aliases   => [qw{ gvrsmg }],
            cmd_flag      => q{gatk_varrecal_snv_max_gau},
            documentation => q{Use hard filtering for snvs},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_variantrecalibration_snv_tsfilter_level} => (
            cmd_aliases => [qw{ gvrstf }],
            cmd_flag    => q{gatk_varrecal_snv_ts_fl},
            cmd_tags    => [q{Defaults: 99.9}],
            documentation =>
              q{Truth sensitivity level for snvs at which to start filtering},
            is  => q{rw},
            isa => Num,
        )
    );

    option(
        q{gatk_variantrecalibration_ts_tranches} => (
            cmd_aliases   => [qw{ gvrtst }],
            documentation => q{Tranches to slice data},
            is            => q{rw},
            isa           => ArrayRef,
        )
    );

    option(
        q{gatk_variantrecalibration_trust_all_polymorphic} => (
            cmd_aliases   => [qw{ gvrtap }],
            cmd_flag      => q{gatk_varrecal_trust_poly},
            documentation => q{Trust all training sites to be polymorphic},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_num_reference_samples_if_no_call} => (
            cmd_aliases => [qw{ gnrsc }],
            cmd_flag    => q{gatk_num_ref_sam_if_ncall},
            cmd_tags    => [q{Defaults: 7854}],
            documentation =>
q{Number of hom-ref genotypes to infer at sites not present in a panel. Connected to option 'gatk_calculate_genotype_call_set'},
            is  => q{rw},
            isa => Int,
        )
    );

    option(
        q{gatk_calculate_genotype_call_set} => (
            cmd_aliases   => [qw{ gcgcs }],
            cmd_flag      => q{gatk_calc_gtype_cs},
            cmd_tags      => [q{Defaults: grch37_gnomad.genomes_-r2.0.1-.vcf.gz}],
            documentation => q{Callset to use in calculating genotype priors},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{gatk_combinevariantcallsets} => (
            cmd_aliases   => [qw{ gcv }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Combine variant call sets},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_combinevariantcallsets_bcf_file} => (
            cmd_aliases   => [qw{ gcvbcf }],
            cmd_flag      => q{gatk_combinevar_bcf_f},
            documentation => q{Produce a bcf from the GATK CombineVariantCallSet vcf},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_combinevariants_genotype_merge_option} => (
            cmd_aliases   => [qw{ gcvgmo }],
            cmd_flag      => q{gatk_combinevar_merge_opt},
            cmd_tags      => [q{Defaults: PRIORITIZE}],
            documentation => q{Type of merge to perform},
            is            => q{rw},
            isa           => enum( [qw{ UNIQUIFY PRIORITIZE UNSORTED REQUIRE_UNIQUE }] ),
        )
    );

    option(
        q{gatk_combinevariants_prioritize_caller} => (
            cmd_aliases   => [qw{ gcvpc }],
            cmd_flag      => q{gatk_combinevar_prio_cal},
            documentation => q{Prioritization order of variant callers},
            is            => q{rw},
            isa           => enum( [qw{ gatk bcftools }] ),
        )
    );

    option(
        q{gatk_variantevalall} => (
            cmd_aliases => [qw{ uvea }],
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
              q{Variant evaluation using GATK varianteval for all variants},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_variantevalexome} => (
            cmd_aliases => [qw{ gvee }],
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
              q{Variant evaluation using GATK varianteval for exonic variants},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_varianteval_dbsnp} => (
            cmd_aliases   => [qw{ gveedbs }],
            cmd_tags      => [q{Default: dbsnp_grch37_138_esa_129.vcf}],
            documentation => q{DbSNP file used in GATK varianteval},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{gatk_varianteval_gold} => (
            cmd_aliases => [qw{ gveedbg }],
            cmd_tags => [q{Default: grch37_mills_and_1000g_indels_-gold_standard-.vcf}],
            documentation => q{Gold indel file used in GATK varianteval},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{prepareforvariantannotationblock} => (
            cmd_aliases => [qw{ pvab }],
            cmd_flag    => q{prep_for_var_ann_bl},
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
q{Prepare for variant annotation block by copying and splitting files per contig},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{rhocall_ar} => (
            cmd_aliases => [qw{ rhc }],
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
              q{Rhocall performs annotation of variants in autozygosity regions},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{rhocall_frequency_file} => (
            cmd_aliases => [qw{ rhcf }],
            cmd_tags    => [q{Default: grch37_anon_swegen_snp_-2016-10-19-.tab.gz; tsv}],
            documentation => q{Frequency file for bcftools roh calculation},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{vt_ar} => (
            cmd_aliases   => [qw{ vt_ar }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Decompose and normalize},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{vt_decompose} => (
            cmd_aliases   => [qw{ vtddec }],
            documentation => q{Split multi allelic records into single records},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{vt_missing_alt_allele} => (
            cmd_aliases   => [qw{ vtmaa }],
            documentation => q{Remove missing alternative alleles '*'},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{vt_normalize} => (
            cmd_aliases   => [qw{ vtdnor }],
            documentation => q{Normalize variants},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{vt_uniq} => (
            cmd_aliases   => [qw{ vtunq }],
            documentation => q{Remove variant duplicates},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{rhocall_viz} => (
            cmd_aliases   => [qw{ rhv }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Create roh files needed for chromograph},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{upd_ar} => (
            cmd_aliases   => [qw{ upd }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Create bed files needed for chromograph},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{chromograph_ar} => (
            cmd_aliases   => [qw{ chgp }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Chromograph},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{frequency_filter} => (
            cmd_aliases   => [qw{ fqf }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Filter variants on frequency},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{fqf_annotations} => (
            cmd_aliases   => [qw{ fqfa }],
            documentation => q{Frequency annotations to use when filtering },
            is            => q{rw},
            isa           => ArrayRef,
        )
    );

    option(
        q{fqf_bcftools_filter_threshold} => (
            cmd_aliases   => [qw{ fqfgft }],
            cmd_flag      => q{freq_bcftools_fil_trh},
            cmd_tags      => [q{Default: 0.10}],
            documentation => q{Threshold for filtering variants},
            is            => q{rw},
            isa           => Num,
        )
    );

    option(
        q{cadd_ar} => (
            cmd_aliases   => [qw{ cad }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Annotate variants with CADD},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{cadd_column_names} => (
            cmd_aliases   => [qw{ cadc }],
            documentation => q{Column names in cadd tsv},
            is            => q{rw},
            isa           => ArrayRef,
        )
    );

    option(
        q{cadd_vcf_header_file} => (
            cmd_aliases   => [qw{ cadvh }],
            documentation => q{},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{varianteffectpredictor} => (
            cmd_aliases   => [qw{ vep }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Annotate variants using VEP},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{vep_custom_annotation} => (
            cmd_aliases   => [qw{ vepcann }],
            documentation => q{VEP custom annotation},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{vep_directory_cache} => (
            cmd_aliases   => [qw{ vepc }],
            documentation => q{Specify the cache directory to use},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{vep_features} => (
            cmd_aliases => [qw{ vepf }],
            cmd_tags    => [
q{Default: hgvs, symbol, numbers, sift, polyphen, humdiv, domains, protein, ccds, uniprot, biotype, regulatory, tsl, canonical, per_gene, appris}
            ],
            documentation => q{VEP features},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{vep_plugins_dir_path} => (
            cmd_aliases   => [qw{ veppldp }],
            documentation => q{Path to directory with VEP plugins},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{vcfparser_ar} => (
            cmd_aliases   => [qw{ vcp }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Parse structural variants using vcfParser.pl},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{vcfparser_add_all_mt_var} => (
            cmd_aliases   => [qw{ vcpamt }],
            cmd_flag      => q{vcfparser_all_mt},
            documentation => q{Add all MT variants in select vcf},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{vcfparser_range_feature_annotation_columns} => (
            cmd_aliases   => [qw{ vcprfa }],
            cmd_flag      => q{vcfparser_fac},
            documentation => q{Range annotations feature columns},
            is            => q{rw},
            isa           => ArrayRef [Int],
        )
    );

    option(
        q{vcfparser_range_feature_file} => (
            cmd_aliases   => [qw{ vcprff }],
            cmd_flag      => q{vcfparser_rff},
            cmd_tags      => [q{Format: tsv}],
            documentation => q{Range annotations file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{vcfparser_select_file} => (
            cmd_aliases   => [qw{ vcpsf }],
            cmd_flag      => q{vcfparser_slt_fl},
            cmd_tags      => [q{Format: tsv; HGNC Symbol required in file}],
            documentation => q{Select file with list of genes to analyse separately},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{vcfparser_select_feature_annotation_columns} => (
            cmd_aliases   => [qw{ vcpsfa }],
            cmd_flag      => q{vcfparser_slt_fac},
            documentation => q{Feature columns to use in annotation},
            is            => q{rw},
            isa           => ArrayRef [Int],
        )
    );

    option(
        q{vcfparser_select_file_matching_column} => (
            cmd_aliases   => [qw{ vcpsfm }],
            cmd_flag      => q{vcfparser_slt_fmc},
            documentation => q{Position of HGNC Symbol column in select file},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{vcfparser_vep_transcripts} => (
            cmd_aliases   => [qw{ vcvt }],
            cmd_flag      => q{vcfparser_vtr},
            documentation => q{Parse VEP transcript specific entries},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{rankvariant} => (
            cmd_aliases   => [qw{ rav }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Ranking of annotated variants},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{genmod_annotate_cadd_files} => (
            cmd_aliases   => [qw{ ravcad }],
            cmd_flag      => q{genmod_ann_cadd_fs},
            documentation => q{CADD score files},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{genmod_annotate_regions} => (
            cmd_aliases => [qw{ ravanr }],
            cmd_flag    => q{genmod_ann_reg},
            documentation =>
              q{Use predefined gene annotation supplied with genmod for defining genes},
            is  => q{rw},
            isa => Bool,
        )
    );

    option(
        q{genmod_annotate_spidex_file} => (
            cmd_aliases   => [qw{ ravspi }],
            cmd_flag      => q{genmod_ann_spidex_f},
            documentation => q{Spidex database for alternative splicing},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{genmod_models_case_type} => (
            cmd_aliases   => [qw{ ravgft }],
            cmd_flag      => q{genmod_mod_fam_typ},
            cmd_tags      => [q{Default: mip}],
            documentation => q{Use one of the known setups},
            is            => q{rw},
            isa           => enum( [qw{ped alt cmms mip}] ),
        )
    );

    option(
        q{genmod_models_reduced_penetrance_file} => (
            cmd_aliases   => [qw{ ravrpf }],
            cmd_flag      => q{genmod_mod_red_pen_f},
            documentation => q{File containing genes with reduced penetrance},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{genmod_models_whole_gene} => (
            cmd_aliases   => [qw{ ravwg }],
            cmd_flag      => q{genmod_mod_whl_gene},
            documentation => q{Allow compound pairs in intronic regions},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{rankvariant_binary_file} => (
            cmd_aliases => [qw{ ravbf }],
            documentation =>
              q{Produce binary file from the rank variant chromosomal sorted vcfs},
            is  => q{rw},
            isa => Bool,
        )
    );

    option(
        q{rank_model_file} => (
            cmd_aliases   => [qw{ ravrm }],
            documentation => q{Rank model config file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{endvariantannotationblock} => (
            cmd_aliases   => [qw{ evab }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{End variant annotation block by concatenating files},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{endvariantannotationblock_remove_genes_file} => (
            cmd_aliases   => [qw{ evabrgf }],
            cmd_flag      => q{endvarannbl_rem_gen_f},
            documentation => q{Remove variants with hgnc_ids from file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{peddy_ar} => (
            cmd_aliases   => [qw{ pedd }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{QC for familial-relationships and sexes},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{plink} => (
            cmd_aliases   => [qw{ plink }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{QC for samples gender and relationship},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{variant_integrity_ar} => (
            cmd_aliases   => [qw{ vai }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{QC for samples relationship},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{rtg_vcfeval} => (
            cmd_aliases   => [qw{ rte }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Compare concordance with benchmark data set},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{nist_call_set_vcf} => (
            cmd_aliases   => [qw{ nist_csv }],
            cmd_tags      => [q{Nist call set vcf information hash}],
            documentation => q{NIST high-confidence variant calls vcf},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{nist_call_set_bed} => (
            cmd_aliases   => [qw{ nist_csb }],
            cmd_tags      => [q{Nist call set bed information hash}],
            documentation => q{NIST high-confidence variant calls bed},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{nist_id} => (
            cmd_aliases   => [qw{ nist_id }],
            cmd_tags      => [q{sample_id=nist_id}],
            documentation => q{Map sample_id to nist_id},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{nist_versions} => (
            cmd_aliases   => [qw{ nist_versions }],
            cmd_tags      => [q{Default: [2.19, 3.3.2]}],
            documentation => q{Map sample_id to nist_id},
            is            => q{rw},
            isa           => ArrayRef,
        )
    );

    option(
        q{qccollect_ar} => (
            cmd_aliases   => [qw{ qcc }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Collect QC metrics from recipes output},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{qccollect_regexp_file} => (
            cmd_aliases => [qw{ qccref }],
            cmd_tags    => [q{Default: qc_regexp_-v1.18-.yaml}],
            documentation =>
q{Regular expression file containing the regular expression to be used for each program},
            is  => q{rw},
            isa => Str,
        )
    );

    option(
        q{qccollect_sampleinfo_file} => (
            cmd_aliases => [qw{ qccsi }],
            cmd_tags =>
              [q{Default: {outdata_dir}/{case_id}/{case_id}_qc_sample_info.yaml}],
            documentation =>
              q{Sample info file containing info on what to parse from this analysis run},
            is  => q{rw},
            isa => Str,
        )
    );

    option(
        q{qccollect_skip_evaluation} => (
            cmd_aliases   => [qw{ qccske }],
            documentation => q{Skip evaluation step in qccollect},
            is            => q{rw},
            isa           => Bool,
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
        q{samtools_subsample_mt} => (
            cmd_aliases   => [qw{ ssmt }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Subsample the mitochondria reads},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{samtools_subsample_mt_depth} => (
            cmd_aliases   => [qw{ ssmtd }],
            cmd_tags      => [q{Default: 60}],
            documentation => q{Set approximate coverage of subsampled bam file},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{varg_ar} => (
            cmd_aliases   => [qw{ varg }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Compare resulting SVs and SNVs with positive controls},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{varg_truth_set_vcf} => (
            cmd_aliases   => [qw{ vts }],
            cmd_tags      => [q{Format: vcf}],
            documentation => q{vcf with expected SVs and SNVs},
            is            => q{rw},
            isa           => Str,
        )
    );

    return;
}

1;
