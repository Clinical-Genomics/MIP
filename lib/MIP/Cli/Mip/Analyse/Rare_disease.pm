package MIP::Cli::Mip::Analyse::Rare_disease;

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

our $VERSION = 0.01;

extends(qw{ MIP::Cli::Mip::Analyse });

command_short_description(q{Rare disease analysis});

command_long_description(
    q{Rare disease analysis on wes, wgs or mixed sequence data});

command_usage(
    q{mip <analyse> <rare_disease> <family_id> --config <config_file> });

## Define, check and get Cli supplied parameters
_build_usage();

sub run {
    my ($arg_href) = @_;

    ## Remove Moose::App extra variable
    delete $arg_href->{extra_argv};

    ## Input from Cli
    my %active_parameter = %{$arg_href};

    use MIP::File::Format::Parameter qw{ parse_definition_file  };
    use MIP::File::Format::Yaml qw{ order_parameter_names };
    use MIP::Get::Analysis qw{ print_program };

    ## Mip analyse rare_disease parameters
    ## The order of files in @definition_files should follow commands inheritance
    my @definition_files = (
        catfile( $Bin, qw{ definitions mip_parameters.yaml } ),
        catfile( $Bin, qw{ definitions rare_disease_parameters.yaml } ),
    );

    ## Non mandatory parameter definition keys to check
    my $non_mandatory_parameter_keys_path =
      catfile( $Bin, qw{ definitions non_mandatory_parameter_keys.yaml } );

    ## Mandatory parameter definition keys to check
    my $mandatory_parameter_keys_path =
      catfile( $Bin, qw{ definitions mandatory_parameter_keys.yaml } );

    ## %parameter holds all defined parameters for MIP
    ## mip analyse rare_disease parameters
    my %parameter;
    foreach my $definition_file (@definition_files) {

        %parameter = (
            %parameter,
            parse_definition_file(
                {
                    define_parameters_path => $definition_file,
                    non_mandatory_parameter_keys_path =>
                      $non_mandatory_parameter_keys_path,
                    mandatory_parameter_keys_path =>
                      $mandatory_parameter_keys_path,
                }
            )
        );
    }

    ## Print programs and exit
    if ( $active_parameter{print_programs} ) {

        print_program(
            {
                define_parameters_files_ref => \@definition_files,
                parameter_href              => \%parameter,
                print_program_mode => $active_parameter{print_program_mode},
            }
        );
        exit;
    }

    ### To add/write parameters in the correct order
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

        # BWA human genome reference file endings
        bwa_build_reference => [qw{ .bwt .ann .amb .pac .sa }],

        exome_target_bed =>
          [qw{ .infile_list .pad100.infile_list .pad100.interval_list }],

        # Human genome meta files
        human_genome_reference_file_endings => [qw{ .dict .fai }],

        # RTG human genome reference file endings
        rtg_vcfeval_reference_genome => [qw{ _sdf_dir }],
    );

    mip_analyse(
        {
            active_parameter_href => \%active_parameter,
            file_info_href        => \%file_info,
            parameter_href        => \%parameter,
            order_parameters_ref  => \@order_parameters,
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
q{gatk_realigner_indel_known_sites, gatk_baserecalibration_known_sites, gatk_haplotypecaller_snp_known_set, gatk_variantrecalibration_resource_snv, gatk_variantrecalibration_resource_indel, frequency_genmod_filter_1000g, sv_vcfanno_config_file, gatk_varianteval_gold, gatk_varianteval_dbsnp, snpsift_annotation_files}
            ],
            documentation =>
              q{Set the references to be decomposed and normalized},
            is  => q{rw},
            isa => ArrayRef [Str],
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
        q{exome_target_bed} => (
            cmd_aliases => [qw{ extb }],
            cmd_tags =>
              [q{file.bed=Sample_id; Default: latest_supported_capturekit.bed}],
            documentation => q{Exome target bed file per sample id},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{supported_capture_kit} => (
            cmd_aliases => [qw{ sck }],
            cmd_tags    => [q{acronym=file.bed}],
            documentation =>
              q{Set the capture kit acronym shortcut in pedigree file},
            is  => q{rw},
            isa => HashRef,
        )
    );

    option(
        q{pbwa_mem} => (
            cmd_aliases   => [qw{ pmem }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Align reads using Bwa Mem},
            is            => q{rw},
            isa           => enum( [ 1, 2 ] ),
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
        q{bwa_mem_cram} => (
            cmd_aliases   => [qw{ memcrm }],
            documentation => q{Use CRAM-format for additional output file},
            is            => q{rw},
            isa           => Bool,
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
        q{ppicardtools_mergesamfiles} => (
            cmd_aliases => [qw{ pptm }],
            cmd_flag    => q{ppicard_mergesamfiles},
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
q{Merge (BAM file(s) ) or rename single samples for downstream processing},
            is  => q{rw},
            isa => enum( [ 1, 2 ] ),
        )
    );

    option(
        q{pmarkduplicates} => (
            cmd_aliases   => [qw{ pmd }],
            cmd_flag      => q{pmarkduplicates},
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
            cmd_aliases => [qw{ mdshts }],
            cmd_flag    => q{sba_mdup_hts},
            cmd_tags    => [q{Default: 262144}],
            documentation =>
              q{Sambamba size of hash table for finding read pairs},
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
        q{pchanjo_sexcheck} => (
            cmd_aliases   => [qw{ pchs }],
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
        q{psambamba_depth} => (
            cmd_aliases   => [qw{ psdt }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Sambamba depth coverage analysis},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
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
        q{sambamba_depth_cutoffs} => (
            cmd_aliases   => [qw{ sdtcut }],
            cmd_flag      => q{sba_depth_co},
            documentation => q{Read depth cutoff},
            is            => q{rw},
            isa           => ArrayRef [Int],
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
        q{sambamba_depth_noduplicates} => (
            cmd_aliases => [qw{ sdtndu }],
            cmd_flag    => q{sba_depth_nod},
            documentation =>
              q{Do not include duplicates in coverage calculation},
            is  => q{rw},
            isa => Bool,
        )
    );

    option(
        q{sambamba_depth_quality_control} => (
            cmd_aliases => [qw{ sdtfqc }],
            cmd_flag    => q{sba_depth_qc},
            documentation =>
              q{Do not include reads with failed quality control},
            is  => q{rw},
            isa => Bool,
        )
    );

    option(
        q{pbedtools_genomecov} => (
            cmd_aliases => [qw{ pbgc }],
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
              q{Genome coverage calculation using bedtools genomecov},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{bedtools_genomecov_max_coverage} => (
            cmd_aliases   => [qw{ bgcmc }],
            cmd_flag      => q{bgc_max_cov},
            cmd_tags      => [q{Default: 30}],
            documentation => q{Max coverage depth},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{ppicardtools_collectmultiplemetrics} => (
            cmd_aliases   => [qw{ pptcmm }],
            cmd_flag      => q{ppt_col_mul_met},
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Qc metrics calculation},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{ppicardtools_collecthsmetrics} => (
            cmd_aliases   => [qw{ pptchs }],
            cmd_flag      => q{ppt_col_hs_met},
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Qc metrics calculation for capture},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{prcovplots} => (
            cmd_aliases   => [qw{ prcp }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Plots of genome coverage using rcovplots},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{pcnvnator} => (
            cmd_aliases   => [qw{ pcnv }],
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
        q{pdelly_call} => (
            cmd_aliases   => [qw{ pdelc }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Structural variant calling using Delly},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{pdelly_reformat} => (
            cmd_aliases   => [qw{ pdel }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Merge, regenotype and filter using Delly},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
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
        q{delly_exclude_file} => (
            cmd_aliases => [qw{ delexc }],
            cmd_tags    => [q{Default: hg19_human_excl_-0.7.6-.tsv}],
            documentation =>
              q{Exclude centomere and telemore regions in delly calling},
            is  => q{rw},
            isa => Str,
        )
    );

    option(
        q{pmanta} => (
            cmd_aliases   => [qw{ pmna }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Structural variant calling using Manta},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{ptiddit} => (
            cmd_aliases   => [qw{ ptid }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Structural variant calling using Tiddit},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
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
        q{tiddit_bin_size} => (
            cmd_aliases   => [qw{ tidbin }],
            cmd_tags      => [q{Default: 500}],
            documentation => q{Size of coverage bins in calculation},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{psv_combinevariantcallsets} => (
            cmd_aliases   => [qw{ psvc }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Combine structural variant call sets},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
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
        q{sv_svdb_merge_prioritize} => (
            cmd_aliases => [qw{ svsvdbmp }],
            documentation =>
              q{Prioritization order of structural variant callers},
            is  => q{rw},
            isa => Str,
        )
    );

    option(
        q{sv_bcftools_view_filter} => (
            cmd_aliases => [qw{ svcbtv }],
            documentation =>
              q{Include structural variants with PASS in FILTER column},
            is  => q{rw},
            isa => Bool,
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
        q{sv_vcfanno} => (
            cmd_aliases   => [qw{ svcvan }],
            documentation => q{Annotate structural variants},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{sv_vcfanno_lua} => (
            cmd_aliases   => [qw{ svcval }],
            documentation => q{VcfAnno lua postscripting file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sv_vcfanno_config} => (
            cmd_aliases   => [qw{ svcvac }],
            documentation => q{VcfAnno toml config},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sv_vcfanno_config_file} => (
            cmd_aliases => [qw{ svcvacf }],
            cmd_tags =>
              [q{Default: GRCh37_all_sv_-phase3_v2.2013-05-02-.vcf.gz}],
            documentation => q{Annotation file within vcfAnno config toml file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sv_vcfannotation_header_lines_file} => (
            cmd_aliases => [qw{ svcvah }],
            cmd_flag    => q{sv_vcfanno_hlf},
            documentation =>
              q{Adjust for postscript by adding required header lines to vcf},
            is  => q{rw},
            isa => Str,
        )
    );

    option(
        q{sv_genmod_filter} => (
            cmd_aliases   => [qw{ svcgmf }],
            documentation => q{Remove common structural variants from vcf},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{sv_genmod_filter_1000g} => (
            cmd_aliases => [qw{ svcgfr }],
            cmd_tags =>
              [q{Default: GRCh37_all_wgs_-phase3_v5b.2013-05-02-.vcf.gz}],
            documentation =>
              q{Genmod annotate structural variants from 1000G reference},
            is  => q{rw},
            isa => Str,
        )
    );

    option(
        q{sv_genmod_filter_threshold} => (
            cmd_aliases   => [qw{ svcgft }],
            cmd_tags      => [q{Default: 0.10}],
            documentation => q{Threshold for filtering structural variants},
            is            => q{rw},
            isa           => Num,
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
        q{pvcf2cytosure} => (
            cmd_aliases => [qw{ pv2cs }],
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
        q{vf2cytosure_no_filter} => (
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
        q{psv_varianteffectpredictor} => (
            cmd_aliases   => [qw{ psvv }],
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
        q{sv_vep_plugins} => (
            cmd_aliases   => [qw{ svvepl }],
            cmd_tags      => [q{Default: UpDownDistance, LoFtool}],
            documentation => q{VEP plugins},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{psv_vcfparser} => (
            cmd_aliases   => [qw{ psvvcp }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Parse structural variants using vcfParser.pl},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{sv_vcfparser_vep_transcripts} => (
            cmd_aliases   => [qw{ svvcpvt }],
            cmd_flag      => q{sv_vcfparser_vtr},
            documentation => q{Parse VEP transcript specific entries},
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
        q{sv_vcfparser_range_feature_annotation_columns} => (
            cmd_aliases   => [qw{ svvcprfa }],
            cmd_flag      => q{sv_vcfparser_fac},
            documentation => q{Range annotations feature columns},
            is            => q{rw},
            isa           => ArrayRef [Int],
        )
    );

    option(
        q{sv_vcfparser_select_file} => (
            cmd_aliases => [qw{ svvcpsf }],
            cmd_flag    => q{sv_vcfparser_slt_fl},
            cmd_tags    => [q{Format: tsv; HGNC Symbol required in file}],
            documentation =>
              q{Select file with list of genes to analyse separately},
            is  => q{rw},
            isa => Str,
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
        q{sv_vcfparser_select_feature_annotation_columns} => (
            cmd_aliases   => [qw{ svvcpsfa }],
            cmd_flag      => q{sv_vcfparser_slt_fac},
            documentation => q{Feature columns to use in annotation},
            is            => q{rw},
            isa           => ArrayRef [Int],
        )
    );

    option(
        q{psv_rankvariant} => (
            cmd_aliases   => [qw{ psvr }],
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
        q{sv_genmod_models_family_type} => (
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
        q{psv_reformat} => (
            cmd_aliases   => [qw{ psvre }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Concatenating files},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
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
        q{sv_reformat_remove_genes_file} => (
            cmd_aliases   => [qw{ svrergf }],
            cmd_flag      => q{sv_reformat_rem_gen_f},
            documentation => q{Remove variants with hgnc_ids from file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{pbcftools_mpileup} => (
            cmd_aliases   => [qw{ pbmp }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Variant calling using bcftools mpileup},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{bcftools_mpileup_filter_variant} => (
            cmd_aliases   => [qw{ pbmpfv }],
            cmd_flag      => q{bcftools_mpileup_fil_var},
            documentation => q{Use standard bcftools filters},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{pfreebayes} => (
            cmd_aliases   => [qw{ pfrb }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Variant calling using Freebayes},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{pgatk_realigner} => (
            cmd_aliases => [qw{ pgra }],
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
q{Realignments of reads using GATK ReAlignerTargetCreator/IndelRealigner},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_realigner_indel_known_sites} => (
            cmd_aliases => [qw{ graks }],
            cmd_flag    => q{gatk_realigner_ind_ks},
            cmd_tags    => [
q{Default: GRCh37_1000g_indels_-phase1-.vcf, GRCh37_mills_and_1000g_indels_-gold_standard-.vcf}
            ],
            documentation =>
              q{GATK ReAlignerTargetCreator/IndelRealigner known indel site},
            is  => q{rw},
            isa => ArrayRef [Str],
        )
    );

    option(
        q{pgatk_baserecalibration} => (
            cmd_aliases => [qw{ pgbr }],
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
        q{gatk_baserecalibration_known_sites} => (
            cmd_aliases => [qw{ gbrkst }],
            cmd_flag    => q{gatk_baserecal_ks},
            cmd_tags    => [
q{Default: GRCh37_dbsnp_-138-.vcf, GRCh37_1000g_indels_-phase1-.vcf, GRCh37_mills_and_1000g_indels_-gold_standard-.vcf}
            ],
            documentation =>
              q{GATK BaseReCalibration known SNV and INDEL sites},
            is  => q{rw},
            isa => ArrayRef [Str],
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
        q{gatk_baserecalibration_disable_indel_qual} => (
            cmd_aliases   => [qw{ gbrdiq }],
            cmd_flag      => q{gatk_baserecal_dis_indel_q},
            documentation => q{Disable indel quality scores},
            is            => q{rw},
            isa           => Bool,
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
        q{pgatk_haplotypecaller} => (
            cmd_aliases   => [qw{ pghc }],
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
        q{gatk_haplotypecaller_snp_known_set} => (
            cmd_aliases => [qw{ ghckse }],
            cmd_flag    => q{gatk_haplotype_snp_ks},
            cmd_tags    => [q{Default: GRCh37_dbsnp_-138-.vcf}],
            documentation =>
              q{GATK HaplotypeCaller dbSNP set for annotating ID columns},
            is  => q{rw},
            isa => Str,
        )
    );

    option(
        q{gatk_haplotypecaller_no_soft_clipped_bases} => (
            cmd_aliases => [qw{ ghcscb }],
            cmd_flag    => q{gatk_haplotype_no_soft_cb},
            documentation =>
              q{Do not include soft clipped bases in the variant calling},
            is  => q{rw},
            isa => Bool,
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
        q{pgatk_genotypegvcfs} => (
            cmd_aliases   => [qw{ pggt }],
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
            cmd_aliases   => [qw{ ggtals }],
            documentation => q{Emit non-variant sites to the output vcf file},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_concatenate_genotypegvcfs_bcf_file} => (
            cmd_aliases => [qw{ ggbcf }],
            cmd_flag    => q{gatk_genotype_bcf_f},
            documentation =>
              q{Produce a bcf from the GATK ConcatenateGenoTypeGVCFs vcf},
            is  => q{rw},
            isa => Bool,
        )
    );

    option(
        q{pgatk_variantrecalibration} => (
            cmd_aliases => [qw{ pgvr }],
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
q{Variant recalibration using GATK VariantRecalibrator/ApplyRecalibration},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_variantrecalibration_annotations} => (
            cmd_aliases => [qw{ gvrann }],
            cmd_flag    => q{gatk_varrecal_ann},
            cmd_tags =>
              [q{Default: QD, MQRankSum, ReadPosRankSum, FS, SOR, DP}],
            documentation =>
              q{Annotations to use with GATK VariantRecalibrator},
            is  => q{rw},
            isa => ArrayRef [Str],
        )
    );

    option(
        q{gatk_variantrecalibration_resource_snv} => (
            cmd_aliases => [qw{ gvrres }],
            cmd_flag    => q{gatk_varrecal_res_snv},
            cmd_tags    => [
q{file.vcf=settings; Default: GRCh37_dbsnp_-138-.vcf="dbsnp,known=true,training=false,truth=false,prior=2.0", GRCh37_hapmap_-3.3-.vcf="hapmap,VCF,known=false,training=true,truth=true,prior=15.0", GRCh37_1000g_omni_-2.5-.vcf="omni,VCF,known=false,training=true,truth=false,prior=12.0", GRCh37_1000g_snps_high_confidence_-phase1-.vcf="1000G,known=false,training=true,truth=false,prior=10.0"}
            ],
            documentation =>
              q{Resource to use with GATK VariantRecalibrator in SNV|BOTH mode},
            is  => q{rw},
            isa => HashRef,
        )
    );

    option(
        q{gatk_variantrecalibration_resource_indel} => (
            cmd_aliases => [qw{ gvrrei }],
            cmd_flag    => q{gatk_varrecal_res_indel},
            cmd_tags    => [
q{file.vcf=settings; Default: GRCh37_dbsnp_-138-.vcf="dbsnp,known=true,training=false,truth=false,prior=2.0", GRCh37_mills_and_1000g_indels_-gold_standard-.vcf="mills,VCF,known=true,training=true,truth=true,prior=12.0"}
            ],
            documentation =>
              q{Resource to use with GATK VariantRecalibrator in INDEL|BOTH},
            is  => q{rw},
            isa => HashRef,
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
        q{gatk_variantrecalibration_dp_annotation} => (
            cmd_aliases   => [qw{ gvrdpa }],
            cmd_flag      => q{gatk_varrecal_dp_ann},
            documentation => q{Use the DP annotation in variant recalibration},
            is            => q{rw},
            isa           => Bool,
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
        q{gatk_variantrecalibration_indel_max_gaussians} => (
            cmd_aliases   => [qw{ gvrimg }],
            cmd_flag      => q{gatk_varrecal_indel_max_gau},
            documentation => q{Use hard filtering for indels},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{gatk_calculategenotypeposteriors_support_set} => (
            cmd_aliases => [qw{ gcgpss }],
            cmd_flag    => q{gatk_calcgenotypepost_ss},
            cmd_tags =>
              [q{Defaults: 1000g_sites_GRCh37_phase3_v4_20130502.vcf}],
            documentation => q{GATK CalculateGenotypePosteriors support set},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{pgatk_combinevariantcallsets} => (
            cmd_aliases   => [qw{ pgcv }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Combine variant call sets},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{gatk_combinevariants_genotype_merge_option} => (
            cmd_aliases   => [qw{ gcvgmo }],
            cmd_flag      => q{gatk_combinevar_merge_opt},
            cmd_tags      => [q{Defaults: PRIORITIZE}],
            documentation => q{Type of merge to perform},
            is            => q{rw},
            isa => enum( [qw{ UNIQUIFY PRIORITIZE UNSORTED REQUIRE_UNIQUE }] ),
        )
    );

    option(
        q{gatk_combinevariants_prioritize_caller} => (
            cmd_aliases   => [qw{ gcvpc }],
            cmd_flag      => q{gatk_combinevar_prio_cal},
            documentation => q{Prioritization order of variant callers},
            is            => q{rw},
            isa           => enum( [qw{ gatk bcftools freebayes }] ),
        )
    );

    option(
        q{gatk_combinevariantcallsets_bcf_file} => (
            cmd_aliases => [qw{ gcvbcf }],
            cmd_flag    => q{gatk_combinevar_bcf_f},
            documentation =>
              q{Produce a bcf from the GATK CombineVariantCallSet vcf},
            is  => q{rw},
            isa => Bool,
        )
    );

    option(
        q{pgatk_variantevalall} => (
            cmd_aliases => [qw{ pgvea }],
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
              q{Variant evaluation using GATK varianteval for all variants},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{pgatk_variantevalexome} => (
            cmd_aliases => [qw{ pgvee }],
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
            cmd_tags      => [q{Default: dbsnp_GRCh37_138_esa_129.vcf}],
            documentation => q{DbSNP file used in GATK varianteval},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{gatk_varianteval_gold} => (
            cmd_aliases => [qw{ gveedbg }],
            cmd_tags =>
              [q{Default: GRCh37_mills_and_1000g_indels_-gold_standard-.vcf}],
            documentation => q{Gold indel file used in GATK varianteval},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{pprepareforvariantannotationblock} => (
            cmd_aliases => [qw{ ppvab }],
            cmd_flag    => q{prep_for_var_ann_bl},
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
q{Prepare for variant annotation block by copying and splitting files per contig},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{prhocall} => (
            cmd_aliases => [qw{ prhc }],
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
            cmd_tags =>
              [q{Default: GRCh37_anon_swegen_snp_-2016-10-19-.tab.gz; tsv}],
            documentation => q{Frequency file for bcftools roh calculation},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{pvt} => (
            cmd_aliases   => [qw{ pvt }],
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
        q{vt_missing_alt_allele} => (
            cmd_aliases   => [qw{ vtmaa }],
            documentation => q{Remove missing alternative alleles '*'},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{pfrequency_filter} => (
            cmd_aliases   => [qw{ pfqf }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Filter variants on frequency},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{frequency_genmod_filter} => (
            cmd_aliases   => [qw{ fqfgmf }],
            documentation => q{Remove common variants from vcf file},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{frequency_genmod_filter_1000g} => (
            cmd_aliases => [qw{ fqfgfr }],
            cmd_flag    => q{freq_genmod_fil_1000g},
            cmd_tags =>
              [q{Default: GRCh37_all_wgs_-phase3_v5b.2013-05-02-.vcf.gz}],
            documentation => q{Genmod annotate 1000G reference},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{frequency_genmod_filter_max_af} => (
            cmd_aliases   => [qw{ fqfmaf }],
            cmd_flag      => q{freq_genmod_fil_maf},
            documentation => q{Annotate MAX_AF from reference},
            is            => q{rw},
            isa           => Num,
        )
    );

    option(
        q{frequency_genmod_filter_threshold} => (
            cmd_aliases   => [qw{ fqfgft }],
            cmd_flag      => q{freq_genmod_fil_trh},
            cmd_tags      => [q{Default: 0.10}],
            documentation => q{Threshold for filtering variants},
            is            => q{rw},
            isa           => Num,
        )
    );

    option(
        q{pvarianteffectpredictor} => (
            cmd_aliases   => [qw{ pvep }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Annotate variants using VEP},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{vep_directory_path} => (
            cmd_aliases   => [qw{ vepp }],
            documentation => q{Path to VEP script directory},
            is            => q{rw},
            isa           => Str,
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
        q{vep_plugins_dir_path} => (
            cmd_aliases   => [qw{ veppldp }],
            documentation => q{Path to directory with VEP plugins},
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
        q{vep_plugins} => (
            cmd_aliases   => [qw{ veppl }],
            cmd_tags      => [q{Default: UpDownDistance, LoFtool}],
            documentation => q{VEP plugins},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{pvcfparser} => (
            cmd_aliases   => [qw{ pvcp }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Parse structural variants using vcfParser.pl},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{vcfparser_vep_transcripts} => (
            cmd_aliases   => [qw{ vcpvt }],
            cmd_flag      => q{vcfparser_vtr},
            documentation => q{Parse VEP transcript specific entries},
            is            => q{rw},
            isa           => Bool,
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
        q{vcfparser_range_feature_annotation_columns} => (
            cmd_aliases   => [qw{ vcprfa }],
            cmd_flag      => q{vcfparser_fac},
            documentation => q{Range annotations feature columns},
            is            => q{rw},
            isa           => ArrayRef [Int],
        )
    );

    option(
        q{vcfparser_select_file} => (
            cmd_aliases => [qw{ vcpsf }],
            cmd_flag    => q{vcfparser_slt_fl},
            cmd_tags    => [q{Format: tsv; HGNC Symbol required in file}],
            documentation =>
              q{Select file with list of genes to analyse separately},
            is  => q{rw},
            isa => Str,
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
        q{vcfparser_select_feature_annotation_columns} => (
            cmd_aliases   => [qw{ vcpsfa }],
            cmd_flag      => q{vcfparser_slt_fac},
            documentation => q{Feature columns to use in annotation},
            is            => q{rw},
            isa           => ArrayRef [Int],
        )
    );

    option(
        q{psnpeff} => (
            cmd_aliases   => [qw{ psne }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Variant annotation using snpEff},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{snpeff_path} => (
            cmd_aliases   => [qw{ snep }],
            documentation => q{Path to snpEff},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{snpeff_ann} => (
            cmd_aliases   => [qw{ sneann }],
            documentation => q{Annotate variants using snpeff},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{snpeff_genome_build_version} => (
            cmd_aliases   => [qw{ snegbv }],
            cmd_tags      => [q{Default: GRCh37.75}],
            documentation => q{Snpeff genome build version},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{snpsift_annotation_files} => (
            cmd_aliases => [qw{ snesaf }],
            cmd_tags    => [
q{Default: GRCh37_all_wgs_-phase3_v5b.2013-05-02-.vcf.gz=AF, GRCh37_exac_reheader_-r0.3.1-.vcf.gz=AF, GRCh37_anon-swegen_snp_-1000samples-.vcf.gz=AF, GRCh37_anon-swegen_indel_-1000samples-.vcf.gz=AF}
            ],
            documentation => q{Annotation files to use with snpsift},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{snpsift_annotation_outinfo_key} => (
            cmd_aliases => [qw{ snesaoi }],
            cmd_flag    => q{snpsift_ann_oik},
            cmd_tags    => [
q{Default: GRCh37_all_wgs_-phase3_v5b.2013-05-02-.vcf=1000G, GRCh37_exac_reheader_-r0.3.1-.vcf.gz=EXAC, GRCh37_anon-swegen_snp_-1000samples-.vcf.gz=SWEREF, GRCh37_anon-swegen_indel_-1000samples-.vcf.gz=SWEREF}
            ],
            documentation => q{Snpsift output INFO key},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{snpsift_dbnsfp_file} => (
            cmd_aliases   => [qw{ snesdbnsfp }],
            cmd_tags      => [q{Default: GRCh37_dbnsfp_-v2.9-.txt.gz}],
            documentation => q{DbNSFP File},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{snpsift_dbnsfp_annotations} => (
            cmd_aliases => [qw{ snesdbnsfpa }],
            cmd_flag    => q{snpsift_dbnsfp_ann},
            cmd_tags    => [
q{Default: SIFT_pred, Polyphen2_HDIV_pred, Polyphen2_HVAR_pred, GERP++_NR, GERP++_RS, phastCons100way_vertebrate}
            ],
            documentation => q{DbNSFP annotations to use with snpsift},
            is            => q{rw},
            isa           => ArrayRef,
        )
    );

    option(
        q{prankvariant} => (
            cmd_aliases   => [qw{ prav }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Ranking of annotated variants},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
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
        q{genmod_annotate_cadd_files} => (
            cmd_aliases   => [qw{ ravcad }],
            cmd_flag      => q{genmod_mod_ann_cadd_fs},
            documentation => q{CADD score files},
            is            => q{rw},
            isa           => ArrayRef [Str],
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
        q{genmod_annotate_spidex_file} => (
            cmd_aliases   => [qw{ ravspi }],
            cmd_flag      => q{genmod_ann_spidex_f},
            documentation => q{Spidex database for alternative splicing},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{genmod_models_family_type} => (
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
        q{rank_model_file} => (
            cmd_aliases   => [qw{ ravrm }],
            documentation => q{Rank model config file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{pendvariantannotationblock} => (
            cmd_aliases => [qw{ pevab }],
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
              q{End variant annotation block by concatenating files},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
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
        q{rankvariant_binary_file} => (
            cmd_aliases => [qw{ ravbf }],
            documentation =>
q{Produce binary file from the rank variant chromosomal sorted vcfs},
            is  => q{rw},
            isa => Bool,
        )
    );

    option(
        q{ppeddy} => (
            cmd_aliases   => [qw{ pped }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{QC for familial-relationships and sexes},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{pplink} => (
            cmd_aliases   => [qw{ pplink }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{QC for samples gender and relationship},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{pvariant_integrity} => (
            cmd_aliases   => [qw{ pvai }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{QC for samples relationship},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{prtg_vcfeval} => (
            cmd_aliases   => [qw{ prte }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Compare concordance with benchmark data set},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{pevaluation} => (
            cmd_aliases   => [qw{ pevl }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Compare concordance with NIST data set},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{nist_id} => (
            cmd_aliases   => [qw{ evlnid }],
            cmd_tags      => [q{Default: NA12878}],
            documentation => q{NIST high-confidence sample_id},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{nist_high_confidence_call_set} => (
            cmd_aliases => [qw{ evlnhc }],
            cmd_flag    => q{nist_hc_call_s},
            cmd_tags    => [q{Default: GRCh37_nist_hg001_-na12878_v2.19-.vcf}],
            documentation => q{NIST high-confidence variant calls},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{nist_high_confidence_call_set_bed} => (
            cmd_aliases => [qw{ evlnil }],
            cmd_flag    => q{nist_hc_call_sb},
            cmd_tags    => [q{Default: GRCh37_nist_hg001_-na12878_v2.19-.bed}],
            documentation =>
              q{NIST high-confidence variant calls interval list},
            is  => q{rw},
            isa => Str,
        )
    );

    option(
        q{pqccollect} => (
            cmd_aliases   => [qw{ pqcc }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Collect QC metrics from programs output},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{qccollect_sampleinfo_file} => (
            cmd_aliases => [qw{ qccsi }],
            cmd_tags    => [
q{Default: {outdata_dir}/{family_id}/{family_id}_qc_sample_info.yaml}
            ],
            documentation =>
q{Sample info file containing info on what to parse from this analysis run},
            is  => q{rw},
            isa => Str,
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
        q{qccollect_skip_evaluation} => (
            cmd_aliases   => [qw{ qccske }],
            documentation => q{Skip evaluation step in qccollect},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{pmultiqc} => (
            cmd_aliases => [qw{ pmqc }],
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
q{Create aggregate bioinformatics analysis report across many samples},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{panalysisrunstatus} => (
            cmd_aliases => [qw{ pars }],
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
q{Check analysis output and sets the analysis run status flag to finished in sample_info_file},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{psacct} => (
            cmd_aliases => [qw{ psac }],
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
        q{psamtools_subsample_mt} => (
            cmd_aliases   => [qw{ pssmt }],
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

    return;
}

1;
