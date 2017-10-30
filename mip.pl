#!/usr/bin/env perl

#### Master script for analysing paired end reads from the Illumina plattform in fastq(.gz) format to annotated ranked disease causing variants. The program performs QC, aligns reads using BWA, performs variant discovery and annotation as well as ranking the found variants according to disease potential.

#### Copyright 2011 Henrik Stranneheim

# Require at least perl 5.18
use 5.018;
use Modern::Perl qw{ 2014 };
use autodie qw{ open close :all };

# Required for autodie :all
use IPC::System::Simple;
use English qw{ -no_match_vars };
use Carp;

## Unicode boilerplate
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };

use Getopt::Long;
use POSIX;
use Params::Check qw{ check allow last_error };
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case
use Cwd;
use Cwd qw{ abs_path };
use File::Basename qw{ dirname basename fileparse };
use File::Spec::Functions qw{ catdir catfile devnull };
use File::Path qw{ make_path };
use File::Copy qw{ copy };
use FindBin qw{ $Bin };
use IPC::Cmd qw{ can_run run};
use Time::Piece;

## Third party module(s)
use Path::Iterator::Rule;
use List::MoreUtils qw { any uniq all };
use Readonly;

##MIPs lib/
# Add MIPs internal lib
use lib catdir( $Bin, q{lib} );
use MIP::Check::Cluster qw{ check_max_core_number };
use MIP::Check::Modules qw{ check_perl_modules };
use MIP::Check::Parameter
  qw{ check_parameter_hash check_allowed_temp_directory };
use MIP::Check::Reference
  qw{ check_references_for_vt check_human_genome_prerequisites check_capture_file_prerequisites };
use MIP::Delete::List qw{ delete_non_wes_contig delete_male_contig };
use MIP::File::Format::Pedigree qw{ create_fam_file };
use MIP::File::Format::Yaml qw{ load_yaml write_yaml order_parameter_names };
use MIP::Get::Analysis qw{ get_overall_analysis_type print_program };
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Script::Utils qw{ help };
use MIP::Set::Contigs qw{ set_contigs };
use MIP::Set::Parameter qw{ set_dynamic_parameter };
use MIP::Update::Contigs qw{ update_contigs_for_run };
use MIP::Update::Programs
  qw{ update_program_mode_with_dry_run_all update_program_mode update_prioritize_flag };

## Recipes
use MIP::Recipes::Bamcalibrationblock qw{ analysis_bamcalibrationblock };
use MIP::Recipes::Bedtools_genomecov qw{ analysis_bedtools_genomecov };
use MIP::Recipes::Bwa_mem qw{ analysis_bwa_mem };
use MIP::Recipes::Build::Human_genome_prerequisites
  qw{ build_human_genome_prerequisites };
use MIP::Recipes::Chanjo_sex_check qw{ analysis_chanjo_sex_check };
use MIP::Recipes::Cnvnator qw{ analysis_cnvnator };
use MIP::Recipes::Fastqc qw{ analysis_fastqc };
use MIP::Recipes::Gatk_baserecalibration qw{ analysis_gatk_baserecalibration };
use MIP::Recipes::Gatk_haplotypecaller qw{ analysis_gatk_haplotypecaller };
use MIP::Recipes::Gatk_realigner qw{ analysis_gatk_realigner };
use MIP::Recipes::Manta qw{ analysis_manta };
use MIP::Recipes::Markduplicates qw{ analysis_markduplicates };
use MIP::Recipes::Picardtools_collecthsmetrics
  qw{ analysis_picardtools_collecthsmetrics };
use MIP::Recipes::Picardtools_collectmultiplemetrics
  qw{ analysis_picardtools_collectmultiplemetrics };
use MIP::Recipes::Picardtools_genotypeconcordance
  qw{ analysis_picardtools_genotypeconcordance };
use MIP::Recipes::Picardtools_mergesamfiles
  qw{ analysis_picardtools_mergesamfiles };
use MIP::Recipes::Pipeline::Wts qw{ pipeline_wts };
use MIP::Recipes::Rcoverageplots qw{ analysis_rcoverageplots };
use MIP::Recipes::Sambamba_depth qw{ analysis_sambamba_depth };
use MIP::Recipes::Split_fastq_file qw{ analysis_split_fastq_file };
use MIP::Recipes::Tiddit qw{ analysis_tiddit };
use MIP::Recipes::Variant_integrity qw{ analysis_variant_integrity };
use MIP::Recipes::Vep qw{ analysis_vep analysis_vep_rio analysis_vep_sv };
use MIP::Recipes::Vt_core qw{ analysis_vt_core analysis_vt_core_rio};

our $USAGE = build_usage( {} );

BEGIN {

    require MIP::Check::Modules;

    my @modules = (
        qw{ YAML Path::Iterator::Rule
          List::Util Log::Log4perl
          MIP::File::Format::Yaml MIP::Set::File
          MIP::Log::MIP_log4perl  MIP::Script::Utils }
    );

    ## Evaluate that all modules required are installed
    check_perl_modules(
        {
            modules_ref  => \@modules,
            program_name => $PROGRAM_NAME,
        }
    );
}

## Constants
Readonly my $DOT       => q{.};
Readonly my $EMPTY_STR => q{};
Readonly my $NEWLINE   => qq{\n};
Readonly my $SPACE     => q{ };
Readonly my $TAB       => qq{\t};

#### Script parameters

# Holds all active parameters after the value has been set
my %active_parameter;

# Holds all set parameters info after add_to_active_parameter
my @broadcasts;

## Add date_time_stamp for later use in log and qc_metrics yaml file
my $date_time       = localtime;
my $date_time_stamp = $date_time->datetime;
my $date            = $date_time->ymd;

# Catches script name and removes ending
my $script = fileparse( basename( $PROGRAM_NAME, $DOT . q{pl} ) );
my $definitions_file =
  catfile( $Bin, qw{ definitions define_parameters.yaml } );
chomp( $date_time_stamp, $date, $script );

#### Set program parameters

### %parameter holds all parameters for MIP
## Loads a YAML file into an arbitrary hash and returns it.
my %parameter = load_yaml( { yaml_file => $definitions_file, } );

### To add/write parameters in the correct order
## Adds the order of first level keys from yaml file to array
my @order_parameters = order_parameter_names(
    {
        file_path => $definitions_file,
    }
);

## Load mandatory keys and values for parameters
my %mandatory_key = load_yaml(
    {
        yaml_file =>
          catfile( $Bin, qw{ definitions mandatory_parameter_keys.yaml } ),
    }
);

## Load non mandatory keys and values for parameters
my %non_mandatory_key = load_yaml(
    {
        yaml_file =>
          catfile( $Bin, qw{ definitions non_mandatory_parameter_keys.yaml } ),

    }
);

## Eval parameter hash
my $error_msg = check_parameter_hash(
    {
        parameter_href         => \%parameter,
        mandatory_key_href     => \%mandatory_key,
        non_mandatory_key_href => \%non_mandatory_key,
        file_path              => $definitions_file,
    }
);
## Broadcast error in parameter key or values
if ($error_msg) {

    croak($error_msg);
}

# Set MIP version
our $VERSION = 'v5.0.11';

## Directories, files, job_ids and sample_info
my ( %infile, %indir_path, %infile_lane_prefix, %lane,
    %infile_both_strands_prefix, %job_id, %sample_info );

#### Staging/Sanity Check Area

my %file_info = (
    bwa_build_reference => $EMPTY_STR,
    exome_target_bed =>
      [qw{ .infile_list .pad100.infile_list .pad100.interval_list }],

    # BWA human genome reference file endings
    bwa_build_reference_file_endings => [qw{ .amb .ann .bwt .pac .sa }],

    # Human genome meta files
    human_genome_reference_file_endings => [qw{ .dict .fai }],
);

## Enables cmd "mip" to print usage help
if ( not @ARGV ) {

    help(
        {
            USAGE     => $USAGE,
            exit_code => 0,
        }
    );
}

###User Options
GetOptions(
    q{ifd|infile_dirs:s}   => \%{ $parameter{infile_dirs}{value} },
    q{rd|reference_dir:s}  => \$parameter{reference_dir}{value},
    q{p|project_id:s}      => \$parameter{project_id}{value},
    q{s|sample_ids:s}      => \@{ $parameter{sample_ids}{value} },
    q{odd|outdata_dir:s}   => \$parameter{outdata_dir}{value},
    q{osd|outscript_dir:s} => \$parameter{outscript_dir}{value},
    q{f|family_id:s}       => \$parameter{family_id}{value},
    q{sck|supported_capture_kit:s} =>
      \%{ $parameter{supported_capture_kit}{value} },
    q{dnr|decompose_normalize_references:s} =>
      \@{ $parameter{decompose_normalize_references}{value} },
    q{ped|pedigree_file:s} => \$parameter{pedigree_file}{value},
    q{hgr|human_genome_reference:s} =>
      \$parameter{human_genome_reference}{value},
    q{al|outaligner_dir:s}    => \$parameter{outaligner_dir}{value},
    q{at|analysis_type:s}     => \%{ $parameter{analysis_type}{value} },
    q{pl|platform:s}          => \$parameter{platform}{value},
    q{ec|expected_coverage:s} => \%{ $parameter{expected_coverage}{value} },
    q{c|config_file:s}        => \$parameter{config_file}{value},
    q{ccp|cluster_constant_path:s} => \$parameter{cluster_constant_path}{value},
    q{acp|analysis_constant_path:s} =>
      \$parameter{analysis_constant_path}{value},
    q{cfa|config_file_analysis:s} => \$parameter{config_file_analysis}{value},
    q{sif|sample_info_file:s}     => \$parameter{sample_info_file}{value},
    q{dra|dry_run_all=i}          => \$parameter{dry_run_all}{value},
    q{jul|java_use_large_pages=n} => \$parameter{java_use_large_pages}{value},
    q{ges|genomic_set:s}          => \$parameter{genomic_set}{value},
    q{rio|reduce_io=n}            => \$parameter{reduce_io}{value},
    q{riu|replace_iupac=n}        => \$parameter{replace_iupac}{value},
    q{ppm|print_program_mode=n}   => \$parameter{print_program_mode}{value},
    q{pp|print_program}           => sub {
        GetOptions( q{ppm|print_program_mode=n} =>
              \$parameter{print_program_mode}{value} )
          ;    #Force ppm to be read before function call
        print_program(
            {
                parameter_href     => \%parameter,
                print_program_mode => $parameter{print_program_mode}{value},
            }
        );
        exit;
    },
    q{l|log_file:s} => \$parameter{log_file}{value},
    q{h|help}       => sub { say STDOUT $USAGE; exit; },
    q{v|version}    => sub {
        say STDOUT $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $VERSION,
          $NEWLINE;
        exit;
    },

    #### Bash
    q{bse|bash_set_errexit:s}   => \$parameter{bash_set_errexit}{value},
    q{bsu|bash_set_nounset:s}   => \$parameter{bash_set_nounset}{value},
    q{bsp|bash_set_pipefail:s}  => \$parameter{bash_set_pipefail}{value},
    q{em|email:s}               => \$parameter{email}{value},
    q{emt|email_types:s}        => \@{ $parameter{email_types}{value} },
    q{mcn|module_core_number:s} => \%{ $parameter{module_core_number}{value} },
    q{mot|module_time:s}        => \%{ $parameter{module_time}{value} },
    q{mcn|max_cores_per_node=n} => \$parameter{max_cores_per_node}{value},
    q{nrm|node_ram_memory=n}    => \$parameter{node_ram_memory}{value},
    q{tmd|temp_directory:s}     => \$parameter{temp_directory}{value},
    q{qos|slurm_quality_of_service=s} =>
      \$parameter{slurm_quality_of_service}{value},
    q{sen|source_environment_commands=s{,}} =>
      \@{ $parameter{source_environment_commands}{value} },
    q{psfq|psplit_fastq_file=n} => \$parameter{psplit_fastq_file}{value},
    q{sfqrdb|split_fastq_file_read_batch=n} =>
      \$parameter{split_fastq_file_read_batch}{value},
    q{pgz|pgzip_fastq=n}         => \$parameter{pgzip_fastq}{value},
    q{pfqc|pfastqc=n}            => \$parameter{pfastqc}{value},
    q{pmem|pbwa_mem=n}           => \$parameter{pbwa_mem}{value},
    q{memhla|bwa_mem_hla=n}      => \$parameter{bwa_mem_hla}{value},
    q{memrdb|bwa_mem_rapid_db:s} => \$parameter{bwa_mem_rapid_db}{value},
    q{memcrm|bwa_mem_cram=n}     => \$parameter{bwa_mem_cram}{value},
    q{memsts|bwa_mem_bamstats=n} => \$parameter{bwa_mem_bamstats}{value},
    q{memssm|bwa_sambamba_sort_memory_limit:s} =>
      \$parameter{bwa_sambamba_sort_memory_limit}{value},
    q{ptp|picardtools_path:s} => \$parameter{picardtools_path}{value},
    q{pptm|ppicardtools_mergesamfiles=n} =>
      \$parameter{ppicardtools_mergesamfiles}{value},
    q{pmd|pmarkduplicates=n} => \$parameter{pmarkduplicates}{value},
    q{mdpmd|markduplicates_picardtools_markduplicates=n} =>
      \$parameter{markduplicates_picardtools_markduplicates}{value},
    q{mdsmd|markduplicates_sambamba_markdup=n} =>
      \$parameter{markduplicates_sambamba_markdup}{value},
    q{mdshts|markduplicates_sambamba_markdup_hash_table_size=n} =>
      \$parameter{markduplicates_sambamba_markdup_hash_table_size}{value},
    q{mdsols|markduplicates_sambamba_markdup_overflow_list_size=n} =>
      \$parameter{markduplicates_sambamba_markdup_overflow_list_size}{value},
    q{mdsibs|markduplicates_sambamba_markdup_io_buffer_size=n} =>
      \$parameter{markduplicates_sambamba_markdup_io_buffer_size}{value},
    q{pchs|pchanjo_sexcheck=n} => \$parameter{pchanjo_sexcheck}{value},
    q{chslle|chanjo_sexcheck_log_level:s} =>
      \$parameter{chanjo_sexcheck_log_level}{value},
    q{psdt|psambamba_depth=n}       => \$parameter{psambamba_depth}{value},
    q{sdtmod|sambamba_depth_mode:s} => \$parameter{sambamba_depth_mode}{value},
    q{sdtcut|sambamba_depth_cutoffs:s} =>
      \@{ $parameter{sambamba_depth_cutoffs}{value} },
    q{sdtbed|sambamba_depth_bed:s} => \$parameter{sambamba_depth_bed}{value},
    q{sdtbaq|sambamba_depth_base_quality=n} =>
      \$parameter{sambamba_depth_base_quality}{value},
    q{sdtmaq|sambamba_depth_mapping_quality=n} =>
      \$parameter{sambamba_depth_mapping_quality}{value},
    q{sdtndu|sambamba_depth_noduplicates=n} =>
      \$parameter{sambamba_depth_noduplicates}{value},
    q{sdtfqc|sambamba_depth_quality_control=n} =>
      \$parameter{sambamba_depth_quality_control}{value},
    q{pbgc|pbedtools_genomecov=n} => \$parameter{pbedtools_genomecov}{value},
    q{bgcmc|bedtools_genomecov_max_coverage=n} =>
      \$parameter{bedtools_genomecov_max_coverage}{value},
    q{pptcmm|ppicardtools_collectmultiplemetrics=n} =>
      \$parameter{ppicardtools_collectmultiplemetrics}{value},
    q{pptchs|ppicardtools_collecthsmetrics=n} =>
      \$parameter{ppicardtools_collecthsmetrics}{value},
    q{extb|exome_target_bed=s}     => \%{ $parameter{exome_target_bed}{value} },
    q{prcp|prcovplots=n}           => \$parameter{prcovplots}{value},
    q{pcnv|pcnvnator=n}            => \$parameter{pcnvnator}{value},
    q{cnvhbs|cnv_bin_size=n}       => \$parameter{cnv_bin_size}{value},
    q{cnvrld|cnv_root_ld_lib=s}    => \$parameter{cnv_root_ld_lib}{value},
    q{pdelc|pdelly_call=n}         => \$parameter{pdelly_call}{value},
    q{pdel|pdelly_reformat=n}      => \$parameter{pdelly_reformat}{value},
    q{deltyp|delly_types:s}        => \@{ $parameter{delly_types}{value} },
    q{delexc|delly_exclude_file:s} => \$parameter{delly_exclude_file}{value},
    q{pmna|pmanta=n}               => \$parameter{pmanta}{value},
    q{ptid|ptiddit=n}              => \$parameter{ptiddit}{value},
    q{tidmsp|tiddit_minimum_number_supporting_pairs=n} =>
      \$parameter{tiddit_minimum_number_supporting_pairs}{value},
    q{psvc|psv_combinevariantcallsets=n} =>
      \$parameter{psv_combinevariantcallsets}{value},
    q{svcvtd|sv_vt_decompose=n} => \$parameter{sv_vt_decompose}{value},
    q{svsvdbmp|sv_svdb_merge_prioritize:s} =>
      \$parameter{sv_svdb_merge_prioritize}{value},
    q{svcbtv|sv_bcftools_view_filter=n} =>
      \$parameter{sv_bcftools_view_filter}{value},
    q{svcdbq|sv_svdb_query=n} => \$parameter{sv_svdb_query}{value},
    q{svcdbqd|sv_svdb_query_db_files:s} =>
      \%{ $parameter{sv_svdb_query_db_files}{value} },
    q{svcvan|sv_vcfanno=n}        => \$parameter{sv_vcfanno}{value},
    q{svcval|sv_vcfanno_lua:s}    => \$parameter{sv_vcfanno_lua}{value},
    q{svcvac|sv_vcfanno_config:s} => \$parameter{sv_vcfanno_config}{value},
    q{svcvacf|sv_vcfanno_config_file:s} =>
      \$parameter{sv_vcfanno_config_file}{value},
    q{svcvah|sv_vcfannotation_header_lines_file:s} =>
      \$parameter{sv_vcfannotation_header_lines_file}{value},
    q{svcgmf|sv_genmod_filter=n} => \$parameter{sv_genmod_filter}{value},
    q{svcgfr|sv_genmod_filter_1000g:s} =>
      \$parameter{sv_genmod_filter_1000g}{value},
    q{svcgft|sv_genmod_filter_threshold:s} =>
      \$parameter{sv_genmod_filter_threshold}{value},
    q{svcbcf|sv_combinevariantcallsets_bcf_file=n} =>
      \$parameter{sv_combinevariantcallsets_bcf_file}{value},
    q{psvv|psv_varianteffectpredictor=n} =>
      \$parameter{psv_varianteffectpredictor}{value},
    q{svvepf|sv_vep_features:s} => \@{ $parameter{sv_vep_features}{value} },
    q{svvepl|sv_vep_plugins:s}  => \@{ $parameter{sv_vep_plugins}{value} },
    q{psvvcp|psv_vcfparser=n}   => \$parameter{psv_vcfparser}{value},
    q{svvcpvt|sv_vcfparser_vep_transcripts=n} =>
      \$parameter{sv_vcfparser_vep_transcripts}{value},
    q{svvcppg|sv_vcfparser_per_gene=n} =>
      \$parameter{sv_vcfparser_per_gene}{value},
    q{svvcprff|sv_vcfparser_range_feature_file:s} =>
      \$parameter{sv_vcfparser_range_feature_file}{value},
    q{svvcprfa|sv_vcfparser_range_feature_annotation_columns:s} =>
      \@{ $parameter{sv_vcfparser_range_feature_annotation_columns}{value} },
    q{svvcpsf|sv_vcfparser_select_file:s} =>
      \$parameter{sv_vcfparser_select_file}{value},
    q{svvcpsfm|sv_vcfparser_select_file_matching_column=n} =>
      \$parameter{sv_vcfparser_select_file_matching_column}{value},
    q{svvcpsfa|sv_vcfparser_select_feature_annotation_columns:s} =>
      \@{ $parameter{sv_vcfparser_select_feature_annotation_columns}{value} },
    q{psvr|psv_rankvariant=n} => \$parameter{psv_rankvariant}{value},
    q{svravanr|sv_genmod_annotate_regions:n} =>
      \$parameter{sv_genmod_annotate_regions}{value},
    q{svravgft|sv_genmod_models_family_type:s} =>
      \$parameter{sv_genmod_models_family_type}{value},
    q{svravrpf|sv_genmod_models_reduced_penetrance_file:s} =>
      \$parameter{sv_genmod_models_reduced_penetrance_file}{value},
    q{svravwg|sv_genmod_models_whole_gene=n} =>
      \$parameter{sv_genmod_models_whole_gene}{value},
    q{svravrm|sv_rank_model_file:s} => \$parameter{sv_rank_model_file}{value},
    q{psvre|psv_reformat=n}         => \$parameter{psv_reformat}{value},
    q{svrevbf|sv_rankvariant_binary_file=n} =>
      \$parameter{sv_rankvariant_binary_file}{value},
    q{svrergf|sv_reformat_remove_genes_file:s} =>
      \$parameter{sv_reformat_remove_genes_file}{value},
    q{psmp|psamtools_mpileup=n} => \$parameter{psamtools_mpileup}{value},
    q{pfrb|pfreebayes=n}        => \$parameter{pfreebayes}{value},
    q{gtp|gatk_path:s}          => \$parameter{gatk_path}{value},
    q{gll|gatk_logging_level:s} => \$parameter{gatk_logging_level}{value},
    q{gbdv|gatk_bundle_download_version:s} =>
      \$parameter{gatk_bundle_download_version}{value},
    q{gdco|gatk_downsample_to_coverage=n} =>
      \$parameter{gatk_downsample_to_coverage}{value},
    q{gdai|gatk_disable_auto_index_and_file_lock=n} =>
      \$parameter{gatk_disable_auto_index_and_file_lock}{value},
    q{pgra|pgatk_realigner=n} => \$parameter{pgatk_realigner}{value},
    q{graks|gatk_realigner_indel_known_sites:s} =>
      \@{ $parameter{gatk_realigner_indel_known_sites}{value} },
    q{pgbr|pgatk_baserecalibration=n} =>
      \$parameter{pgatk_baserecalibration}{value},
    q{gbrcov|gatk_baserecalibration_covariates:s} =>
      \@{ $parameter{gatk_baserecalibration_covariates}{value} },
    q{gbrkst|gatk_baserecalibration_known_sites:s} =>
      \@{ $parameter{gatk_baserecalibration_known_sites}{value} },
    q{gbrrf|gatk_baserecalibration_read_filters=s} =>
      \@{ $parameter{gatk_baserecalibration_read_filters}{value} },
    q{gbrdiq|gatk_baserecalibration_disable_indel_qual=n} =>
      \$parameter{gatk_baserecalibration_disable_indel_qual}{value},
    q{gbrsqq|gatk_baserecalibration_static_quantized_quals:s} =>
      \@{ $parameter{gatk_baserecalibration_static_quantized_quals}{value} },
    q{pghc|pgatk_haplotypecaller=n} =>
      \$parameter{pgatk_haplotypecaller}{value},
    q{ghcann|gatk_haplotypecaller_annotation:s} =>
      \@{ $parameter{gatk_haplotypecaller_annotation}{value} },
    q{ghckse|gatk_haplotypecaller_snp_known_set:s} =>
      \$parameter{gatk_haplotypecaller_snp_known_set}{value},
    q{ghcscb|gatk_haplotypecaller_no_soft_clipped_bases=n} =>
      \$parameter{gatk_haplotypecaller_no_soft_clipped_bases}{value},
    q{ghcpim|gatk_haplotypecaller_pcr_indel_model:s} =>
      \$parameter{gatk_haplotypecaller_pcr_indel_model}{value},
    q{pggt|pgatk_genotypegvcfs=n} => \$parameter{pgatk_genotypegvcfs}{value},
    q{ggtgrl|gatk_genotypegvcfs_ref_gvcf:s} =>
      \$parameter{gatk_genotypegvcfs_ref_gvcf}{value},
    q{ggtals|gatk_genotypegvcfs_all_sites=n} =>
      \$parameter{gatk_genotypegvcfs_all_sites}{value},
    q{ggbcf|gatk_concatenate_genotypegvcfs_bcf_file=n} =>
      \$parameter{gatk_concatenate_genotypegvcfs_bcf_file}{value},
    q{pgvr|pgatk_variantrecalibration=n} =>
      \$parameter{pgatk_variantrecalibration}{value},
    q{gvrann|gatk_variantrecalibration_annotations:s} =>
      \@{ $parameter{gatk_variantrecalibration_annotations}{value} },
    q{gvrres|gatk_variantrecalibration_resource_snv=s} =>
      \%{ $parameter{gatk_variantrecalibration_resource_snv}{value} },
    q{gvrrei|gatk_variantrecalibration_resource_indel=s} =>
      \%{ $parameter{gatk_variantrecalibration_resource_indel}{value} },
    q{gvrstf|gatk_variantrecalibration_snv_tsfilter_level:s} =>
      \$parameter{gatk_variantrecalibration_snv_tsfilter_level}{value},
    q{gvritf|gatk_variantrecalibration_indel_tsfilter_level:s} =>
      \$parameter{gatk_variantrecalibration_indel_tsfilter_level}{value},
    q{gvrdpa|gatk_variantrecalibration_dp_annotation=n} =>
      \$parameter{gatk_variantrecalibration_dp_annotation}{value},
    q{gvrsmg|gatk_variantrecalibration_snv_max_gaussians=n} =>
      \$parameter{gatk_variantrecalibration_snv_max_gaussians}{value},
    q{gvrimg|gatk_variantrecalibration_indel_max_gaussians=n} =>
      \$parameter{gatk_variantrecalibration_indel_max_gaussians}{value},
    q{gvrevf|gatk_variantrecalibration_exclude_nonvariants_file=n} =>
      \$parameter{gatk_variantrecalibration_exclude_nonvariants_file}{value},
    q{gvrbcf|gatk_variantrecalibration_bcf_file=n} =>
      \$parameter{gatk_variantrecalibration_bcf_file}{value},
    q{gcgpss|gatk_calculategenotypeposteriors_support_set:s} =>
      \$parameter{gatk_calculategenotypeposteriors_support_set}{value},
    q{pgcv|pgatk_combinevariantcallsets=n} =>
      \$parameter{pgatk_combinevariantcallsets}{value},
    q{gcvgmo|gatk_combinevariants_genotype_merge_option:s} =>
      \$parameter{gatk_combinevariants_genotype_merge_option}{value},
    q{gcvpc|gatk_combinevariants_prioritize_caller:s} =>
      \$parameter{gatk_combinevariants_prioritize_caller}{value},
    q{gcvbcf|gatk_combinevariantcallsets_bcf_file=n} =>
      \$parameter{gatk_combinevariantcallsets_bcf_file}{value},
    q{pgvea|pgatk_variantevalall=n} => \$parameter{pgatk_variantevalall}{value},
    q{pgvee|pgatk_variantevalexome=n} =>
      \$parameter{pgatk_variantevalexome}{value},
    q{gveedbs|gatk_varianteval_dbsnp:s} =>
      \$parameter{gatk_varianteval_dbsnp}{value},
    q{gveedbg|gatk_varianteval_gold:s} =>
      \$parameter{gatk_varianteval_gold}{value},
    q{ppvab|pprepareforvariantannotationblock=n} =>
      \$parameter{pprepareforvariantannotationblock}{value},
    q{prhc|prhocall=n} => \$parameter{prhocall}{value},
    q{rhcf|rhocall_frequency_file:s} =>
      \$parameter{rhocall_frequency_file}{value},
    q{pvt|pvt=n}             => \$parameter{pvt}{value},
    q{vtddec|vt_decompose=n} => \$parameter{vt_decompose}{value},
    q{vtdnor|vt_normalize=n} => \$parameter{vt_normalize}{value},
    q{vtunq|vt_uniq=n}       => \$parameter{vt_uniq}{value},
    q{vtmaa|vt_missing_alt_allele=n} =>
      \$parameter{vt_missing_alt_allele}{value},
    q{vtgmf|vt_genmod_filter=n} => \$parameter{vt_genmod_filter}{value},
    q{vtgfr|vt_genmod_filter_1000g:s} =>
      \$parameter{vt_genmod_filter_1000g}{value},
    q{vtmaf|vt_genmod_filter_max_af=n} =>
      \$parameter{vt_genmod_filter_max_af}{value},
    q{vtgft|vt_genmod_filter_threshold:s} =>
      \$parameter{vt_genmod_filter_threshold}{value},
    q{pvep|pvarianteffectpredictor=n} =>
      \$parameter{pvarianteffectpredictor}{value},
    q{vepp|vep_directory_path:s}  => \$parameter{vep_directory_path}{value},
    q{vepc|vep_directory_cache:s} => \$parameter{vep_directory_cache}{value},
    q{vepr|vep_reference:n}       => \$parameter{vep_reference}{value},
    q{vepf|vep_features:s}        => \@{ $parameter{vep_features}{value} },
    q{veppl|vep_plugins:s}        => \@{ $parameter{vep_plugins}{value} },
    q{pvcp|pvcfparser=n}          => \$parameter{pvcfparser}{value},
    q{vcpvt|vcfparser_vep_transcripts=n} =>
      \$parameter{vcfparser_vep_transcripts}{value},
    q{vcprff|vcfparser_range_feature_file:s} =>
      \$parameter{vcfparser_range_feature_file}{value},
    q{vcprfa|vcfparser_range_feature_annotation_columns:s} =>
      \@{ $parameter{vcfparser_range_feature_annotation_columns}{value} },
    q{vcpsf|vcfparser_select_file:s} =>
      \$parameter{vcfparser_select_file}{value},
    q{vcpsfm|vcfparser_select_file_matching_column=n} =>
      \$parameter{vcfparser_select_file_matching_column}{value},
    q{vcpsfa|vcfparser_select_feature_annotation_columns:s} =>
      \@{ $parameter{vcfparser_select_feature_annotation_columns}{value} },
    q{psne|psnpeff=n}      => \$parameter{psnpeff}{value},
    q{snep|snpeff_path:s}  => \$parameter{snpeff_path}{value},
    q{sneann|snpeff_ann=n} => \$parameter{snpeff_ann}{value},
    q{snegbv|snpeff_genome_build_version:s} =>
      \$parameter{snpeff_genome_build_version}{value},
    q{snesaf|snpsift_annotation_files=s} =>
      \%{ $parameter{snpsift_annotation_files}{value} },
    q{snesaoi|snpsift_annotation_outinfo_key=s} =>
      \%{ $parameter{snpsift_annotation_outinfo_key}{value} },
    q{snesdbnsfp|snpsift_dbnsfp_file:s} =>
      \$parameter{snpsift_dbnsfp_file}{value},
    q{snesdbnsfpa|snpsift_dbnsfp_annotations:s} =>
      \@{ $parameter{snpsift_dbnsfp_annotations}{value} },
    q{prav|prankvariant=n} => \$parameter{prankvariant}{value},
    q{ravgft|genmod_models_family_type:s} =>
      \$parameter{genmod_models_family_type}{value},
    q{ravanr|genmod_annotate_regions:n} =>
      \$parameter{genmod_annotate_regions}{value},
    q{ravcad|genmod_annotate_cadd_files:s} =>
      \@{ $parameter{genmod_annotate_cadd_files}{value} },
    q{ravspi|genmod_annotate_spidex_file:s} =>
      \$parameter{genmod_annotate_spidex_file}{value},
    q{ravwg|genmod_models_whole_gene=n} =>
      \$parameter{genmod_models_whole_gene}{value},
    q{ravrpf|genmod_models_reduced_penetrance_file:s} =>
      \$parameter{genmod_models_reduced_penetrance_file}{value},
    q{ravrm|rank_model_file:s} => \$parameter{rank_model_file}{value},
    q{pevab|pendvariantannotationblock=n} =>
      \$parameter{pendvariantannotationblock}{value},
    q{evabrgf|endvariantannotationblock_remove_genes_file:s} =>
      \$parameter{endvariantannotationblock_remove_genes_file}{value},
    q{ravbf|rankvariant_binary_file=n} =>
      \$parameter{rankvariant_binary_file}{value},
    q{pped|ppeddy=n}             => \$parameter{ppeddy}{value},
    q{pplink|pplink=n}           => \$parameter{pplink}{value},
    q{pvai|pvariant_integrity=n} => \$parameter{pvariant_integrity}{value},
    q{pevl|pevaluation=n}        => \$parameter{pevaluation}{value},
    q{evlnid|nist_id:s}          => \$parameter{nist_id}{value},
    q{evlnhc|nist_high_confidence_call_set:s} =>
      \$parameter{nist_high_confidence_call_set}{value},
    q{evlnil|nist_high_confidence_call_set_bed:s} =>
      \$parameter{nist_high_confidence_call_set_bed}{value},
    q{pqcc|pqccollect=n} => \$parameter{pqccollect}{value},
    q{qccsi|qccollect_sampleinfo_file:s} =>
      \$parameter{qccollect_sampleinfo_file}{value},
    q{qccref|qccollect_regexp_file:s} =>
      \$parameter{qccollect_regexp_file}{value},
    q{qccske|qccollect_skip_evaluation} =>
      \$parameter{qccollect_skip_evaluation}{value},
    q{pmqc|pmultiqc=n} => \$parameter{pmultiqc}{value},
    q{prem|premoveredundantfiles=n} =>
      \$parameter{premoveredundantfiles}{value},
    q{pars|panalysisrunstatus=n} => \$parameter{panalysisrunstatus}{value},
    q{psac|psacct=n}             => \$parameter{psacct}{value},
    q{sacfrf|sacct_format_fields:s} =>
      \@{ $parameter{sacct_format_fields}{value} },
  )
  or help(
    {
        USAGE     => $USAGE,
        exit_code => 1,
    }
  );

## Change relative path to absolute path for certain parameters
update_to_absolute_path( { parameter_href => \%parameter, } );

##Special case:Enable/activate MIP. Cannot be changed from cmd or config
$active_parameter{mip} = $parameter{mip}{default};

### Config file
## Input from cmd
if ( defined $parameter{config_file}{value} ) {

    ## Loads a YAML file into an arbitrary hash and returns it.
    %active_parameter =
      load_yaml( { yaml_file => $parameter{config_file}{value}, } );

    ## Special case: Enable/activate MIP. Cannot be changed from cmd or config
    $active_parameter{mip} = $parameter{mip}{default};

    ## Special case: Required when turning of vcfParser to know how many files should be analysed (.select.vcf or just .vcf)
    $active_parameter{vcfparser_outfile_count} =
      $parameter{vcfparser_outfile_count}{default};

    ## Compare keys from config and definitions file
    check_config_vs_definition_file(
        {
            reference_href  => \%active_parameter,
            comparison_href => \%parameter,
        }
    );

    my @config_dynamic_parameters =
      (qw(cluster_constant_path analysis_constant_path outaligner_dir));

    ## Replace config parameter with cmd info for config dynamic parameter
    replace_config_parameters_with_cmd_info(
        {
            parameter_href        => \%parameter,
            active_parameter_href => \%active_parameter,
            parameter_names_ref   => \@config_dynamic_parameters,
        }
    );

    ## Loop through all parameters and update info
  PARAMETER:
    foreach my $order_parameter_element (@order_parameters) {

        ## Updates the config file to particular user/cluster for entries following specifications. Leaves other entries untouched.
        update_config_file(
            {
                active_parameter_href => \%active_parameter,
                parameter_name_ref    => \$order_parameter_element,
                family_id_ref         => \$parameter{family_id}{value},
            }
        );
    }

    ## Remove previous analysis specific info not relevant for current run e.g. log file, sample_ids which are read from pedigree or cmd
    my @remove_keys = (qw{ log_file sample_ids });

  KEY:
    foreach my $key (@remove_keys) {

        delete $active_parameter{$key};
    }
}

###Populate active_parameters{parameter_name} => 'Value'
foreach my $order_parameter_element (@order_parameters) {

    ## Type of variables: mip, path or program/program_parameters each is handled in the add_to_active_parameter subroutine.
    ## Checks and sets user input or default values to active_parameters
    add_to_active_parameter(
        {
            parameter_href        => \%parameter,
            active_parameter_href => \%active_parameter,
            sample_info_href      => \%sample_info,
            file_info_href        => \%file_info,
            broadcasts_ref        => \@broadcasts,
            associated_programs_ref =>
              \@{ $parameter{$order_parameter_element}{associated_program} },
            parameter_name => $order_parameter_element,
        }
    );

    ## Special case for parameters that are dependent on other parameters values
    if ( $order_parameter_element eq 'outdata_dir' )
    { #Set defaults depending on $active_parameter{outdata_dir} value that now has been set

        $parameter{sample_info_file}{default} = catfile(
            $active_parameter{outdata_dir},
            $active_parameter{family_id},
            $active_parameter{family_id} . '_qc_sample_info.yaml'
        );

        ## Set the default Log4perl file using supplied dynamic parameters.
        $parameter{log_file}{default} = deafult_log4perl_file(
            {
                active_parameter_href => \%active_parameter,
                cmd_input_ref         => \$parameter{log_file}{value},
                script_ref            => \$script,
                date_ref              => \$date,
                date_time_stamp_ref   => \$date_time_stamp,
            }
        );

        $parameter{qccollect_sampleinfo_file}{default} =
          $parameter{sample_info_file}{default};
    }
    if ( $order_parameter_element eq q{log_file} ) {

        ## Creates log object
        my $log = initiate_logger(
            {
                file_path => $active_parameter{log_file},
                log_name  => q{MIP},
            }
        );
    }
    if ( $order_parameter_element eq 'pedigree_file' )
    {    #Write QC for only pedigree data used in analysis

        if ( defined( $active_parameter{pedigree_file} ) ) {

            ## Retrieve logger object now that log_file has been set
            my $log = Log::Log4perl->get_logger(q{MIP});

            make_path(
                catdir(
                    $active_parameter{outdata_dir},
                    $active_parameter{family_id}
                )
            );    #Create family directory
            my $yaml_file = catfile(
                $active_parameter{outdata_dir},
                $active_parameter{family_id},
                'qc_pedigree.yaml'
            );

            ## Writes a YAML hash to file
            write_yaml(
                {
                    yaml_href      => \%sample_info,
                    yaml_file_path => $yaml_file,
                }
            );
            $log->info( 'Wrote: ' . $yaml_file, "\n" );

            ## Removes all elements at hash third level except keys in allowed_entries
            remove_pedigree_elements( { hash_ref => \%sample_info, } );
        }
    }
    if ( $order_parameter_element eq q{analysis_type} ) {

        ## Detect if all samples has the same sequencing type and return consensus if reached
        $parameter{dynamic_parameter}{consensus_analysis_type} =
          get_overall_analysis_type(
            { analysis_type_href => \%{ $active_parameter{analysis_type} }, } );
    }
}

## Retrieve logger object now that log_file has been set
my $log = Log::Log4perl->get_logger(q{MIP});

###Checks

## Check Existance of files and directories
foreach my $parameter_name ( keys %parameter ) {

    if ( exists( $parameter{$parameter_name}{exists_check} ) ) {

        check_parameter_files(
            {
                parameter_href        => \%parameter,
                active_parameter_href => \%active_parameter,
                sample_info_href      => \%sample_info,
                file_info_href        => \%file_info,
                broadcasts_ref        => \@broadcasts,
                associated_programs_ref =>
                  \@{ $parameter{$parameter_name}{associated_program} },
                parameter_name => $parameter_name,
                parameter_exists_check =>
                  $parameter{$parameter_name}{exists_check},
            }
        );
    }
}

## Detect family constellation based on pedigree file
$parameter{dynamic_parameter}{trio} = detect_trio(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
);

## Detect number of founders (i.e. parents ) based on pedigree file
detect_founders(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
);

## Check email adress format
if ( exists( $active_parameter{email} ) ) {    #Allow no malformed email adress

    check_email_address( { email_ref => \$active_parameter{email}, } );
}

if (
    not check_allowed_temp_directory(
        { temp_directory => $active_parameter{temp_directory}, }
    )
  )
{

    $log->fatal( q{'--temp_directory }
          . $active_parameter{temp_directory}
          . q{' is not allowed because MIP will remove the temp directory after processing.}
          . "\n" );
    exit 1;
}

## Parameters that have keys as MIP program names
my @parameter_keys_to_check = (qw(module_time module_core_number));
foreach my $parameter_name (@parameter_keys_to_check) {

    ## Test if key from query hash exists truth hash
    check_key_exists_in_hash(
        {
            truth_href     => \%parameter,
            query_href     => \%{ $active_parameter{$parameter_name} },
            parameter_name => $parameter_name,
        }
    );
}

## Check that the module core number do not exceed the maximum per node
foreach my $program_name ( keys %{ $active_parameter{module_core_number} } ) {

    ## Limit number of cores requested to the maximum number of cores available per node
    $active_parameter{module_core_number}{$program_name} =
      check_max_core_number(
        {
            max_cores_per_node => $active_parameter{max_cores_per_node},
            core_number_requested =>
              $active_parameter{module_core_number}{$program_name},
        }
      );
}

## Parameters that have elements as MIP program names
my @parameter_elements_to_check =
  (qw(associated_program decompose_normalize_references));
foreach my $parameter_name (@parameter_elements_to_check) {

    ## Test if element from query array exists truth hash
    check_element_exists_in_hash(
        {
            truth_href     => \%parameter,
            queryies       => \@{ $active_parameter{$parameter_name} },
            parameter_name => $parameter_name,
        }
    );
}

## Parameters with key(s) that have elements as MIP program names
my @parameter_key_to_check = qw(associated_program);
PARAMETER:
foreach my $parameter ( keys %parameter ) {

  KEY:
    foreach my $parameter_name (@parameter_key_to_check) {

        if ( exists $parameter{$parameter}{$parameter_name} ) {

            ## Test if element from query array exists truth hash
            check_element_exists_in_hash(
                {
                    truth_href => \%parameter,
                    queryies   => \@{ $parameter{$parameter}{$parameter_name} },
                    parameter_name => $parameter_name,
                }
            );
        }
    }
}

## Check programs in path, and executable
check_command_in_path(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter,
    }
);

## Test that the family_id and the sample_id(s) exists and are unique. Check if id sample_id contains "_".
check_unique_ids(
    {
        active_parameter_href => \%active_parameter,
        sample_ids_ref        => \@{ $active_parameter{sample_ids} },
    }
);

## Check sample_id provided in hash parameter is included in the analysis and only represented once
check_sample_id_in_parameter(
    {
        active_parameter_href => \%active_parameter,
        sample_ids_ref        => \@{ $active_parameter{sample_ids} },
        parameter_names_ref   => [qw{ expected_coverage analysis_type }],
    }
);

## Check sample_id provided in hash path parameter is included in the analysis and only represented once
check_sample_id_in_parameter_path(
    {
        active_parameter_href => \%active_parameter,
        sample_ids_ref        => \@{ $active_parameter{sample_ids} },
        parameter_names_ref   => [qw{ infile_dirs exome_target_bed }],
    }
);

## Check that VEP directory and VEP cache match
check_vep_directories(
    {
        vep_directory_path_ref  => \$active_parameter{vep_directory_path},
        vep_directory_cache_ref => \$active_parameter{vep_directory_cache},
    }
);

## Check that the supplied vcfanno toml frequency file match record 'file=' within toml config file
if (   ( $active_parameter{psv_combinevariantcallsets} > 0 )
    && ( $active_parameter{sv_vcfanno} > 0 ) )
{

    check_vcfanno_toml(
        {
            vcfanno_file_toml => $active_parameter{sv_vcfanno_config},
            vcfanno_file_freq => $active_parameter{sv_vcfanno_config_file},
        }
    );
}

check_snpsift_keys(
    {
        snpsift_annotation_files_href =>
          \%{ $active_parameter{snpsift_annotation_files} },
        snpsift_annotation_outinfo_key_href =>
          \%{ $active_parameter{snpsift_annotation_outinfo_key} },
    }
);

## Adds dynamic aggregate information from definitions to parameter hash
set_dynamic_parameter(
    {
        parameter_href => \%parameter,
        aggregates_ref => [
            ## Collects all programs that MIP can handle
            q{type:program},
            ## Collects all variant_callers
            q{program_type:variant_callers},
            ## Collects all structural variant_callers
            q{program_type:structural_variant_callers},
            ## Collect all aligners
            q{program_type:aligners},
            ## Collects all references in that are supposed to be in reference directory
            q{reference:reference_dir},
            ## Collect all programs outfiles that are redundant
            q{remove_redundant_file:yes},
        ],
    }
);

## Check correct value for program mode in MIP
check_program_mode(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter
    }
);

## Update program mode depending on dry_run_all flag
update_program_mode_with_dry_run_all(
    {
        active_parameter_href => \%active_parameter,
        programs_ref          => \@{ $parameter{dynamic_parameter}{program} },
        dry_run_all           => $active_parameter{dry_run_all},
    }
);

## Check that the correct number of aligners is used in MIP and sets the aligner flag accordingly
check_aligner(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter,
        broadcasts_ref        => \@broadcasts,
    }
);

## Check that all active variant callers have a prioritization order and that the prioritization elements match a supported variant caller.
my %priority_call_parameter = (
    variant_callers            => 'gatk_combinevariants_prioritize_caller',
    structural_variant_callers => 'sv_svdb_merge_prioritize',
);
while ( my ( $variant_caller_type, $prioritize_parameter_name ) =
    each %priority_call_parameter )
{

    ## Check if we have any active callers
    my $activate_caller_tracker = 0;
    foreach my $variant_caller (
        @{ $parameter{dynamic_parameter}{$variant_caller_type} } )
    {

        if ( $active_parameter{$variant_caller} > 0 ) {

            $activate_caller_tracker++;
        }
    }
    if ( $activate_caller_tracker > 0 ) {

        check_prioritize_variant_callers(
            {
                parameter_href        => \%parameter,
                active_parameter_href => \%active_parameter,
                variant_callers_ref =>
                  \@{ $parameter{dynamic_parameter}{$variant_caller_type} },
                parameter_names_ref => \$prioritize_parameter_name,
            }
        );
    }
}

## Broadcast set parameters info
foreach my $parameter_info (@broadcasts) {

    $log->info($parameter_info);
}

## Update program mode depending on analysis run value as some programs are not applicable for e.g. wes
my @warning_msgs = update_program_mode(
    {
        active_parameter_href => \%active_parameter,
        programs_ref => [qw{ cnvnator delly_call delly_reformat tiddit }],
        consensus_analysis_type =>
          $parameter{dynamic_parameter}{consensus_analysis_type},
    }
);

## Broadcast
if (@warning_msgs) {

    foreach my $warning_msg (@warning_msgs) {

        $log->warn($warning_msg);
    }
}
## Update prioritize flag depending on analysis run value as some programs are not applicable for e.g. wes
$active_parameter{sv_svdb_merge_prioritize} = update_prioritize_flag(
    {
        prioritize_key => $active_parameter{sv_svdb_merge_prioritize},
        programs_ref   => [qw{ cnvnator delly_call delly_reformat tiddit }],
        consensus_analysis_type =>
          $parameter{dynamic_parameter}{consensus_analysis_type},
    }
);

## Write config file for family
if ( $active_parameter{config_file_analysis} ne 0 ) {

    ## Create directory unless it already exists
    make_path( dirname( $active_parameter{config_file_analysis} ) );

    ## Remove previous analysis specific info not relevant for current run e.g. log file, sample_ids which are read from pedigree or cmd
    my @remove_keys = (qw{ associated_program });

  KEY:
    foreach my $key (@remove_keys) {

        delete $active_parameter{$key};
    }

    ## Writes a YAML hash to file
    write_yaml(
        {
            yaml_href      => \%active_parameter,
            yaml_file_path => $active_parameter{config_file_analysis},
        }
    );
    $log->info( 'Wrote: ' . $active_parameter{config_file_analysis}, "\n" );

    ## Add to qc_sample_info
    $sample_info{config_file_analysis} =
      $active_parameter{config_file_analysis};
}

## Detect the gender included in current analysis
(

    $active_parameter{found_male},
    $active_parameter{found_female},
    $active_parameter{found_other},
    $active_parameter{found_other_count},
  )
  = detect_sample_id_gender(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
  );

### Contigs
## Set contig prefix and contig names depending on reference used
set_contigs(
    {
        file_info_href         => \%file_info,
        human_genome_reference => $active_parameter{human_genome_reference},
    }
);

## Update contigs depending on settings in run (wes or if only male samples)
update_contigs_for_run(
    {
        file_info_href     => \%file_info,
        analysis_type_href => \%{ $active_parameter{analysis_type} },
        found_male         => $active_parameter{found_male},
    }
);

## Sorts array depending on reference array. NOTE: Only entries present in reference array will survive in sorted array.
@{ $file_info{sorted_select_file_contigs} } = size_sort_select_file_contigs(
    {
        file_info_href => \%file_info,
        consensus_analysis_type_ref =>
          \$parameter{dynamic_parameter}{consensus_analysis_type},
        hash_key_to_sort        => 'select_file_contigs',
        hash_key_sort_reference => 'contigs_size_ordered',
    }
);

## Write CMD to MIP log file
write_cmd_mip_log(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter,
        order_parameters_ref  => \@order_parameters,
        script_ref            => \$script,
        log_file_ref          => \$active_parameter{log_file},
        mip_version_ref       => \$VERSION,
    }
);

## Collects the ".fastq(.gz)" files from the supplied infiles directory. Checks if any of the files exist
collect_infiles(
    {
        active_parameter_href => \%active_parameter,
        indir_path_href       => \%indir_path,
        infile_href           => \%infile,
    }
);

## Reformat files for MIP output, which have not yet been created into, correct format so that a sbatch script can be generated with the correct filenames
my $uncompressed_file_switch = infiles_reformat(
    {
        active_parameter_href           => \%active_parameter,
        sample_info_href                => \%sample_info,
        file_info_href                  => \%file_info,
        infile_href                     => \%infile,
        indir_path_href                 => \%indir_path,
        infile_lane_prefix_href         => \%infile_lane_prefix,
        infile_both_strands_prefix_href => \%infile_both_strands_prefix,
        lane_href                       => \%lane,
        job_id_href                     => \%job_id,
        outaligner_dir_ref              => \$active_parameter{outaligner_dir},
        program_name                    => 'infiles_reformat',
    }
);

## Creates all fileendings as the samples is processed depending on the chain of modules activated
create_file_endings(
    {
        parameter_href          => \%parameter,
        active_parameter_href   => \%active_parameter,
        file_info_href          => \%file_info,
        infile_lane_prefix_href => \%infile_lane_prefix,
        order_parameters_ref    => \@order_parameters,
    }
);

## Create .fam file to be used in variant calling analyses
create_fam_file(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
        execution_mode        => 'system',
        fam_file_path         => catfile(
            $active_parameter{outdata_dir},
            $active_parameter{family_id},
            $active_parameter{family_id} . '.fam'
        ),
    }
);

## Add to SampleInfo
add_to_sample_info(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
        file_info_href        => \%file_info,
    }
);

############
####MAIN####
############

if ( $active_parameter{dry_run_all} == 0 ) {

    $sample_info{analysis_date} = $date_time_stamp;
    $sample_info{mip_version}   = $VERSION;
}

### Build recipes
$log->info(q{[Reference check - Reference prerequisites]});
## Check capture file prerequistes exists
foreach
  my $program_name ( @{ $parameter{exome_target_bed}{associated_program} } )
{

    if ( $active_parameter{$program_name} > 0 ) {

        ## Remove initial "p" from program_name
        substr( $program_name, 0, 1 ) = $EMPTY_STR;

        check_capture_file_prerequisites(
            {
                parameter_href              => \%parameter,
                active_parameter_href       => \%active_parameter,
                sample_info_href            => \%sample_info,
                infile_lane_prefix_href     => \%infile_lane_prefix,
                job_id_href                 => \%job_id,
                infile_list_suffix          => $file_info{exome_target_bed}[0],
                padded_infile_list_suffix   => $file_info{exome_target_bed}[1],
                padded_interval_list_suffix => $file_info{exome_target_bed}[2],
                program_name                => $program_name,
                log                         => $log,
            }
        );
    }
}

## Check human genome prerequistes exists
PROGRAM:
foreach my $program_name (
    @{ $parameter{human_genome_reference}{associated_program} } )
{

    next PROGRAM if ( $program_name eq q{mip} );

    if ( $active_parameter{$program_name} > 0 ) {

        ## Remove initial "p" from program_name
        substr( $program_name, 0, 1 ) = $EMPTY_STR;

        my $is_finished = check_human_genome_prerequisites(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                program_name            => $program_name,
                log                     => $log,
            }
        );
        last PROGRAM if ($is_finished);
    }
}
$log->info( $TAB . q{Reference check: Reference prerequisites checked} );

## Check if vt has processed references, if not try to reprocesses them before launcing modules
$log->info(q{[Reference check - Reference processed by VT]});
if (   $active_parameter{vt_decompose}
    || $active_parameter{vt_normalize} )
{

    my @to_process_references = check_references_for_vt(
        {
            parameter_href        => \%parameter,
            active_parameter_href => \%active_parameter,
            vt_references_ref =>
              \@{ $active_parameter{decompose_normalize_references} },
            log => $log,
        }
    );

  REFERENCE:
    foreach my $reference_file_path (@to_process_references) {

        $log->info(q{[VT - Normalize and decompose]});
        $log->info( $TAB . q{File: } . $reference_file_path );

        ## Split multi allelic records into single records and normalize
        analysis_vt_core(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                infile_path             => $reference_file_path,
                program_directory       => q{vt},
                decompose               => 1,
                normalize               => 1,
            }
        );
    }
}

## Split of fastq files in batches
if ( $active_parameter{psplit_fastq_file} > 0 ) {

    $log->info(q{[Split fastq files in batches]});

    my $program_name = lc q{split_fastq_file};

    foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

        ## Split input fastq files into batches of reads, versions and compress. Moves original file to subdirectory
        analysis_split_fastq_file(
            {
                parameter_href        => \%parameter,
                active_parameter_href => \%active_parameter,
                infile_href           => \%infile,
                job_id_href           => \%job_id,
                insample_directory    => $indir_path{$sample_id},
                outsample_directory   => $indir_path{$sample_id},
                sample_id             => $sample_id,
                program_name          => $program_name,
                sequence_read_batch =>
                  $active_parameter{split_fastq_file_read_batch},
            }
        );
    }

    # End here if this module is turned on
    exit;
}

### WTS
if ( $parameter{dynamic_parameter}{consensus_analysis_type} eq q{wts} ) {

    $log->info( q{Pipeline analysis type: }
          . $parameter{dynamic_parameter}{consensus_analysis_type} );

    ## Pipeline recipe for wts data
    pipeline_wts(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            indir_path_href         => \%indir_path,
            infile_href             => \%infile,
            infile_lane_prefix_href => \%infile_lane_prefix,
            lane_href               => \%lane,
            job_id_href             => \%job_id,
            outaligner_dir          => $active_parameter{outaligner_dir},
            log                     => $log,
        }
    );

    exit;
}

### WES|WGS

## GZip of fastq files
if (   ( $active_parameter{pgzip_fastq} > 0 )
    && ( $uncompressed_file_switch eq q{uncompressed} ) )
{

    $log->info(q{[Gzip for fastq files]});

    my $program_name = lc q{gzip_fastq};

  SAMPLES:
    foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

        ## Determine which sample id had the uncompressed files
      INFILES:
        foreach my $infile ( @{ $infile{$sample_id} } ) {

            my $infile_suffix = $parameter{pgzip_fastq}{infile_suffix};

            if ( $infile =~ /$infile_suffix$/ ) {

                ## Automatically gzips fastq files
                analysis_gzip_fastq(
                    {
                        parameter_href          => \%parameter,
                        active_parameter_href   => \%active_parameter,
                        sample_info_href        => \%sample_info,
                        infile_href             => \%infile,
                        infile_lane_prefix_href => \%infile_lane_prefix,
                        job_id_href             => \%job_id,
                        insample_directory      => $indir_path{$sample_id},
                        sample_id               => $sample_id,
                        program_name            => $program_name,
                    }
                );

                # Call once per sample_id
                last INFILES;
            }
        }
    }
}

# Run FastQC
if ( $active_parameter{pfastqc} > 0 ) {

    $log->info(q{[Fastqc]});

    foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

        my $program_name = lc q{fastqc};
        my $outsample_directory =
          catdir( $active_parameter{outdata_dir}, $sample_id, $program_name );
        analysis_fastqc(
            {
                parameter_href        => \%parameter,
                active_parameter_href => \%active_parameter,
                sample_info_href      => \%sample_info,
                infiles_ref           => \@{ $infile{$sample_id} },

                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                insample_directory      => $indir_path{$sample_id},
                outsample_directory     => $outsample_directory,
                sample_id               => $sample_id,
                program_name            => $program_name,
            }
        );
    }
}

# Run BWA Mem
if ( $active_parameter{pbwa_mem} > 0 ) {

    $log->info(q{[BWA Mem]});

    my $program_name = lc q{bwa_mem};

    if ( $active_parameter{dry_run_all} != 1 ) {

        if (   ( $parameter{human_genome_reference}{build_file} eq 1 )
            || ( $parameter{bwa_build_reference}{build_file} eq 1 ) )
        {

            build_bwa_prerequisites(
                {
                    parameter_href          => \%parameter,
                    active_parameter_href   => \%active_parameter,
                    sample_info_href        => \%sample_info,
                    file_info_href          => \%file_info,
                    infile_lane_prefix_href => \%infile_lane_prefix,
                    job_id_href             => \%job_id,
                    bwa_build_reference_file_endings_ref =>
                      \@{ $file_info{bwa_build_reference_file_endings} },
                    program_name => $program_name,
                }
            );
        }
    }

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

        my $outsample_directory = catdir( $active_parameter{outdata_dir},
            $sample_id, $active_parameter{outaligner_dir} );

        analysis_bwa_mem(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infiles_ref             => \@{ $infile{$sample_id} },
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                insample_directory      => $indir_path{$sample_id},
                outsample_directory     => $outsample_directory,
                sample_id               => $sample_id,
                program_name            => $program_name,
            }
        );
    }
}

## Run consecutive models
if ( $active_parameter{reduce_io} ) {

    ## Enable as program
    $active_parameter{pbamcalibrationblock} = 1;
    $log->info(q{[Bamcalibrationblock]});

    my $program_name = lc q{bamcalibrationblock};

    analysis_bamcalibrationblock(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            lane_href               => \%lane,
            job_id_href             => \%job_id,
            outaligner_dir          => $active_parameter{outaligner_dir},
            program_name            => q{bamcalibrationblock},
            log                     => $log,
        }
    );

}
else {

    ## Always run even for single samples to rename them correctly for standardised downstream processing.
    ## Will also split alignment per contig and copy to temporary directory for '-rio 1' block to enable selective removal of block submodules.
    if ( $active_parameter{ppicardtools_mergesamfiles} > 0 ) {

        $log->info(q{[Picardtools mergesamfiles]});

        my $program_name = lc q{picardtools_mergesamfiles};

        foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

            my $insample_directory = catdir( $active_parameter{outdata_dir},
                $sample_id, $active_parameter{outaligner_dir} );
            my $outsample_directory = catdir( $active_parameter{outdata_dir},
                $sample_id, $active_parameter{outaligner_dir} );

            analysis_picardtools_mergesamfiles(
                {
                    parameter_href          => \%parameter,
                    active_parameter_href   => \%active_parameter,
                    sample_info_href        => \%sample_info,
                    file_info_href          => \%file_info,
                    infile_lane_prefix_href => \%infile_lane_prefix,
                    lane_href               => \%lane,
                    job_id_href             => \%job_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    sample_id               => $sample_id,
                    program_name            => $program_name,
                }
            );
        }
    }

    # Markduplicates
    if ( $active_parameter{pmarkduplicates} > 0 ) {

        $log->info(q{[Markduplicates]});

        my $program_name = lc q{markduplicates};

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

            ## Assign directories
            my $insample_directory = catdir( $active_parameter{outdata_dir},
                $sample_id, $active_parameter{outaligner_dir} );
            my $outsample_directory = catdir( $active_parameter{outdata_dir},
                $sample_id, $active_parameter{outaligner_dir} );

            analysis_markduplicates(
                {
                    parameter_href          => \%parameter,
                    active_parameter_href   => \%active_parameter,
                    sample_info_href        => \%sample_info,
                    file_info_href          => \%file_info,
                    infile_lane_prefix_href => \%infile_lane_prefix,
                    job_id_href             => \%job_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    sample_id               => $sample_id,
                    program_name            => $program_name,
                }
            );
        }
    }

    #Run GATK realignertargetcreator/indelrealigner
    if ( $active_parameter{pgatk_realigner} > 0 ) {

        $log->info(q{[GATK realignertargetcreator/indelrealigner]});

        my $program_name = lc q{gatk_realigner};

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

            ## Assign directories
            my $insample_directory = catdir( $active_parameter{outdata_dir},
                $sample_id, $active_parameter{outaligner_dir} );
            my $outsample_directory = catdir( $active_parameter{outdata_dir},
                $sample_id, $active_parameter{outaligner_dir} );

            analysis_gatk_realigner(
                {
                    parameter_href          => \%parameter,
                    active_parameter_href   => \%active_parameter,
                    sample_info_href        => \%sample_info,
                    file_info_href          => \%file_info,
                    infile_lane_prefix_href => \%infile_lane_prefix,
                    job_id_href             => \%job_id,
                    sample_id               => $sample_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    program_name            => $program_name,
                }
            );
        }
    }

    ## Run GATK baserecalibrator/printreads
    if ( $active_parameter{pgatk_baserecalibration} > 0 ) {

        $log->info(q{[GATK baserecalibrator/printreads]});

        my $program_name = lc q{gatk_baserecalibration};

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

            ## Assign directories
            my $insample_directory = catdir( $active_parameter{outdata_dir},
                $sample_id, $active_parameter{outaligner_dir} );
            my $outsample_directory = catdir( $active_parameter{outdata_dir},
                $sample_id, $active_parameter{outaligner_dir} );

            analysis_gatk_baserecalibration(
                {
                    parameter_href          => \%parameter,
                    active_parameter_href   => \%active_parameter,
                    sample_info_href        => \%sample_info,
                    file_info_href          => \%file_info,
                    infile_lane_prefix_href => \%infile_lane_prefix,
                    job_id_href             => \%job_id,
                    sample_id               => $sample_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    program_name            => $program_name,
                }
            );
        }
    }
}

if ( $active_parameter{pchanjo_sexcheck} > 0 ) {

    $log->info(q{[Chanjo sexcheck]});

    my $program_name = lc q{chanjo_sexcheck};

  SAMPLE_IDS:
    foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

        my $insample_directory = catdir( $active_parameter{outdata_dir},
            $sample_id, $active_parameter{outaligner_dir} );
        my $outsample_directory = catdir(
            $active_parameter{outdata_dir},    $sample_id,
            $active_parameter{outaligner_dir}, q{coveragereport}
        );

        analysis_chanjo_sex_check(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                sample_id               => $sample_id,
                insample_directory      => $insample_directory,
                outsample_directory     => $outsample_directory,
                program_name            => $program_name,
            }
        );
    }
}

if ( $active_parameter{psambamba_depth} > 0 ) {

    $log->info(q{[Sambamba depth]});

    my $program_name = lc q{sambamba_depth};

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

        my $insample_directory = catdir( $active_parameter{outdata_dir},
            $sample_id, $active_parameter{outaligner_dir} );
        my $outsample_directory = catdir(
            $active_parameter{outdata_dir},    $sample_id,
            $active_parameter{outaligner_dir}, q{coveragereport}
        );

        analysis_sambamba_depth(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                sample_id               => $sample_id,
                insample_directory      => $insample_directory,
                outsample_directory     => $outsample_directory,
                program_name            => $program_name,
            }
        );
    }
}

# Run bedtools genomecov
if ( $active_parameter{pbedtools_genomecov} > 0 ) {

    $log->info(q{[Bedtools genomecov]});

    my $program_name = lc q{bedtools_genomecov};

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

        ## Assign directories
        my $insample_directory = catdir( $active_parameter{outdata_dir},
            $sample_id, $active_parameter{outaligner_dir} );
        my $outsample_directory = catdir(
            $active_parameter{outdata_dir},    $sample_id,
            $active_parameter{outaligner_dir}, q{coveragereport}
        );

        analysis_bedtools_genomecov(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                sample_id               => $sample_id,
                insample_directory      => $insample_directory,
                outsample_directory     => $outsample_directory,
                program_name            => $program_name,
            }
        );
    }
}

## Run picardtools_collectmultiplemetrics
if ( $active_parameter{ppicardtools_collectmultiplemetrics} > 0 ) {

    $log->info(q{[Picardtools collectmultiplemetrics]});

    my $program_name = lc q{picardtools_collectmultiplemetrics};

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

        ## Assign directories
        my $insample_directory = catdir( $active_parameter{outdata_dir},
            $sample_id, $active_parameter{outaligner_dir} );
        my $outsample_directory = catdir(
            $active_parameter{outdata_dir},    $sample_id,
            $active_parameter{outaligner_dir}, q{coveragereport}
        );

        analysis_picardtools_collectmultiplemetrics(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                sample_id               => $sample_id,
                insample_directory      => $insample_directory,
                outsample_directory     => $outsample_directory,
                program_name            => $program_name,
            }
        );
    }
}

## Run Picardtools_collecthsmetrics
if ( $active_parameter{ppicardtools_collecthsmetrics} > 0 ) {

    $log->info(q{[Picardtools collecthsmetrics]});

    my $program_name = lc q{picardtools_collecthsmetrics};

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

        ## Assign directories
        my $insample_directory = catdir( $active_parameter{outdata_dir},
            $sample_id, $active_parameter{outaligner_dir} );
        my $outsample_directory = catdir(
            $active_parameter{outdata_dir},    $sample_id,
            $active_parameter{outaligner_dir}, q{coveragereport}
        );

        analysis_picardtools_collecthsmetrics(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                sample_id               => $sample_id,
                insample_directory      => $insample_directory,
                outsample_directory     => $outsample_directory,
                program_name            => $program_name,
            }
        );
    }
}

## Run Rcovplot scripts
if ( $active_parameter{prcovplots} > 0 ) {

    if ( $active_parameter{pbedtools_genomecov} > 0 ) {

        $log->info(q{[Rcovplots]});

        my $program_name = lc q{rcovplots};

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

            ## Assign directories
            my $insample_directory = catdir( $active_parameter{outdata_dir},
                $sample_id, $active_parameter{outaligner_dir} );
            my $outsample_directory = catdir(
                $active_parameter{outdata_dir},    $sample_id,
                $active_parameter{outaligner_dir}, q{coveragereport}
            );
            analysis_rcoverageplots(
                {
                    parameter_href          => \%parameter,
                    active_parameter_href   => \%active_parameter,
                    sample_info_href        => \%sample_info,
                    file_info_href          => \%file_info,
                    lane_href               => \%lane,
                    infile_lane_prefix_href => \%infile_lane_prefix,
                    job_id_href             => \%job_id,
                    sample_id               => $sample_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    program_name            => $program_name,
                }
            );
        }
    }
}

#Run CNVnator
if ( $active_parameter{pcnvnator} > 0 ) {

    $log->info(q{[CNVnator]});

    my $program_name = lc q{cnvnator};

    foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

        ## Assign directories
        my $insample_directory = catdir( $active_parameter{outdata_dir},
            $sample_id, $active_parameter{outaligner_dir} );

        my $outsample_directory = catdir(
            $active_parameter{outdata_dir},    $sample_id,
            $active_parameter{outaligner_dir}, $program_name
        );

        analysis_cnvnator(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                sample_id               => $sample_id,
                insample_directory      => $insample_directory,
                outsample_directory     => $outsample_directory,
                program_name            => $program_name,
            }
        );
    }
}

if ( $active_parameter{pdelly_call} > 0 ) {    #Run delly_call

    $log->info("[Delly_call]\n");

    foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

        delly_call(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                sample_id_ref           => \$sample_id,
                program_name            => "delly_call",
            }
        );
    }
}

if ( $active_parameter{pdelly_reformat} > 0 )
{    #Run Delly merge, regenotype, bcftools merge

    $log->info("[Delly_reformat]\n");

    delly_reformat(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => "delly_reformat",
        }
    );
}

if ( $active_parameter{pmanta} > 0 ) {    #Run Manta

    $log->info(q{[Manta]});
    my $program_name = lc q{manta};

    my $outfamily_directory = catfile(
        $active_parameter{outdata_dir},    $active_parameter{family_id},
        $active_parameter{outaligner_dir}, $program_name,
    );

    analysis_manta(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => $program_name,
            outfamily_directory     => $outfamily_directory,
        }
    );
}

if ( $active_parameter{ptiddit} > 0 ) {    #Run Tiddit

    $log->info(q{[Tiddit]});
    my $program_name = lc q{tiddit};

    my $outfamily_directory = catfile(
        $active_parameter{outdata_dir},    $active_parameter{family_id},
        $active_parameter{outaligner_dir}, $program_name,
    );

    analysis_tiddit(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => $program_name,
            outfamily_directory     => $outfamily_directory,
        }
    );
}

if ( $active_parameter{psv_combinevariantcallsets} > 0 )
{   #Run combinevariantcallsets. For all Sample_ids and StructuralVariantCallers

    $log->info("[SV combinevariantcallsets]\n");

    sv_combinevariantcallsets(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => "sv_combinevariantcallsets",
        }
    );
}

# Run sv_varianteffectpredictor. Family-level
if ( $active_parameter{psv_varianteffectpredictor} > 0 ) {

    $log->info(q{[SV varianteffectpredictor]});

    my $program_name = lc q{sv_varianteffectpredictor};

    analysis_vep_sv(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            contigs_ref             => \@{ $file_info{contigs_sv} },
            program_name            => $program_name,
        }
    );
}

if ( $active_parameter{psv_vcfparser} > 0 ) { #Run sv_vcfparser. Done per family

    $log->info("[SV vcfparser]\n");

    sv_vcfparser(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => "sv_vcfparser",
        }
    );
}

if ( $active_parameter{psv_rankvariant} > 0 )
{    #Run sv_rankvariant. Done per family

    $log->info("[SV rankvariant]\n");

    sv_rankvariant(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => "sv_rankvariant",
        }
    );
}

if ( $active_parameter{psv_reformat} > 0 ) {   #Run sv_reformat. Done per family

    $log->info("[SV reformat]\n");

    sv_reformat(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => "sv_reformat",
        }
    );
}

if ( $active_parameter{psamtools_mpileup} > 0 ) {    #Run samtools mpileup

    $log->info("[Samtools mpileup]\n");

    msamtools_mpileup(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => "samtools_mpileup",
        }
    );
}

if ( $active_parameter{pfreebayes} > 0 ) {    #Run Freebayes

    $log->info("[Freebayes]\n");

    freebayes(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => "freebayes",
        }
    );
}

if ( $active_parameter{pgatk_haplotypecaller} > 0 ) {  #Run GATK haplotypecaller

    $log->info(q{[GATK haplotypecaller]});
    my $program_name = lc q{gatk_haplotypecaller};

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

        my $insample_directory = catdir( $active_parameter{outdata_dir},
            $sample_id, $active_parameter{outaligner_dir} );
        my $outsample_directory = catdir(
            $active_parameter{outdata_dir},    $sample_id,
            $active_parameter{outaligner_dir}, q{gatk}
        );

        analysis_gatk_haplotypecaller(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                sample_id               => $sample_id,
                insample_directory      => $insample_directory,
                outsample_directory     => $outsample_directory,
                program_name            => $program_name,
            }
        );
    }
}

if ( $active_parameter{pgatk_genotypegvcfs} > 0 )
{    #Run GATK genotypegvcfs. Done per family

    $log->info("[GATK genotypegvcfs]\n");

    mgatk_genotypegvcfs(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => "gatk_genotypegvcfs",
        }
    );

    gatk_concatenate_genotypegvcfs(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => "gatk_genotypegvcfs",
        }
    );
}

if ( $active_parameter{pgatk_variantrecalibration} > 0 )
{    #Run GATK VariantRecalibrator/ApplyRecalibration. Done per family

    $log->info("[GATK variantrecalibrator/applyrecalibration]\n");

    gatk_variantrecalibration(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => "gatk_variantrecalibration",
        }
    );
}

if ( $active_parameter{pgatk_combinevariantcallsets} > 0 )
{    #Run gatk_combinevariantcallsets. Done per family

    $log->info("[GATK combinevariantcallsets]\n");

    gatk_combinevariantcallsets(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => "gatk_combinevariantcallsets",
        }
    );
}

if ( $active_parameter{ppeddy} > 0 ) {    #Run plink. Done per family

    $log->info("[Peddy]\n");

    mpeddy(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => "peddy",
        }
    );
}

if ( $active_parameter{pplink} > 0 ) {    #Run plink. Done per family

    $log->info("[Plink]\n");

    mplink(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => "plink",
        }
    );
}

if ( $active_parameter{pvariant_integrity} > 0 ) {

    #Run variant_integrity. Done per family
    $log->info(q{[Variant_integrity]});

    my $program_name = lc q{variant_integrity};

    my $infamily_directory = catdir(
        $active_parameter{outdata_dir},
        $active_parameter{family_id},
        $active_parameter{outaligner_dir}
    );
    my $outfamily_directory = catfile(
        $active_parameter{outdata_dir},
        $active_parameter{family_id},
        $active_parameter{outaligner_dir},
        q{casecheck}, $program_name
    );

    analysis_variant_integrity(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => $program_name,
            infamily_directory      => $infamily_directory,
            outfamily_directory     => $outfamily_directory,
        }
    );
}

if ( $active_parameter{pevaluation} > 0 ) {    #Run evaluation. Done per family

    foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

        if ( $sample_id =~ /$active_parameter{nist_id}/ ) {

            $log->info("[Evaluation]\n");

            my $program_name = lc q{evaluation};

            ## Assign directories
            my $infamily_directory = catdir(
                $active_parameter{outdata_dir},
                $active_parameter{family_id},
                $active_parameter{outaligner_dir}
            );
            my $outfamily_directory = catfile(
                $active_parameter{outdata_dir},
                $active_parameter{family_id},
                $active_parameter{outaligner_dir},
                $program_name
            );
            analysis_picardtools_genotypeconcordance(
                {
                    parameter_href          => \%parameter,
                    active_parameter_href   => \%active_parameter,
                    sample_info_href        => \%sample_info,
                    file_info_href          => \%file_info,
                    infile_lane_prefix_href => \%infile_lane_prefix,
                    job_id_href             => \%job_id,
                    sample_id               => $sample_id,
                    call_type               => q{BOTH},
                    infamily_directory      => $infamily_directory,
                    outfamily_directory     => $outfamily_directory,
                    program_name            => $program_name,
                }
            );
        }

    }
}

if ( $active_parameter{pgatk_variantevalall} > 0 )
{    #Run GATK varianteval for all variants. Done per sample_id

    $log->info("[GATK variantevalall]\n");

    foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

        gatk_variantevalall(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                sample_id_ref           => \$sample_id,
                program_name            => "gatk_variantevalall",
            }
        );
    }
}

### If no males or other remove contig Y from all downstream analysis
my @file_info_contig_keys = ( "contigs_size_ordered", "contigs" );

foreach my $key (@file_info_contig_keys) {

    ## Removes contig_names from contigs array if no male or 'other' found
    @{ $file_info{$key} } = delete_male_contig(
        {
            contigs_ref => \@{ $file_info{$key} },
            found_male  => $active_parameter{found_male},
        }
    );
}

if ( $active_parameter{reduce_io} ) {    #Run consecutive models

    $active_parameter{pvariantannotationblock} = 1;    #Enable as program
    $log->info("[Variantannotationblock]\n");

    variantannotationblock(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            outaligner_dir_ref      => \$active_parameter{outaligner_dir},
            call_type               => q{BOTH},
            program_name            => "variantannotationblock",
        }
    );
}
else {

    if ( $active_parameter{pprepareforvariantannotationblock} > 0 ) {

        $log->info("[Prepareforvariantannotationblock]\n");

        prepareforvariantannotationblock(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                call_type               => q{BOTH},
                program_name            => "prepareforvariantannotationblock",
            }
        );
    }
    if ( $active_parameter{prhocall} > 0 ) {    #Run rhocall. Done per family

        $log->info("[Rhocall]\n");

        rhocall(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                call_type               => q{BOTH},
                program_name            => "rhocall",
            }
        );
    }
    if ( $active_parameter{pvt} > 0 ) {    #Run vt. Done per family

        $log->info("[Vt]\n");

        vt(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                call_type               => q{BOTH},
                program_name            => "vt",
            }
        );
    }
    ## Run varianteffectpredictor {family-level}
    if ( $active_parameter{pvarianteffectpredictor} > 0 ) {

        $log->info(q{[Varianteffectpredictor]});

        my $program_name = lc q{varianteffectpredictor};

        analysis_vep(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                call_type               => q{BOTH},
                program_name            => $program_name,
            }
        );
    }
    if ( $active_parameter{pvcfparser} > 0 ) {  #Run pvcfparser. Done per family

        $log->info("[Vcfparser]\n");

        mvcfparser(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                call_type               => q{BOTH},
                program_name            => "vcfparser",
            }
        );
    }

    if ( $active_parameter{psnpeff} > 0 ) {    #Run snpEff. Done per family

        $log->info("[Snpeff]\n");

        snpeff(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                call_type               => q{BOTH},
                program_name            => "snpeff",
            }
        );
    }
    if ( $active_parameter{prankvariant} > 0 )
    {    #Run rankvariant. Done per family

        $log->info("[Rankvariant]\n");

        rankvariant(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                call_type               => q{BOTH},
                program_name            => "rankvariant",
            }
        );
    }
    if ( $active_parameter{pendvariantannotationblock} > 0 )
    {    #Run endvariantannotationblock. Done per family

        $log->info("[Endvariantannotationblock]\n");

        endvariantannotationblock(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                call_type               => q{BOTH},
                program_name            => "endvariantannotationblock",
            }
        );
    }
}

if ( $active_parameter{pgatk_variantevalexome} > 0 )
{    #Run GATK varianteval for exome variants. Done per sample_id

    $log->info("[GATK variantevalexome]\n");

    foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

        gatk_variantevalexome(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                infile_lane_prefix_href => \%infile_lane_prefix,
                job_id_href             => \%job_id,
                sample_id_ref           => \$sample_id,
                program_name            => "gatk_variantevalexome",
            }
        );
    }
}

if ( $active_parameter{pqccollect} > 0 ) {    #Run qccollect. Done per family

    $log->info("[Qccollect]\n");

    mqccollect(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => "qccollect",
        }
    );
}

if ( $active_parameter{pmultiqc} > 0 ) {

    $log->info("[Multiqc]\n");

    mmultiqc(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => "multiqc",
        }
    );
}

if ( $active_parameter{premoveredundantfiles} > 0 )
{    #Sbatch generation of removal of alignment files

    $log->info("[Removal of redundant files]\n");

    removeredundantfiles(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            lane_href               => \%lane,
            program_name            => "removeredundantfiles",
        }
    );
}

if (   ( $active_parameter{panalysisrunstatus} == 1 )
    && ( !$active_parameter{dry_run_all} ) )
{

    $sample_info{analysisrunstatus} =
      "not_finished";    #Add analysis run status flag.
}

if ( $active_parameter{panalysisrunstatus} > 0 ) {

    $log->info("[Analysis run status]\n");

    analysisrunstatus(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => "analysisrunstatus",
        }
    );
}

if (   ( $active_parameter{psacct} > 0 )
    && ( $active_parameter{dry_run_all} == 0 ) )
{    #Sbatch generation of sacct job_id data info

    $log->info("[Sacct]\n");

    msacct(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => "sacct",
        }
    );
}

#Write QC for programs used in analysis
if ( $active_parameter{sample_info_file} ne 0 ) { #Write SampleInfo to yaml file

    ## Writes a YAML hash to file
    write_yaml(
        {
            yaml_href      => \%sample_info,
            yaml_file_path => $active_parameter{sample_info_file},
        }
    );
    $log->info( "Wrote: " . $active_parameter{sample_info_file}, "\n" );
}

######################
####Sub routines######
######################

sub build_usage {

##Function : Build the USAGE instructions
##Returns  :
##Arguments: $program_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $program_name;

    my $tmpl = {
        program_name => {
            default     => basename($PROGRAM_NAME),
            strict_type => 1,
            store       => \$program_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options] -ifd [infile_dirs=sample_id] -rd [reference_dir] -p [project_id] -s [sample_ids,.,.,.,n] -em [email] -osd [outscript_dir] -odd [outdata_dir] -f [family_id] -p[program] -at [sample_id=analysis_type]

    ####MIP
    -ifd/--infile_dirs Infile directory(s) (Hash infile_dirs=sample_id; mandatory)
    -rd/--reference_dir Reference(s) directory (mandatory)
    -p/--project_id The project ID (mandatory)
    -s/--sample_ids The sample ID(s)(comma sep; mandatory)
    -odd/--outdata_dir The data output directory (mandatory)
    -osd/--outscript_dir The script files (.sh) output directory (mandatory)
    -f/--family_id Group id of samples to be compared (defaults to "", (Ex: 1 for IDN 1-1-1A))
    -sck/--supported_capture_kit Set the capture kit acronym shortcut in pedigree file
    -dnr/--decompose_normalize_references Set the references to be decomposed and normalized (defaults: "gatk_realigner_indel_known_sites", "gatk_baserecalibration_known_sites","gatk_haplotypecaller_snp_known_set", "gatk_variantrecalibration_resource_snv", "gatk_variantrecalibration_resource_indel", "vt_genmod_filter_1000g", "sv_vcfanno_config_file", "gatk_varianteval_gold", "gatk_varianteval_dbsnp","snpsift_annotation_files")
    -ped/--pedigree_file Meta data on samples (defaults to "")
    -hgr/--human_genome_reference Fasta file for the human genome reference (defaults to "GRCh37_homo_sapiens_-d5-.fasta;1000G decoy version 5")
    -ald/--outaligner_dir Setting which aligner out directory was used for alignment in previous analysis (defaults to "{outdata_dir}{outaligner_dir}")
    -at/--analysis_type Type of analysis to perform (sample_id=analysis_type, defaults to "wgs";Valid entries: "wgs", "wes", "wts")
    -pl/--platform Platform/technology used to produce the reads (defaults to "ILLUMINA")
    -ec/--expected_coverage Expected mean target coverage for analysis (sample_id=expected_coverage, defaults to "")
    -c/--config_file YAML config file for analysis parameters (defaults to "")
    -ccp/--cluster_constant_path Set the cluster constant path (defaults to "")
    -acp/--analysis_constant_path Set the analysis constant path (defaults to "analysis")
    -cfa/--config_file_analysis Write YAML configuration file for analysis parameters (defaults to "")
    -sif/--sample_info_file YAML file for sample info used in the analysis (defaults to "{outdata_dir}/{family_id}/{family_id}_qc_sample_info.yaml")
    -dra/--dry_run_all Sets all programs to dry run mode i.e. no sbatch submission (defaults to "0" (=no))
    -jul/--java_use_large_pages Use large page memory. (defaults to "0" (=no))
    -ges/--genomic_set Selection of relevant regions post alignment (Format=sorted BED; defaults to "")
    -rio/--reduce_io Run consecutive models at nodes (defaults to "0" (=no))
    -riu/--replace_iupac Replace IUPAC code in alternative alleles with N (defaults to "0" (=no))
    -pp/--print_program Print all programs that are supported
    -ppm/--print_program_mode Print all programs that are supported in: 0 (off mode), 1 (on mode), 2 (dry run mode; defaults to "2")
    -l/--log_file Mip log file (defaults to "{outdata_dir}/{family_id}/mip_log/{date}/{scriptname}_{timestamp}.log")
    -h/--help Display this help message
    -v/--version Display version of MIP

    ####Bash
    -bse/--bash_set_errexit Set errexit in bash scripts (defaults to "0")
    -bsu/--bash_set_nounset Set nounset in bash scripts (defaults to "0")
    -bsp/--bash_set_pipefail Set pipefail in bash scripts (defaults to "0")
    -mot/--module_time Set the time allocation for each module (Format: module "program name"=time(Hours))
    -mcn/--module_core_number Set the number of cores for each module (Format: module "program_name"=X(cores))
    -mcn/--max_cores_per_node The maximum number of processor cores per node used in the analysis (defaults to "16")
    -nrm/--node_ram_memory The RAM memory size of the node(s) in GigaBytes (Defaults to 24)
    -tmd/--temp_directory Set the temporary directory for all programs (defaults to "/scratch/SLURM_JOB_ID";supply whole path)
    -em/--email E-mail (defaults to "")
    -emt/--email_types E-mail type (defaults to FAIL (=FAIL);Options: BEGIN (=BEGIN) and/or F (=FAIL) and/or END=(END))
    -qos/--slurm_quality_of_service SLURM quality of service command in sbatch scripts (defaults to "normal")
    -sen/--source_environment_commands Source environment command in sbatch scripts (defaults to "")

    ####Programs
    -psfq/--psplit_fastq_file Split fastq files in batches of X reads and exits (defaults to "0" (=no))
      -sfqrdb/--split_fastq_file_read_batch The number of sequence reads to place in each batch (defaults to "25,000,000")
    -pgz/--pgzip_fastq Gzip fastq files (defaults to "0" (=no))
    -pfqc/--pfastqc Sequence quality analysis using FastQC (defaults to "0" (=no))

    ##BWA
    -pmem/--pbwa_mem Align reads using Bwa Mem (defaults to "0" (=no))
      -memhla/--bwa_mem_hla Apply HLA typing (defaults to "1" (=yes))
      -memrdb/--bwa_mem_rapid_db Selection of relevant regions post alignment (defaults to "")
      -memcrm/--bwa_mem_cram Use CRAM-format for additional output file (defaults to "1" (=yes))
      -memsts/--bwa_mem_bamstats Collect statistics from BAM files (defaults to "1" (=yes))
      -memssm/--bwa_sambamba_sort_memory_limit Set the memory limit for Sambamba sort after bwa alignment (defaults to "32G")

    ##Picardtools
    -ptp/--picardtools_path Path to Picardtools. Mandatory for use of Picardtools (defaults to "")
    -pptm/--ppicardtools_mergesamfiles Merge (BAM file(s) ) using Picardtools mergesamfiles or rename single samples for downstream processing (defaults to "0" (=no))

    ##Markduplicates
    -pmd/--pmarkduplicates Markduplicates using either Picardtools markduplicates or sambamba markdup (defaults to "0" (=no))
    -mdpmd/--markduplicates_picardtools_markduplicates Markduplicates using Picardtools markduplicates (defaults to "1" (=yes))
    -mdsmd/--markduplicates_sambamba_markdup Markduplicates using Sambamba markduplicates (defaults to "0" (=no))
      -mdshts/--markduplicates_sambamba_markdup_hash_table_size Sambamba size of hash table for finding read pairs (defaults to "262144")
      -mdsols/--markduplicates_sambamba_markdup_overflow_list_size Sambamba size of the overflow list (defaults to "200000")
      -mdsibs/--markduplicates_sambamba_markdup_io_buffer_size Sambamba size of the io buffer for reading and writing BAM during the second pass (defaults to "2048")

    ###Coverage calculations
    -pchs/--pchanjo_sexcheck Predicts gender from sex chromosome coverage (defaults to "0" (=no))
      -chslle/--chanjo_sexcheck_log_level Set chanjo sex log level (defaults to "DEBUG")
    -psdt/--psambamba_depth Sambamba depth coverage analysis (defaults to "0" (=no))
      -sdtmod/--sambamba_depth_mode Mode unit to print the statistics on (defaults to "region")
      -sdtcut/--sambamba_depth_cutoffs Read depth cutoff (comma sep; defaults to "10", "20", "30", "50", "100")
      -sdtbed/--sambamba_depth_bed Reference database (defaults to "CCDS.current.bed")
      -sdtbaq/--sambamba_depth_base_quality Do not count bases with lower base quality (defaults to "10")
      -stdmaq/--sambamba_depth_mapping_quality  Do not count reads with lower mapping quality (defaults to "10")
      -stdndu/--sambamba_depth_noduplicates Do not include duplicates in coverage calculation (defaults to "1" (=yes))
      -stdfqc/--sambamba_depth_quality_control Do not include reads with failed quality control (defaults to "1" (=yes))
    -pbgc/--pbedtools_genomecov Genome coverage calculation using bedtools genomecov (defaults to "0" (=no))
     -bgcmc/--bedtools_genomecov_max_coverage Max coverage depth when using '-pbedtools_genomecov' (defaults to "30")
    -pptcmm/--ppicardtools_collectmultiplemetrics Metrics calculation using Picardtools CollectMultipleMetrics (defaults to "0" (=no))
    -pptchs/--ppicardtools_collecthsmetrics Capture calculation using Picardtools Collecthsmetrics (defaults to "0" (=no))
      -extb/--exome_target_bed Exome target bed file per sample_id (defaults to "latest_supported_capturekit.bed"; -extb file.bed=Sample_idX,Sample_idY -extb file.bed=Sample_idZ)
    -prcp/--prcovplots Plots of genome coverage using rcovplots (defaults to "0" (=no))

    ###Structural variant callers
    -pcnv/--pcnvnator Structural variant calling using CNVnator (defaults to "0" (=no))
      -cnvhbs/--cnv_bin_size CNVnator bin size (defaults to "1000")
      -cnvrld/--cnv_root_ld_lib Root ld library path (defaults to "")
    -pdelc/--pdelly_call Structural variant calling using Delly (defaults to "0" (=no))
    -pdel/--pdelly_reformat Merge, regenotype and filter using Delly (defaults to "0" (=no))
      -deltyp/--delly_types Type of SV to call (defaults to "DEL,DUP,INV,TRA"; comma sep)
      -delexc/--delly_exclude_file Exclude centomere and telemore regions in delly calling (defaults to "hg19_human_excl_-0.7.6-.tsv")
    -pmna/--pmanta Structural variant calling using Manta (defaults to "0" (=no))
    -ptid/--ptiddit Structural variant calling using Tiddit (defaults to "0" (=no))
      -tidmsp/--tiddit_minimum_number_supporting_pairs The minimum number of supporting reads (defaults to "6")
    -psvc/--psv_combinevariantcallsets Combine variant call sets (defaults to "0" (=no))
      -svcvtd/--sv_vt_decompose Split multi allelic records into single records (defaults to "1" (=yes))
      -svsvdbmp/--sv_svdb_merge_prioritize The prioritization order of structural variant callers.(defaults to ""; comma sep; Options: manta|delly|cnvnator|tiddit)
      -svcbtv/--sv_bcftools_view_filter Include structural variants with PASS in FILTER column (defaults to "1" (=yes))
      -svcdbq/--sv_svdb_query Annotate structural variants using svdb query (defaults to "1" (=yes))
      -svcdbqd/--sv_svdb_query_db_files Database file(s) for annotation (defaults to "")
      -svcvan/--sv_vcfanno Annotate structural variants (defaults to "1" (=yes)
      -svcval/--sv_vcfanno_lua vcfAnno lua postscripting file (defaults to "")
      -svcvac/--sv_vcfanno_config vcfAnno toml config (defaults to "")
      -svcvacf/--sv_vcfanno_config_file Annotation file within vcfAnno config toml file (defaults to "GRCh37_all_sv_-phase3_v2.2013-05-02-.vcf.gz")
      -svcvah/--sv_vcfannotation_header_lines_file Adjust for postscript by adding required header lines to vcf (defaults to "")
      -svcgmf/--sv_genmod_filter Remove common structural variants from vcf (defaults to "1" (=yes))
      -svcgfr/--sv_genmod_filter_1000g Genmod annotate structural variants from 1000G reference (defaults to "GRCh37_all_wgs_-phase3_v5b.2013-05-02-.vcf.gz")
      -svcgft/--sv_genmod_filter_threshold Threshold for filtering structural variants (defaults to "0.10")
      -svcbcf/--sv_combinevariantcallsets_bcf_file Produce a bcf from the CombineStructuralVariantCallSet vcf (defaults to "1" (=yes))
    -psvv/--psv_varianteffectpredictor Annotate SV variants using VEP (defaults to "0" (=no))
      -svvepf/--sv_vep_features VEP features (defaults to ("hgvs","symbol","numbers","sift","polyphen","humdiv","domains","protein","ccds","uniprot","biotype","regulatory", "tsl", "canonical", "per_gene", "appris"); comma sep)
      -svveppl/--sv_vep_plugins VEP plugins (defaults to ("UpDownDistance, LoFtool, LoF"); comma sep)
    -psvvcp/--psv_vcfparser Parse structural variants using vcfParser.pl (defaults to "0" (=no))
      -svvcpvt/--sv_vcfparser_vep_transcripts Parse VEP transcript specific entries (defaults to "0" (=no))
      -vcppg/--vcfparser_per_gene Keep only most severe consequence per gene (defaults to "1" (=yes))
      -svvcprff/--sv_vcfparser_range_feature_file Range annotations file (defaults to ""; tab-sep)
      -svvcprfa/--sv_vcfparser_range_feature_annotation_columns Range annotations feature columns (defaults to ""; comma sep)
      -svvcpsf/--sv_vcfparser_select_file File containging list of genes to analyse seperately (defaults to "";tab-sep file and HGNC Symbol required)
      -svvcpsfm/--sv_vcfparser_select_file_matching_column Position of HGNC Symbol column in select file (defaults to "")
      -svvcpsfa/--sv_vcfparser_select_feature_annotation_columns Feature columns to use in annotation (defaults to ""; comma sep)
    -psvr/--psv_rankvariant Ranking of annotated SV variants (defaults to "0" (=no))
      -svravanr/--sv_genmod_annotate_regions Use predefined gene annotation supplied with genmod for defining genes (defaults to "1" (=yes))
      -svravgft/--sv_genmod_models_family_type Use one of the known setups (defaults to "mip")
      -svravwg/--sv_genmod_models_whole_gene Allow compound pairs in intronic regions (defaults to "0" (=yes))
      -svravrpf/--sv_genmod_models_reduced_penetrance_file File containg genes with reduced penetrance (defaults to "")
      -svravrm/--sv_rank_model_file Rank model config file (defaults to "")
    -psvre/--psv_reformat Concatenating files (defaults to "0" (=no))
      -svrevbf/--sv_rankvariant_binary_file Produce binary file from the rank variant chromosome sorted vcfs (defaults to "1" (=yes))
      -svrergf/--sv_reformat_remove_genes_file Remove variants in hgnc_ids (defaults to "")

    ##Samtools
    -psmp/--psamtools_mpileup Variant calling using samtools mpileup and bcftools (defaults to "0" (=no))

    ##Freebayes
    -pfrb/--pfreebayes Variant calling using Freebayes and bcftools (defaults to "0" (=no))

    ##GATK
    -gtp/--gatk_path  Path to GATK. Mandatory for use of GATK (defaults to "")
    -gll/--gatk_logging_level Set the GATK log level (defaults to "INFO")
    -gbdv/--gatk_bundle_download_version  GATK FTP bundle download version.(defaults to "2.8")
    -gdco/--gatk_downsample_to_coverage Coverage to downsample to at any given locus (defaults to "1000")
    -gdai/--gatk_disable_auto_index_and_file_lock Disable auto index creation and locking when reading rods (defaults to "0" (=no))
    -pgra/--pgatk_realigner Realignments of reads using GATK ReAlignerTargetCreator/IndelRealigner (defaults to "0" (=no))
      -graks/--gatk_realigner_indel_known_sites GATK ReAlignerTargetCreator/IndelRealigner known indel site (defaults to "GRCh37_1000g_indels_-phase1-.vcf", "GRCh37_mills_and_1000g_indels_-gold_standard-.vcf")
    -pgbr/--pgatk_baserecalibration Recalibration of bases using GATK BaseReCalibrator/PrintReads (defaults to "0" (=no))
      -gbrcov/--gatk_baserecalibration_covariates GATK BaseReCalibration covariates (defaults to "ReadGroupCovariate", "ContextCovariate", "CycleCovariate", "QualityScoreCovariate")
      -gbrkst/--gatk_baserecalibration_known_sites GATK BaseReCalibration known SNV and INDEL sites (defaults to "GRCh37_dbsnp_-138-.vcf", "GRCh37_1000g_indels_-phase1-.vcf", "GRCh37_mills_and_1000g_indels_-gold_standard-.vcf")
      -gbrrf/--gatk_baserecalibration_read_filters Filter out reads according to set filter (defaults to "1" (=yes))
      -gbrdiq/--gatk_baserecalibration_disable_indel_qual Disable indel quality scores (defaults to "1" (=yes))
      -gbrsqq/--gatk_baserecalibration_static_quantized_quals Static binning of base quality scores (defaults to "10,20,30,40"; comma sep)
    -pghc/--pgatk_haplotypecaller Variant discovery using GATK HaplotypeCaller (defaults to "0" (=no))
      -ghcann/--gatk_haplotypecaller_annotation GATK HaploTypeCaller annotations (defaults to "BaseQualityRankSumTest", "ChromosomeCounts", "Coverage", "DepthPerAlleleBySample", "FisherStrand", "MappingQualityRankSumTest", "QualByDepth", "RMSMappingQuality", "ReadPosRankSumTest", "StrandOddsRatio")
      -ghckse/--gatk_haplotypecaller_snp_known_set GATK HaplotypeCaller dbSNP set for annotating ID columns (defaults to "GRCh37_dbsnp_-138-.vcf")
      -ghcscb/--gatk_haplotypecaller_no_soft_clipped_bases Do not include soft clipped bases in the variant calling (defaults to "1" (=yes))
      -ghcpim/--gatk_haplotypecaller_pcr_indel_model The PCR indel model to use (defaults to "None"; Set to "0" to disable)
    -pggt/--pgatk_genotypegvcfs Merge gVCF records using GATK GenotypeGVCFs (defaults to "0" (=no))
      -ggtgrl/--gatk_genotypegvcfs_ref_gvcf GATK GenoTypeGVCFs gVCF reference infile list for joint genotyping (defaults to "")
      -ggtals/--gatk_genotypegvcfs_all_sites Emit non-variant sites to the output vcf file (defaults to "0" (=no))
      -ggbcf/gatk_concatenate_genotypegvcfs_bcf_file Produce a bcf from the GATK ConcatenateGenoTypeGVCFs vcf (defaults to "1" (=yes))
    -pgvr/--pgatk_variantrecalibration Variant recalibration using GATK VariantRecalibrator/ApplyRecalibration (defaults to "0" (=no))
      -gvrann/--gatk_variantrecalibration_annotations Annotations to use with GATK VariantRecalibrator (defaults to "QD", "MQRankSum", "ReadPosRankSum", "FS", "SOR", "DP")
      -gvrres/gatk_variantrecalibration_resource_snv Resource to use with GATK VariantRecalibrator in SNV|BOTH mode (defaults to "GRCh37_dbsnp_-138-.vcf: dbsnp,known=true,training=false,truth=false,prior=2.0, GRCh37_hapmap_-3.3-.vcf: hapmap,VCF,known=false,training=true,truth=true,prior=15.0, GRCh37_1000g_omni_-2.5-.vcf: omni,VCF,known=false,training=true,truth=false,prior=12.0, GRCh37_1000g_snps_high_confidence_-phase1-.vcf: 1000G,known=false,training=true,truth=false,prior=10.0")
      -gvrrei/gatk_variantrecalibration_resource_indel Resource to use with GATK VariantRecalibrator in INDEL|BOTH (defaults to "GRCh37_dbsnp_-138-.vcf: dbsnp,known=true,training=false,truth=false,prior=2.0, GRCh37_mills_and_1000g_indels_-gold_standard-.vcf: mills,VCF,known=true,training=true,truth=true,prior=12.0")
      -gvrstf/--gatk_variantrecalibration_snv_tsfilter_level The truth sensitivity level for snvs at which to start filtering used in GATK VariantRecalibrator (defaults to "99.9")
      -gvritf/--gatk_variantrecalibration_indel_tsfilter_level The truth sensitivity level for indels at which to start filtering used in GATK VariantRecalibrator (defaults to "99.9")
      -gvrdpa/--gatk_variantrecalibration_dp_annotation Use the DP annotation in variant recalibration. (defaults to "1" (=yes))
      -gvrsmg/--gatk_variantrecalibration_snv_max_gaussians Use hard filtering for snvs (defaults to "0" (=no))
      -gvrimg/--gatk_variantrecalibration_indel_max_gaussians Use hard filtering for indels (defaults to "1" (=yes))
      -gvrevf/--gatk_variantrecalibration_exclude_nonvariants_file Produce a vcf containing non-variant loci alongside the vcf only containing non-variant loci after GATK VariantRecalibrator (defaults to "0" (=no))
      -gvrbcf/--gatk_variantrecalibration_bcf_file Produce a bcf from the GATK VariantRecalibrator vcf (defaults to "1" (=yes))
      -gcgpss/--gatk_calculategenotypeposteriors_support_set GATK CalculateGenotypePosteriors support set (defaults to "1000g_sites_GRCh37_phase3_v4_20130502.vcf")
    -pgcv/--pgatk_combinevariantcallsets Combine variant call sets (defaults to "0" (=no))
      -gcvbcf/--gatk_combinevariantcallsets_bcf_file Produce a bcf from the GATK CombineVariantCallSet vcf (defaults to "1" (=yes))
      -gcvgmo/--gatk_combinevariants_genotype_merge_option Type of merge to perform (defaults to "PRIORITIZE")
      -gcvpc/--gatk_combinevariants_prioritize_caller The prioritization order of variant callers.(defaults to ""; comma sep; Options: gatk|samtools|freebayes)
    -pgvea/--pgatk_variantevalall Variant evaluation using GATK varianteval for all variants  (defaults to "0" (=no))
    -pgvee/--pgatk_variantevalexome Variant evaluation using GATK varianteval for exonic variants  (defaults to "0" (=no))
      -gveedbs/--gatk_varianteval_dbsnp DbSNP file used in GATK varianteval (defaults to "dbsnp_GRCh37_138_esa_129.vcf")
      -gveedbg/--gatk_varianteval_gold Gold indel file used in GATK varianteval (defaults to "GRCh37_mills_and_1000g_indels_-gold_standard-.vcf")

    ###Annotation
    -ppvab/--pprepareforvariantannotationblock Prepare for variant annotation block by copying and splitting files per contig (defaults to "0" (=no))
    -prhc/--prhocall Rhocall performs annotation of variants in autozygosity regions (defaults to "0" (=no))
      -rhcf/--rhocall_frequency_file Frequency file for bcftools roh calculation (defaults to "GRCh37_anon_swegen_snp_-2016-10-19-.tab.gz", tab sep)
    -pvt/--pvt VT decompose and normalize (defaults to "0" (=no))
      -vtdec/--vt_decompose Split multi allelic records into single records (defaults to "1" (=yes))
      -vtnor/--vt_normalize Normalize variants (defaults to "1" (=yes))
      -vtunq/--vt_uniq Remove variant duplicates (defaults to "1" (=yes))
      -vtmaa/--vt_missing_alt_allele Remove missing alternative alleles '*' (defaults to "1" (=yes))
      -vtgmf/--vt_genmod_filter Remove common variants from vcf file (defaults to "1" (=yes))
      -vtgfr/--vt_genmod_filter_1000g Genmod annotate 1000G reference (defaults to "GRCh37_all_wgs_-phase3_v5b.2013-05-02-.vcf.gz")
      -vtmaf/--vt_genmod_filter_max_af Annotate MAX_AF from reference (defaults to "")
      -vtgft/--vt_genmod_filter_threshold Threshold for filtering variants (defaults to "0.10")
    -pvep/--pvarianteffectpredictor Annotate variants using VEP (defaults to "0" (=no))
      -vepp/--vep_directory_path Path to VEP script directory (defaults to "")
      -vepc/--vep_directory_cache Specify the cache directory to use (defaults to "")
      -vepr/--vep_reference Use Human reference file with VEP (defaults to "0" (=no))
      -vepf/--vep_features VEP features (defaults to ("hgvs","symbol","numbers","sift","polyphen","humdiv","domains","protein","ccds","uniprot","biotype","regulatory", "tsl", "canonical", "appris"); comma sep)
      -veppl/--vep_plugins VEP plugins (defaults to ("UpDownDistance, LoFtool, LoF"); comma sep)
    -pvcp/--pvcfparser Parse variants using vcfParser.pl (defaults to "0" (=no))
      -vcpvt/--vcfparser_vep_transcripts Parse VEP transcript specific entries (defaults to "0" (=no))
      -vcprff/--vcfparser_range_feature_file Range annotations file (defaults to ""; tab-sep)
      -vcprfa/--vcfparser_range_feature_annotation_columns Range annotations feature columns (defaults to ""; comma sep)
      -vcpsf/--vcfparser_select_file File containging list of genes to analyse seperately (defaults to "";tab-sep file and HGNC Symbol required)
      -vcpsfm/--vcfparser_select_file_matching_column Position of HGNC symbol column in SelectFile (defaults to "")
      -vcpsfa/--vcfparser_select_feature_annotation_columns Feature columns to use in annotation (defaults to ""; comma sep)
    -psne/--psnpeff Variant annotation using snpEff (defaults to "0" (=no))
#snpEffAnn
      -snep/--snpeff_path Path to snpEff. Mandatory for use of snpEff (defaults to "")
      -sneann/--snpeff_ann Annotate variants using snpeff (defaults to "1" (=yes))
      -snegbv/--snpeff_genome_build_version snpeff genome build version (defaults to "GRCh37.75")
      -snesaf/--snpsift_annotation_files Annotation files to use with snpsift (default to (GRCh37_all_wgs_-phase3_v5b.2013-05-02-.vcf.gz=AF GRCh37_exac_reheader_-r0.3.1-.vcf.gz=AF GRCh37_anon-swegen_snp_-1000samples-.vcf.gz=AF GRCh37_anon-swegen_indel_-1000samples-.vcf.gz=AF); Hash flag i.e. --Flag key=value)
      -snesaoi/--snpsift_annotation_outinfo_key snpsift output INFO key (default to (GRCh37_all_wgs_-phase3_v5b.2013-05-02-.vcf=1000G GRCh37_exac_reheader_-r0.3.1-.vcf.gz=EXAC GRCh37_anon-swegen_snp_-1000samples-.vcf.gz=SWEREF GRCh37_anon-swegen_indel_-1000samples-.vcf.gz=SWEREF); Hash flag i.e. --Flag key=value)
      -snesdbnsfp/--snpsift_dbnsfp_file DbNSFP File (defaults to "GRCh37_dbnsfp_-v2.9-.txt.gz")
      -snesdbnsfpa/--snpsift_dbnsfp_annotations DbNSFP annotations to use with snpsift (defaults to ("SIFT_pred","Polyphen2_HDIV_pred","Polyphen2_HVAR_pred","GERP++_NR","GERP++_RS","phastCons100way_vertebrate"); comma sep)

    ##Rankvariant
    -prav/--prankvariant Ranking of annotated variants (defaults to "0" (=no))
      -ravgft/--genmod_models_family_type Use one of the known setups (defaults to "mip")
      -ravanr/--genmod_annotate_regions Use predefined gene annotation supplied with genmod for defining genes (defaults to "1" (=yes))
      -ravcad/--genmod_annotate_cadd_files CADD score files (defaults to ""; comma sep)
      -ravspi/--genmod_annotate_spidex_file Spidex database for alternative splicing (defaults to "")
      -ravwg/--genmod_models_whole_gene Allow compound pairs in intronic regions (defaults to "1" (=yes))
      -ravrpf/--genmod_models_reduced_penetrance_file File containg genes with reduced penetrance (defaults to "")
      -ravrm/--rank_model_file Rank model config file (defaults to "")

    -pevab/--pendvariantannotationblock End variant annotation block by concatenating files (defaults to "0" (=no))
      -ravbf/--rankvariant_binary_file Produce binary file from the rank variant chromosomal sorted vcfs (defaults to "1" (=yes))
      -evabrgf/--endvariantannotationblock_remove_genes_file Remove variants in hgnc_ids (defaults to "")

    ###Utility
    -pped/--ppeddy QC for familial-relationships and sexes (defaults to "0" (=no) )
    -pplink/--pplink QC for samples gender and relationship (defaults to "0" (=no) )
    -pvai/--pvariant_integrity QC for samples relationship (defaults to "0" (=no) )
    -pevl/--pevaluation Compare concordance with NIST data set (defaults to "0" (=no) )
      -evlnid/--nist_id NIST high-confidence sample_id (defaults to "NA12878")
      -evlnhc/--nist_high_confidence_call_set NIST high-confidence variant calls (defaults to "GRCh37_nist_hg001_-na12878_v2.19-.vcf")
      -evlnil/--nist_high_confidence_call_set_bed NIST high-confidence variant calls interval list (defaults to "GRCh37_nist_hg001_-na12878_v2.19-.bed")
    -pqcc/--pqccollect Collect QC metrics from programs processed (defaults to "0" (=no) )
      -qccsi/--qccollect_sampleinfo_file SampleInfo file containing info on what to parse from this analysis run (defaults to "{outdata_dir}/{family_id}/{family_id}_qc_sample_info.yaml")
      -qccref/--qccollect_regexp_file Regular expression file containing the regular expression to be used for each program (defaults to "qc_regexp_-v1.13-.yaml")
      -qccske/--qccollect_skip_evaluation Skip evaluation step in qccollect (boolean)
    -pmqc/--pmultiqc Create aggregate bioinformatics analysis report across many samples (defaults to "0" (=no))
    -prem/--premoveredundantfiles Generating sbatch script for deletion of redundant files (defaults to "0" (=no);Note: Must be submitted manually to SLURM)
    -pars/--panalysisrunstatus Sets the analysis run status flag to finished in sample_info_file (defaults to "0" (=no))
    -psac/--psacct Generating sbatch script for SLURM info on each submitted job (defaults to "0" (=no))
    -sacfrf/--sacct_format_fields Format and fields of sacct output (defaults to "jobid", "jobname%50", "account", "partition", "alloccpus", "TotalCPU", "elapsed", "start", "end", "state", "exitcode")
END_USAGE
}

sub msacct {

##msacct

##Function : Output SLURM info on each job via sacct command
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, $family_id_ref
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $family_id_ref              => Family id {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use MIP::Workloadmanager::Slurm qw(slurm_sacct);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_chain_job_ids_dependency_add_to_path);

    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory     => lc($program_name),
            core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
        }
    );

    slurm_sacct(
        {
            fields_format_ref =>
              \@{ $active_parameter_href->{sacct_format_fields} },
            job_ids_ref => \@{ $job_id_href->{PAN}{PAN} },
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say $FILEHANDLE "\n";

    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        slurm_submit_chain_job_ids_dependency_add_to_path(
            {
                job_id_href         => $job_id_href,
                path                => $job_id_chain,
                log                 => $log,
                sbatch_file_name    => $file_path,
                job_dependency_type => q{afterany},
            }
        );
    }
}

sub analysisrunstatus {

##analysisrunstatus

##Function : Execute last in MAIN chain, tests that all recorded files exists, have a file sixe greater than zero, checks QC-metrics for PASS or FAIL and sets analysis run status flag to finished.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, $family_id_ref,
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $family_id_ref              => Family id {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_chain_job_ids_dependency_add_to_path);

    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory     => lc($program_name),
            core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
        }
    );

    say $FILEHANDLE q?STATUS="0"?
      ;  #Set status flagg so that perl not_finished remains in sample_info_file

    ###Test all file that are supposed to exists as they are present in the sample_info file
    my @paths_ref;

    ## Collects all programs file path(s) created by MIP located in %sample_info
    collect_path_entries(
        {
            sample_info_href => $sample_info_href,
            paths_ref        => \@paths_ref,
        }
    );

    print {$FILEHANDLE} q?readonly FILES=(?;    #Create bash array
    foreach my $path (@paths_ref) {

        if ( defined($path) )
        { #First analysis and dry run will otherwise cause try to print uninitialized values

            print {$FILEHANDLE} q?"? . $path . q?" ?;    #Add to array
        }
    }
    say $FILEHANDLE ")";                           #Close bash array
    say $FILEHANDLE q?for file in "${FILES[@]}"?;  #loop over files
    say $FILEHANDLE "do ";                         #for each element in array do
    say $FILEHANDLE "\t"
      . q?if [ -s "$file" ]; then?;    #file exists and is larger than zero
    say $FILEHANDLE "\t\t" . q?echo "Found file $file"?;    #Echo
    say $FILEHANDLE "\t" . q?else?;
    say $FILEHANDLE "\t\t"
      . q?echo "Could not find $file" >&2?;                 #Redirect to STDERR
    say $FILEHANDLE "\t\t"
      . q?STATUS="1"?
      ;   #Set status flagg so that perl notFinished remains in sample_info_file
    say $FILEHANDLE "\t" . q?fi?;
    say $FILEHANDLE q?done ?, "\n";

    ## Test varianteffectpredictor fork status. If varianteffectpredictor is unable to fork it will prematurely end the analysis and we will lose variants.
    if (
        defined(
            $sample_info_href->{program}{varianteffectpredictor}{stderrfile}
              {path}
        )
      )
    {

        my $variant_effect_predictor_file = catfile(
            $sample_info_href->{program}{varianteffectpredictor}{stderrfile}
              {path},
        );

        print {$FILEHANDLE}
          q?if grep -q "WARNING Unable to fork" ?
          ;    #not output the matched text only return the exit status code
        say $FILEHANDLE $variant_effect_predictor_file . q?; then?;    #Infile
        say $FILEHANDLE "\t" . q?STATUS="1"?;    #Found pattern
        say $FILEHANDLE "\t"
          . q?echo "variant_effector_predictor fork status=FAILED for file: ?
          . $variant_effect_predictor_file
          . q?" >&2?;                            #Echo
        say $FILEHANDLE q?else?;                 #Infile is clean
        say $FILEHANDLE "\t"
          . q?echo "variant_effector_predictor fork status=PASSED for file: ?
          . $variant_effect_predictor_file
          . q?" >&2?;                            #Echo
        say $FILEHANDLE q?fi?, "\n";
    }

    ## Test if FAIL exists in qccollect file i.e. issues with samples e.g. Sex and seq data correlation, relationship etc
    if ( !$active_parameter_href->{qccollect_skip_evaluation} ) {

        if ( defined( $sample_info_href->{program}{qccollect}{outfile} ) ) {

            my $qccollect_file = catfile(
                $sample_info_href->{program}{qccollect}{outdirectory},
                $sample_info_href->{program}{qccollect}{outfile}
            );

            print {$FILEHANDLE}
              q?if grep -q "FAIL" ?
              ;    #not output the matched text only return the exit status code
            say $FILEHANDLE $qccollect_file . q?; then?;    #Infile
            say $FILEHANDLE "\t" . q?STATUS="1"?;           #Found pattern
            say $FILEHANDLE "\t"
              . q?echo "qccollect status=FAILED for file: ?
              . $qccollect_file
              . q?" >&2?;                                   #Echo
            say $FILEHANDLE q?else?;                        #Infile is clean
            say $FILEHANDLE "\t"
              . q?echo "qccollect status=PASSED for file: ?
              . $qccollect_file
              . q?" >&2?;                                   #Echo
            say $FILEHANDLE q?fi?, "\n";
        }
    }

    ## Test integrity of vcf data keys in header and body
    if (   ( defined( $sample_info_href->{vcf_file}{clinical}{path} ) )
        || ( defined( $sample_info_href->{vcf_file}{research}{path} ) )
        || ( defined( $sample_info_href->{sv_vcf_file}{clinical}{path} ) )
        || ( defined( $sample_info_href->{sv_vcf_file}{research}{path} ) ) )
    {

        print {$FILEHANDLE} q?perl -MTest::Harness -e ' ?;    #Execute on cmd
        print {$FILEHANDLE} q?my %args = (?; #Adjust arguments to harness object
        print {$FILEHANDLE}
          q?verbosity => 1, ?;    #Print individual test results to STDOUT
        print {$FILEHANDLE} q?test_args => { ?;    #Argument to test script

        if ( defined( $sample_info_href->{vcf_file}{clinical}{path} ) ) {

            print {$FILEHANDLE}
              q?"test select file" => [ ?; #Add test for select file using alias
            print {$FILEHANDLE} q?"?
              . $sample_info_href->{vcf_file}{clinical}{path}
              . q?", ?;                    #Infile
            print {$FILEHANDLE} q?"?
              . $active_parameter_href->{config_file_analysis}
              . q?", ?;                    #ConfigFile
            print {$FILEHANDLE} q?], ?;
        }

        if ( defined( $sample_info_href->{vcf_file}{research}{path} ) ) {

            print {$FILEHANDLE}
              q?"test research file" => [ ?; #Add test research file using alias
            print {$FILEHANDLE} q?"?
              . $sample_info_href->{vcf_file}{research}{path}
              . q?", ?;                      #Infile
            print {$FILEHANDLE} q?"?
              . $active_parameter_href->{config_file_analysis}
              . q?", ?;                      #ConfigFile
            print {$FILEHANDLE} q?], ?;
        }
        if ( defined( $sample_info_href->{sv_vcf_file}{clinical}{path} ) ) {

            print {$FILEHANDLE}
              q?"test sv select file" => [ ?
              ;    #Add test for select file using alias
            print {$FILEHANDLE} q?"?
              . $sample_info_href->{vcf_file}{clinical}{path}
              . q?", ?;    #Infile
            print {$FILEHANDLE} q?"?
              . $active_parameter_href->{config_file_analysis}
              . q?", ?;    #ConfigFile
            print {$FILEHANDLE} q?], ?;
        }

        if ( defined( $sample_info_href->{sv_vcf_file}{research}{path} ) ) {

            print {$FILEHANDLE}
              q?"test sv research file" => [ ?
              ;            #Add test research file using alias
            print {$FILEHANDLE} q?"?
              . $sample_info_href->{vcf_file}{research}{path}
              . q?", ?;    #Infile
            print {$FILEHANDLE} q?"?
              . $active_parameter_href->{config_file_analysis}
              . q?", ?;    #ConfigFile
            print {$FILEHANDLE} q?], ?;
        }

        print {$FILEHANDLE} q?}); ?;
        print {$FILEHANDLE}
          q?my $harness = TAP::Harness->new( \%args ); ?
          ;                #Create harness using arguments provided
        print {$FILEHANDLE} q?$harness->runtests( ?;    #Execute test(s)

        if ( defined( $sample_info_href->{vcf_file}{clinical}{path} ) ) {

            print {$FILEHANDLE} q?["?
              . catfile( $Bin, "t", "mip_analysis.t" )
              . q?", "test select file"], ?;
        }

        if ( defined( $sample_info_href->{vcf_file}{research}{path} ) ) {

            print {$FILEHANDLE} q?["?
              . catfile( $Bin, "t", "mip_analysis.t" )
              . q?", "test research file"], ?;
        }
        if ( defined( $sample_info_href->{sv_vcf_file}{clinical}{path} ) ) {

            print {$FILEHANDLE} q?["?
              . catfile( $Bin, "t", "mip_analysis.t" )
              . q?", "test sv select file"], ?;
        }

        if ( defined( $sample_info_href->{sv_vcf_file}{research}{path} ) ) {

            print {$FILEHANDLE} q?["?
              . catfile( $Bin, "t", "mip_analysis.t" )
              . q?", "test sv research file"], ?;
        }
        print {$FILEHANDLE} q?)'?;
        say $FILEHANDLE "\n";
    }

    say $FILEHANDLE q?if [ $STATUS -ne 1 ]; then?;    #eval status flag
    say $FILEHANDLE "\t"
      . q?perl -i -p -e 'if($_=~/analysisrunstatus\:/) { s/not_finished/finished/g }' ?
      . $active_parameter_href->{sample_info_file} . q? ?;
    say $FILEHANDLE q?else?;    #Found discrepancies - exit
    say $FILEHANDLE "\t" . q?exit 1?;
    say $FILEHANDLE q?fi?, "\n";

    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        slurm_submit_chain_job_ids_dependency_add_to_path(
            {
                job_id_href      => $job_id_href,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
    return;
}

sub removeredundantfiles {

##removeredundantfiles

##Function : Generates a sbatch script, which removes redundant files.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $lane_href, $program_name, family_id_ref, $outaligner_dir_ref, $call_type
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $lane_href                  => The lane info hash {REF}
##         : $program_name               => Program name
##         : $family_id_ref              => Family id {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $outaligner_dir_ref;
    my $call_type;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $lane_href    = $arg_href->{lane_href};
    my $program_name = $arg_href->{program_name};

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        lane_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$lane_href
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);

    my $reduce_io_ref = \$active_parameter_href->{reduce_io};
    my $xargs_file_path_prefix;

    ## Filehandles
    my $FILEHANDLE      = IO::Handle->new();    #Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory     => $$outaligner_dir_ref,
            core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
        }
    );

    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Sample files
        ##Removes intermediate files from the MIP analysis depending on set MIP parameters
        remove_redundant_files(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                lane_href               => $lane_href,
                FILEHANDLE              => $FILEHANDLE,
                sample_id               => $sample_id,
                reduce_io_ref           => $reduce_io_ref,
                outaligner_dir_ref      => $outaligner_dir_ref,
            }
        );
    }

    ## Family files
    ##Removes intermediate files from the MIP analysis depending on set MIP parameters
    remove_redundant_files(
        {
            parameter_href          => $parameter_href,
            active_parameter_href   => $active_parameter_href,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            sample_info_href        => $sample_info_href,
            file_info_href          => $file_info_href,
            lane_href               => $lane_href,
            FILEHANDLE              => $FILEHANDLE,
            reduce_io_ref           => $reduce_io_ref,
            outaligner_dir_ref      => $outaligner_dir_ref,
        }
    );
    close $FILEHANDLE;
}

sub mmultiqc {

##mmultiqc

##Function : Aggregate bioinforamtics reports per case
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, $family_id_ref,
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $family_id_ref              => Family id {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use Program::Qc::Multiqc qw(multiqc);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_chain_job_ids_dependency_add_to_path);

    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory     => lc($program_name),
            core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
        }
    );

    ## Assign directories
    my $program_outdirectory_name =
      $parameter_href->{ "p" . $program_name }{outdir_name};

    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Assign directories
        my $insample_directory =
          catdir( $active_parameter_href->{outdata_dir}, $sample_id );
        my $outsample_directory = catdir( $active_parameter_href->{outdata_dir},
            $sample_id, $program_outdirectory_name );

        multiqc(
            {
                indir_path  => $insample_directory,
                outdir_path => $outsample_directory,
                force       => 1,
                FILEHANDLE  => $FILEHANDLE,
            }
        );
        say $FILEHANDLE "\n";
    }

    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        slurm_submit_chain_job_ids_dependency_add_to_path(
            {
                job_id_href      => $job_id_href,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
}

sub mqccollect {

##mqccollect

##Function : Collect qc metrics for this analysis run.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, family_id_ref, $call_type,
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $family_id_ref              => Family id {REF}
##         : $call_type                  => Variant call type

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $call_type;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use Program::Qc::Mip qw(qccollect);
    use MIP::QC::Record qw(add_program_outfile_to_sample_info);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_chain_job_ids_dependency_add_to_path);

    my $reduce_io_ref = \$active_parameter_href->{reduce_io};
    my $job_id_chain  = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory     => lc($program_name),
            core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
        }
    );

    ## Assign directories
    my $outfamily_directory =
      catdir( $active_parameter_href->{outdata_dir}, $$family_id_ref );

    qccollect(
        {
            infile_path  => $active_parameter_href->{qccollect_sampleinfo_file},
            outfile_path => catfile(
                $outfamily_directory, $$family_id_ref . "_qc_metrics.yaml"
            ),
            log_file_path => catfile(
                $outfamily_directory, $$family_id_ref . "_qccollect.log"
            ),
            regexp_file_path => $active_parameter_href->{qccollect_regexp_file},
            skip_evaluation =>
              $active_parameter_href->{qccollect_skip_evaluation},
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say $FILEHANDLE "\n";

    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        ## Collect QC metadata info for later use
        my $qc_metric_outfile = $$family_id_ref . q{_qc_metrics.yaml};
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => 'qccollect',
                outdirectory     => $outfamily_directory,
                outfile          => $qc_metric_outfile,
                path => catfile( $outfamily_directory, $qc_metric_outfile ),
            }
        );

        slurm_submit_chain_job_ids_dependency_add_to_path(
            {
                job_id_href      => $job_id_href,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
}

sub endvariantannotationblock {

##endvariantannotationblock

##Function : Concatenate ouput from variant annotation block.
##Returns  : "|$xargs_file_counter"
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, $program_info_path, $file_path, $FILEHANDLE, family_id_ref, $temp_directory_ref, $reference_dir_ref, $outaligner_dir_ref, $call_type, $xargs_file_counter
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $program_info_path          => The program info path
##         : $file_path                  => File path
##         : $FILEHANDLE                 => Sbatch filehandle to write to
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $reference_dir_ref          => MIP reference directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type
##         : $xargs_file_counter         => The xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $reference_dir_ref;
    my $outaligner_dir_ref;
    my $call_type;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;
    my $program_info_path;
    my $file_path;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        file_path         => { strict_type => 1, store => \$file_path },
        FILEHANDLE    => { store => \$FILEHANDLE, },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        reference_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{reference_dir},
            strict_type => 1,
            store       => \$reference_dir_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Set::File qw{set_file_suffix};
    use MIP::Get::File qw{get_file_suffix};
    use MIP::Delete::List qw{ delete_contig_elements };
    use MIP::IO::Files qw{ xargs_migrate_contig_files };
    use Program::Htslib qw(bgzip tabix);
    use MIP::Gnu::Software::Gnu_grep qw(gnu_grep);
    use MIP::QC::Record qw(add_program_metafile_to_sample_info);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_family);

    my $reduce_io_ref = \$active_parameter_href->{reduce_io};
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $xargs_file_path_prefix;
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};
    my $vcfparser_analysis_type  = "";
    my $contigs_size_ordered_ref = \@{ $file_info_href->{contigs_size_ordered} }
      ;    #Set default for size ordered contigs
    my @contigs = @{ $file_info_href->{contigs} };    #Set default for contigs

    ## Set the number of cores
    my $core_number =
      $active_parameter_href->{module_core_number}{ "p" . $program_name };

    ## Filehandles
    my $XARGSFILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    unless ( defined($FILEHANDLE) ) {           #Run as individual sbatch script

        $FILEHANDLE = IO::Handle->new();        #Create anonymous filehandle

        use MIP::Script::Setup_script qw(setup_script);
        use MIP::IO::Files qw(migrate_file);

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        ( $file_path, $program_info_path ) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                job_id_href           => $job_id_href,
                FILEHANDLE            => $FILEHANDLE,
                directory_id          => $$family_id_ref,
                program_name          => $program_name,
                program_directory     => catfile( lc($$outaligner_dir_ref) ),
                core_number           => $core_number,
                process_time =>
                  $active_parameter_href->{module_time}{ "p" . $program_name },
                temp_directory => $$temp_directory_ref
            }
        );
    }

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    my $outfamily_file_directory =
      catfile( $active_parameter_href->{outdata_dir}, $$family_id_ref );

    ## Assign file_tags
    my $infile_tag = $file_info_href->{$$family_id_ref}{prankvariant}{file_tag};
    my $infile_prefix = $$family_id_ref . $infile_tag . $call_type;
    my $file_path_prefix = catfile( $$temp_directory_ref, $infile_prefix );
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{prankvariant}{file_tag};
    my $outfile_prefix      = $$family_id_ref . $outfile_tag . $call_type;
    my $outfile_path_prefix = catfile( $$temp_directory_ref,
        $$family_id_ref . $outfile_tag . $call_type );
    my $final_path_prefix = catfile( $outfamily_directory, $outfile_prefix );

    ## Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            jobid_chain    => $job_id_chain,
        }
    );

    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    for (
        my $vcfparser_outfile_counter = 0 ;
        $vcfparser_outfile_counter <
        $active_parameter_href->{vcfparser_outfile_count} ;
        $vcfparser_outfile_counter++
      )
    {

        if ( $vcfparser_outfile_counter == 1 ) {

            $vcfparser_analysis_type = ".selected";    #SelectFile variants
            $contigs_size_ordered_ref =
              \@{ $file_info_href->{sorted_select_file_contigs} }
              ;                                        #Selectfile contigs
            @contigs = @{ $file_info_href->{select_file_contigs} };

            if ( $consensus_analysis_type eq "wes" )
            {    #Remove MT|M since no exome kit so far has mitochondrial probes

                ## Removes an element from array and return new array while leaving orginal elements_ref untouched
                @contigs = delete_contig_elements(
                    {
                        elements_ref =>
                          \@{ $file_info_href->{select_file_contigs} },
                        remove_contigs_ref => [qw{ MT M }],
                    }
                );
            }
        }

        if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

            ## Copy file(s) to temporary directory
            say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
            $xargs_file_counter = xargs_migrate_contig_files(
                {
                    FILEHANDLE         => $FILEHANDLE,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    contigs_ref        => $contigs_size_ordered_ref,
                    file_path          => $file_path,
                    program_info_path  => $program_info_path,
                    core_number        => $core_number,
                    xargs_file_counter => $xargs_file_counter,
                    infile             => $infile_prefix,
                    file_ending        => $vcfparser_analysis_type
                      . $infile_suffix . "*",
                    indirectory    => $infamily_directory,
                    temp_directory => $active_parameter_href->{temp_directory},
                }
            );
        }

        ## Writes sbatch code to supplied filehandle to concatenate variants in vcf format. Each array element is combined with the infile prefix and postfix.
        concatenate_variants(
            {
                active_parameter_href => $active_parameter_href,
                FILEHANDLE            => $FILEHANDLE,
                elements_ref          => \@contigs,
                infile_prefix         => $file_path_prefix . "_",
                infile_postfix => $vcfparser_analysis_type . $infile_suffix,
                outfile        => $outfile_path_prefix
                  . $vcfparser_analysis_type
                  . $outfile_suffix,
            }
        );

        ## Remove variants in hgnc_id list from vcf
        if ( $active_parameter_href
            ->{endvariantannotationblock_remove_genes_file} )
        {

            ## Removes contig_names from contigs array if no male or other found
            gnu_grep(
                {
                    filter_file_path => catfile(
                        $$reference_dir_ref,
                        $active_parameter_href->{sv_reformat_remove_genes_file}
                    ),
                    infile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $outfile_suffix,
                    outfile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . "_filtered"
                      . $outfile_suffix,
                    invert_match => 1,
                    FILEHANDLE   => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} "\n";

            if ( $vcfparser_outfile_counter == 1 ) {

                $sample_info_href->{program}{$program_name}
                  {reformat_remove_genes_file}{clinical}{path} =
                    $final_path_prefix
                  . $vcfparser_analysis_type
                  . "_filtered"
                  . $outfile_suffix;    #Save filtered file
            }
            else {

                $sample_info_href->{program}{$program_name}
                  {reformat_remove_genes_file}{research}{path} =
                    $final_path_prefix
                  . $vcfparser_analysis_type
                  . "_filtered"
                  . $outfile_suffix;    #Save filtered file
            }

            ## Copies file from temporary directory.
            say {$FILEHANDLE} q{## Copy file from temporary directory};
            migrate_file(
                {
                    infile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . q{_filtered}
                      . $outfile_suffix,
                    outfile_path => $outfamily_directory,
                    FILEHANDLE   => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} q{wait}, "\n";
        }

        if ( $active_parameter_href->{rankvariant_binary_file} ) {

            ## Compress or decompress original file or stream to outfile (if supplied)
            bgzip(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $outfile_suffix,
                    outfile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $outfile_suffix . ".gz",
                    write_to_stdout => 1,
                }
            );
            say {$FILEHANDLE} "\n";

            ## Index file using tabix
            tabix(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $outfile_suffix . ".gz",
                    force  => 1,
                    preset => substr( $outfile_suffix, 1 ),
                }
            );
            say {$FILEHANDLE} "\n";
        }

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                infile_path => $outfile_path_prefix
                  . $vcfparser_analysis_type
                  . $outfile_suffix . q{*},
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} q{wait}, "\n";

        ## Adds the most complete vcf file to sample_info
        add_most_complete_vcf(
            {
                active_parameter_href => $active_parameter_href,
                sample_info_href      => $sample_info_href,
                program_name          => $program_name,
                path                  => $final_path_prefix
                  . $vcfparser_analysis_type
                  . $outfile_suffix,
                vcfparser_outfile_counter => $vcfparser_outfile_counter,
            }
        );

        if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

            if ( $vcfparser_outfile_counter == 1 ) {

                # Save clinical candidate list path
                my $clinical_candidate_path =
                    $final_path_prefix
                  . $vcfparser_analysis_type
                  . $outfile_suffix;
                add_program_metafile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        program_name     => $program_name,
                        metafile_tag     => q{clinical},
                        path             => $clinical_candidate_path,
                    }
                );

                if ( $active_parameter_href->{rankvariant_binary_file} ) {

                    $sample_info_href->{vcf_binary_file}{clinical}{path} =
                        $final_path_prefix
                      . $vcfparser_analysis_type
                      . $outfile_suffix . ".gz";
                }
            }
            else {

                # Save research candidate list path
                my $research_candidate_path =
                    $final_path_prefix
                  . $vcfparser_analysis_type
                  . $outfile_suffix;
                add_program_metafile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        program_name     => $program_name,
                        metafile_tag     => q{research},
                        path             => $research_candidate_path,
                    }
                );

                if ( $active_parameter_href->{rankvariant_binary_file} ) {

                    $sample_info_href->{vcf_binary_file}{research}{path} =
                        $final_path_prefix
                      . $vcfparser_analysis_type
                      . $outfile_suffix . ".gz";
                }
            }
        }
    }

    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $$family_id_ref,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
    if ($$reduce_io_ref) {

        return
          $xargs_file_counter
          ; #Track the number of created xargs scripts per module for Block algorithm
    }
}

sub rankvariant {

##rankvariant

##Function : Annotate and score variants depending on mendelian inheritance, frequency and phenotype etc.
##Returns  : "|$xargs_file_counter"
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, $program_info_path, $file_path, $FILEHANDLE, family_id_ref, $temp_directory_ref, $outaligner_dir_ref, $call_type, $xargs_file_counter
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $program_info_path          => The program info path
##         : $file_path                  => File path
##         : $FILEHANDLE                 => Sbatch filehandle to write to
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type
##         : $xargs_file_counter         => The xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
    my $call_type;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;
    my $program_info_path;
    my $file_path;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        file_path         => { strict_type => 1, store => \$file_path },
        FILEHANDLE    => { store => \$FILEHANDLE, },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw(get_core_number);
    use MIP::Set::File qw{set_file_suffix};
    use MIP::Get::File qw{get_file_suffix};
    use MIP::Delete::List qw{ delete_contig_elements };
    use MIP::IO::Files qw{ xargs_migrate_contig_files };
    use MIP::Recipes::Xargs qw{ xargs_command };
    use Program::Variantcalling::Genmod qw(annotate models score compound);
    use MIP::QC::Record
      qw(add_program_outfile_to_sample_info add_program_metafile_to_sample_info);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_family);

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $reduce_io_ref = \$active_parameter_href->{reduce_io};
    my $xargs_file_path_prefix;
    my $vcfparser_analysis_type = "";
    my @contigs_size_ordered    = @{ $file_info_href->{contigs_size_ordered} }
      ;    #Set default for size ordered contigs
    my @contigs = @{ $file_info_href->{contigs} };    #Set default for contigs
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $XARGSFILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Get core number depending on user supplied input exists or not and max number of cores
    my $core_number = get_core_number(
        {
            module_core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            modifier_core_number => scalar(@contigs),
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    ### Detect the number of cores to use per genmod process.
    ## Limit number of cores requested to the maximum number of cores available per node
    my $genmod_core_number = check_max_core_number(
        {
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
            core_number_requested => 16,
        }
    );

    unless ( defined($FILEHANDLE) ) {    #Run as individual sbatch script

        $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

        use MIP::Script::Setup_script qw(setup_script);
        use MIP::IO::Files qw(migrate_file);

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        ( $file_path, $program_info_path ) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                job_id_href           => $job_id_href,
                FILEHANDLE            => $FILEHANDLE,
                directory_id          => $$family_id_ref,
                program_name          => $program_name,
                program_directory     => catfile( lc($$outaligner_dir_ref) ),
                core_number           => $core_number,
                process_time =>
                  $active_parameter_href->{module_time}{ "p" . $program_name },
                temp_directory => $$temp_directory_ref
            }
        );
    }

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    my $outfamily_file_directory =
      catfile( $active_parameter_href->{outdata_dir}, $$family_id_ref );

    ## Assign file_tags
    my $infile_tag    = $file_info_href->{$$family_id_ref}{psnpeff}{file_tag};
    my $infile_prefix = $$family_id_ref . $infile_tag . $call_type;
    my $file_path_prefix = catfile( $$temp_directory_ref, $infile_prefix );
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{ "p" . $program_name }{file_tag};
    my $outfile_prefix = $$family_id_ref . $outfile_tag . $call_type;
    my $outfile_path_prefix = catfile( $$temp_directory_ref, $outfile_prefix );

    ## Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            jobid_chain    => $job_id_chain,
        }
    );

    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    my $family_file =
      catfile( $outfamily_file_directory, $$family_id_ref . ".fam" );

    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            parameter_href        => $parameter_href,
            active_parameter_href => $active_parameter_href,
            sample_info_href      => $sample_info_href,
            FILEHANDLE            => $FILEHANDLE,
            fam_file_path         => $family_file,
        }
    );

    for (
        my $vcfparser_outfile_counter = 0 ;
        $vcfparser_outfile_counter <
        $active_parameter_href->{vcfparser_outfile_count} ;
        $vcfparser_outfile_counter++
      )
    {

        if ( $vcfparser_outfile_counter == 1 ) {

            $vcfparser_analysis_type = ".selected";    #SelectFile variants
            @contigs_size_ordered =
              @{ $file_info_href->{sorted_select_file_contigs} }
              ;                                        #Selectfile contigs
            @contigs = @{ $file_info_href->{select_file_contigs} };

            if ( $consensus_analysis_type eq "wes" )
            {    #Remove MT|M since no exome kit so far has mitochondrial probes

                ## Removes an element from array and return new array while leaving orginal elements_ref untouched
                @contigs = delete_contig_elements(
                    {
                        elements_ref =>
                          \@{ $file_info_href->{select_file_contigs} },
                        remove_contigs_ref => [q{ MT M }],
                    }
                );

                ## Removes an element from array and return new array while leaving orginal elements_ref untouched
                @contigs_size_ordered = delete_contig_elements(
                    {
                        elements_ref =>
                          \@{ $file_info_href->{sorted_select_file_contigs} },
                        remove_contigs_ref => [qw{ MT M }],
                    }
                );
            }
        }

        if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

            ## Copy file(s) to temporary directory
            say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
            $xargs_file_counter = xargs_migrate_contig_files(
                {
                    FILEHANDLE         => $FILEHANDLE,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    contigs_ref        => \@contigs_size_ordered,
                    file_path          => $file_path,
                    program_info_path  => $program_info_path,
                    core_number        => $core_number,
                    xargs_file_counter => $xargs_file_counter,
                    infile             => $infile_prefix,
                    file_ending        => $vcfparser_analysis_type
                      . $infile_suffix . "*",
                    indirectory    => $infamily_directory,
                    temp_directory => $active_parameter_href->{temp_directory},
                }
            );
        }

        ## Genmod
        say {$FILEHANDLE} "## GeneMod";

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                FILEHANDLE         => $FILEHANDLE,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                file_path          => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $genmod_core_number,
                xargs_file_counter => $xargs_file_counter,
            }
        );

        my $genmod_module = "";   #Track which genmod modules has been processed

        ## Process per contig
        while ( my ( $contig_index, $contig ) = each(@contigs_size_ordered) ) {

            ## Get parameters
            $genmod_module = "";    #Restart for next contig
            my $genmod_indata =
                $file_path_prefix . "_"
              . $contig
              . $vcfparser_analysis_type
              . $infile_suffix;     #InFile

            my $genmod_outfile_path;
            ## Check affected/unaffected status
            if (
                ( defined( $parameter_href->{dynamic_parameter}{unaffected} ) )
                && ( @{ $parameter_href->{dynamic_parameter}{unaffected} } eq
                    @{ $active_parameter_href->{sample_ids} } )
              )
            {                       #Only unaffected

                if ( ( !$contig_index ) && ( !$vcfparser_outfile_counter ) ) {

                    $log->warn(
"Only unaffected sample in pedigree - skipping genmod 'models', 'score' and 'compound'"
                    );
                }

                ## Output file
                $genmod_outfile_path =
                    $outfile_path_prefix . "_"
                  . $contig
                  . $vcfparser_analysis_type
                  . $infile_suffix;

            }
            else {

                ## Output stream
                $genmod_outfile_path =
                  catfile( dirname( devnull() ), "stdout" );
            }

            ## Genmod Annotate
            $genmod_module = "_annotate";

            Program::Variantcalling::Genmod::annotate(
                {
                    cadd_file_paths_ref =>
                      \@{ $active_parameter_href->{genmod_annotate_cadd_files}
                      },
                    infile_path     => $genmod_indata,
                    outfile_path    => $genmod_outfile_path,
                    stderrfile_path => $xargs_file_path_prefix . "."
                      . $contig
                      . $genmod_module
                      . ".stderr.txt",
                    spidex_file_path =>
                      $active_parameter_href->{genmod_annotate_spidex_file},
                    verbosity           => "v",
                    temp_directory_path => $$temp_directory_ref,
                    annotate_region =>
                      $active_parameter_href->{genmod_annotate_regions},
                    FILEHANDLE => $XARGSFILEHANDLE,
                }
            );

            if (
                ( defined( $parameter_href->{dynamic_parameter}{unaffected} ) )
                && ( @{ $parameter_href->{dynamic_parameter}{unaffected} } eq
                    @{ $active_parameter_href->{sample_ids} } )
              )
            {    #Only unaffected - do nothing
            }
            else {

                print $XARGSFILEHANDLE "| ";    #Pipe

                $genmod_indata = "-";           #Preparation for next module

                ## Genmod Models
                $genmod_module .= "_models";

                my $use_vep;
                if (   ( $active_parameter_href->{pvarianteffectpredictor} > 0 )
                    && ( !$active_parameter_href->{genmod_annotate_regions} ) )
                {    #Use VEP annotations in compound models

                    $use_vep = 1;
                }
                models(
                    {
                        infile_path => $genmod_indata,
                        outfile_path =>
                          catfile( dirname( devnull() ), "stdout" ),
                        stderrfile_path => $xargs_file_path_prefix . "."
                          . $contig
                          . $genmod_module
                          . ".stderr.txt",
                        family_file => $family_file,
                        family_type =>
                          $active_parameter_href->{genmod_models_family_type},
                        reduced_penetrance_file_path => $active_parameter_href
                          ->{genmod_models_reduced_penetrance_file},
                        thread_number => 4,
                        vep           => $use_vep,
                        whole_gene =>
                          $active_parameter_href->{genmod_models_whole_gene},
                        verbosity           => "v",
                        temp_directory_path => $$temp_directory_ref,
                        FILEHANDLE          => $XARGSFILEHANDLE,
                    }
                );
                print $XARGSFILEHANDLE "| ";    #Pipe

                ## Genmod Score
                $genmod_module .= "_score";
                score(
                    {
                        infile_path => $genmod_indata,
                        outfile_path =>
                          catfile( dirname( devnull() ), "stdout" ),
                        stderrfile_path => $xargs_file_path_prefix . "."
                          . $contig
                          . $genmod_module
                          . ".stderr.txt ",
                        verbosity   => "v",
                        family_file => $family_file,
                        family_type =>
                          $active_parameter_href->{genmod_models_family_type},
                        rank_result => 1,
                        rank_model_file_path =>
                          $active_parameter_href->{rank_model_file},
                        FILEHANDLE => $XARGSFILEHANDLE,
                    }
                );
                print $XARGSFILEHANDLE "| ";    #Pipe

                ##Genmod Compound
                $genmod_module .= "_compound";

                compound(
                    {
                        infile_path  => $genmod_indata,
                        outfile_path => $outfile_path_prefix . "_"
                          . $contig
                          . $vcfparser_analysis_type
                          . $infile_suffix,
                        stderrfile_path => $xargs_file_path_prefix . "."
                          . $contig
                          . $genmod_module
                          . ".stderr.txt",
                        verbosity           => "v",
                        temp_directory_path => $$temp_directory_ref,
                        vep                 => $use_vep,
                        FILEHANDLE          => $XARGSFILEHANDLE,
                    }
                );
            }
            say {$XARGSFILEHANDLE} "\n";
        }

        if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

            ## Copies file from temporary directory. Per contig
            say {$FILEHANDLE} q{## Copy file from temporary directory};
            ($xargs_file_counter) = xargs_migrate_contig_files(
                {
                    FILEHANDLE        => $FILEHANDLE,
                    XARGSFILEHANDLE   => $XARGSFILEHANDLE,
                    contigs_ref       => \@contigs_size_ordered,
                    file_path         => $file_path,
                    program_info_path => $program_info_path,
                    core_number => $active_parameter_href->{max_cores_per_node},
                    xargs_file_counter => $xargs_file_counter,
                    outfile            => $outfile_prefix,
                    file_ending        => $vcfparser_analysis_type
                      . $outfile_suffix . "*",
                    outdirectory   => $outfamily_directory,
                    temp_directory => $$temp_directory_ref,
                }
            );
        }
        else {

            ## QC Data File(s)
            migrate_file(
                {
                    infile_path => $outfile_path_prefix . q{_}
                      . $file_info_href->{contigs_size_ordered}[0]
                      . $vcfparser_analysis_type
                      . $outfile_suffix,
                    outfile_path => $outfamily_directory,
                    FILEHANDLE   => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} q{wait}, "\n";
        }

        if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

            my $qc_genmod_outfile =
                $outfile_prefix . q{_}
              . $file_info_href->{contigs_size_ordered}[0]
              . $vcfparser_analysis_type
              . $outfile_suffix;
            add_program_outfile_to_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    program_name     => q{genmod},
                    outdirectory     => $outfamily_directory,
                    outfile          => $qc_genmod_outfile,
                }
            );

            # Add to Sample_info
            if ( defined( $active_parameter_href->{rank_model_file} ) ) {

                my $rank_model_version;
                if ( $active_parameter_href->{rank_model_file} =~
                    /v(\d+\.\d+.\d+|\d+\.\d+)/ )
                {

                    $rank_model_version = $1;
                }
                add_program_metafile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        program_name     => q{genmod},
                        metafile_tag     => q{rank_model},
                        file =>
                          basename( $active_parameter_href->{rank_model_file} ),
                        path    => $active_parameter_href->{rank_model_file},
                        version => $rank_model_version,
                    }
                );
            }
        }
    }
    if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

        close $FILEHANDLE;

        if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

            slurm_submit_job_sample_id_dependency_add_to_family(
                {
                    job_id_href             => $job_id_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    sample_ids_ref =>
                      \@{ $active_parameter_href->{sample_ids} },
                    family_id        => $$family_id_ref,
                    path             => $job_id_chain,
                    log              => $log,
                    sbatch_file_name => $file_path,
                }
            );
        }
    }
    if ($$reduce_io_ref) {

        return
          $xargs_file_counter
          ; #Track the number of created xargs scripts per module for Block algorithm
    }
}

sub gatk_variantevalexome {

##gatk_variantevalexome

##Function : GATK varianteval for exome variants.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $sample_id_ref, $program_name, family_id_ref, $temp_directory_ref, $outaligner_dir_ref, $call_type,
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $sample_id_ref              => Sample id {REF}
##         : $program_name               => Program name {REF}
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
    my $call_type;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $sample_id_ref;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        sample_id_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$sample_id_ref
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use MIP::IO::Files qw(migrate_file);
    use MIP::Language::Java qw{java_core};
    use MIP::Set::File qw{set_file_suffix };
    use MIP::Get::File qw{get_file_suffix get_merged_infile_prefix };
    use MIP::Gnu::Coreutils qw(gnu_cat gnu_sort);
    use Program::Variantcalling::Bedtools qw(intersectbed);
    use MIP::Program::Variantcalling::Gatk qw(gatk_varianteval);
    use MIP::QC::Record qw(add_program_outfile_to_sample_info);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_family_dead_end);

    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$sample_id_ref,
            program_name          => $program_name,
            program_directory =>
              catfile( lc($$outaligner_dir_ref), lc($program_name) ),
            call_type   => $call_type,
            core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
            temp_directory => $$temp_directory_ref
        }
    );

    ## Assign directories
    my $outsample_directory = catfile( $active_parameter_href->{outdata_dir},
        $$sample_id_ref, $$outaligner_dir_ref, lc($program_name) );
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $$sample_id_ref,
        }
    );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$$family_id_ref}{pgatk_combinevariantcallsets}
      {file_tag};
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{pgatk_combinevariantcallsets}
      {file_tag};
    my $infile_prefix       = $$family_id_ref . $infile_tag . $call_type;
    my $file_path_prefix    = catfile( $$temp_directory_ref, $infile_prefix );
    my $outfile_prefix      = $merged_infile_prefix . $outfile_tag;
    my $outfile_path_prefix = catfile( $$temp_directory_ref, $outfile_prefix );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            jobid_chain =>
              $parameter_href->{pgatk_combinevariantcallsets}{chain},
        }
    );
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_eval_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    my $extract_exonic_regexp =
      q?perl -ne ' if ( ($_=~/exonic/) || ($_=~/splicing/) ) {print $_;}' ?;

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => catfile(
                $infamily_directory, $infile_prefix . $infile_suffix . q{*}
            ),
            outfile_path => $$temp_directory_ref
        }
    );

    say {$FILEHANDLE} q{wait}, "\n";

    ## Select sample id from family id vcf file

    ## GATK SelectVariants
    say {$FILEHANDLE} "## GATK SelectVariants";

    gatk_selectvariants(
        {
            memory_allocation => q{Xmx2g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $$temp_directory_ref,
            java_jar       => catfile(
                $active_parameter_href->{gatk_path},
                "GenomeAnalysisTK.jar"
            ),
            sample_names_ref => [$$sample_id_ref],
            logging_level    => $active_parameter_href->{gatk_logging_level},
            referencefile_path =>
              $active_parameter_href->{human_genome_reference},
            infile_path  => $file_path_prefix . $infile_suffix,
            outfile_path => $outfile_path_prefix
              . $call_type . "_temp"
              . $infile_suffix,
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

    ## Update file_tags
    $infile_tag    = $file_info_href->{$$family_id_ref}{prankvariant}{file_tag};
    $infile_prefix = $$family_id_ref . $infile_tag . $call_type;
    $file_path_prefix = catfile( $$temp_directory_ref, $infile_prefix );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => catfile(
                $infamily_directory, $infile_prefix . $infile_suffix . q{*}
            ),
            outfile_path => $$temp_directory_ref
        }
    );
    say {$FILEHANDLE} q{wait}, "\n";

    ## Extract exonic variants
    say   {$FILEHANDLE} "## Extract exonic variants";
    print {$FILEHANDLE} $extract_exonic_regexp;
    print {$FILEHANDLE} $file_path_prefix . $infile_suffix . " ";    #InFile
    say   {$FILEHANDLE} "> "
      . catfile(
        $$temp_directory_ref,
        $$sample_id_ref
          . $infile_tag
          . $call_type
          . "_exonic_variants"
          . $infile_suffix
      ),
      "\n";                                                          #OutFile

    ## Include potential select file variants
    if ( $active_parameter_href->{vcfparser_outfile_count} == 2 ) {

        my $vcfparser_analysis_type = ".selected";    #SelectFile variants

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        migrate_file(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => catfile(
                    $infamily_directory,
                    $infile_prefix
                      . $vcfparser_analysis_type
                      . $infile_suffix . q{*}
                ),
                outfile_path => $$temp_directory_ref
            }
        );
        say {$FILEHANDLE} q{wait}, "\n";

        ## Extract exonic variants
        say   {$FILEHANDLE} "## Extract exonic variants";
        print {$FILEHANDLE} $extract_exonic_regexp;
        print {$FILEHANDLE} $file_path_prefix
          . $vcfparser_analysis_type
          . $infile_suffix
          . " ";    #InFile
        say {$FILEHANDLE} "> "
          . catfile(
            $$temp_directory_ref,
            $$sample_id_ref
              . $infile_tag
              . $call_type
              . "_exonic_variants"
              . $vcfparser_analysis_type
              . $infile_suffix
          ),
          "\n";     #OutFile

        ## Merge orphans and selectfiles
        say {$FILEHANDLE} "## Merge orphans and selectfile(s)";
        gnu_cat(
            {
                infile_paths_ref => [
                    catfile(
                        $$temp_directory_ref,
                        $$sample_id_ref
                          . $infile_tag
                          . $call_type
                          . "_exonic_variants"
                          . $infile_suffix
                    ),
                    catfile(
                        $$temp_directory_ref,
                        $$sample_id_ref
                          . $infile_tag
                          . $call_type
                          . "_exonic_variants"
                          . $vcfparser_analysis_type
                          . $infile_suffix
                    )
                ],
                outfile_path => catfile(
                    $$temp_directory_ref,
                    $$sample_id_ref
                      . $infile_tag
                      . $call_type
                      . "_exonic_variants_combined"
                      . $infile_suffix
                ),
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} "\n";

        ## Sort combined file
        say {$FILEHANDLE} "## Sort combined file";
        gnu_sort(
            {
                keys_ref    => [ "1,1", "2,2n" ],
                infile_path => catfile(
                    $$temp_directory_ref,
                    $$sample_id_ref
                      . $infile_tag
                      . $call_type
                      . "_exonic_variants_combined"
                      . $infile_suffix
                ),
                outfile_path => catfile(
                    $$temp_directory_ref,
                    $$sample_id_ref
                      . $infile_tag
                      . $call_type
                      . "_exonic_variants"
                      . $infile_suffix
                ),
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} "\n";
    }

    print {$FILEHANDLE} q?perl -ne ' if ($_=~/^#/) {print $_;}' ?;
    print {$FILEHANDLE} $file_path_prefix . $infile_suffix . " ";    #InFile
    print {$FILEHANDLE} "| ";                                        #Pipe
    gnu_cat(
        {
            infile_paths_ref => [
                "-",
                catfile(
                    $$temp_directory_ref,
                    $$sample_id_ref
                      . $infile_tag
                      . $call_type
                      . "_exonic_variants"
                      . $infile_suffix
                )
            ],
            outfile_path => catfile(
                $$temp_directory_ref,
                $$sample_id_ref
                  . $infile_tag
                  . $call_type
                  . "_exonic_variants_head"
                  . $infile_suffix
            ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

    ## Intersect exonic variants from created sample_id vcf file (required for gatk_varianteval for exonic variants)
    say {$FILEHANDLE}
      "## Intersect exonic variants from created sample_id vcf file";
    intersectbed(
        {
            with_header => 1,
            infile_path => $outfile_path_prefix
              . $call_type . "_temp"
              . $infile_suffix,
            intersectfile_path => catfile(
                $$temp_directory_ref,
                $$sample_id_ref
                  . $infile_tag
                  . $call_type
                  . "_exonic_variants_head"
                  . $infile_suffix
            ),
            outfile_path => $outfile_path_prefix
              . $call_type
              . "_exome"
              . $infile_suffix,
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

    ## VariantEval
    say {$FILEHANDLE} "## GATK varianteval";

    gatk_varianteval(
        {
            memory_allocation => q{Xmx2g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $$temp_directory_ref,
            java_jar       => catfile(
                $active_parameter_href->{gatk_path},
                "GenomeAnalysisTK.jar"
            ),
            infile_paths_ref =>
              [ $outfile_path_prefix . $call_type . "_exome" . $infile_suffix ],
            outfile_path => $outfile_path_prefix
              . $call_type
              . "_exome"
              . $outfile_suffix,
            logging_level => $active_parameter_href->{gatk_logging_level},
            referencefile_path =>
              $active_parameter_href->{human_genome_reference},
            dbsnp_file_path => $active_parameter_href->{gatk_varianteval_dbsnp},
            indel_gold_standard_file_path =>
              $active_parameter_href->{gatk_varianteval_gold},
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            infile_path => $outfile_path_prefix
              . $call_type
              . q{_exome}
              . $outfile_suffix,
            outfile_path => $outsample_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, "\n";

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        ## Collect QC metadata info for later use
        my $qc_exome_outfile =
          $outfile_prefix . $call_type . q{_exome} . $outfile_suffix;
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                sample_id        => $$sample_id_ref,
                program_name     => 'variantevalexome',
                infile           => $merged_infile_prefix,
                outdirectory     => $outsample_directory,
                outfile          => $qc_exome_outfile,
            }
        );
    }
    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        slurm_submit_job_sample_id_dependency_family_dead_end(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $$family_id_ref,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
}

sub gatk_variantevalall {

##gatk_variantevalall

##Function : GATK varianteval for all variants.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $sample_id_ref, $program_name,
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $sample_id_ref              => Sample id {REF}
##         : $program_name               => Program name

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
    my $call_type;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $sample_id_ref;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        sample_id_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$sample_id_ref
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use MIP::IO::Files qw(migrate_file);
    use MIP::Language::Java qw{java_core};
    use MIP::Set::File qw{set_file_suffix};
    use MIP::Get::File qw{get_file_suffix get_merged_infile_prefix };
    use MIP::Program::Variantcalling::Gatk qw(gatk_varianteval);
    use MIP::QC::Record qw(add_program_outfile_to_sample_info);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_family_dead_end);

    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$sample_id_ref,
            program_name          => $program_name,
            program_directory =>
              catfile( lc($$outaligner_dir_ref), lc($program_name) ),
            call_type   => $call_type,
            core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
            temp_directory => $$temp_directory_ref,
        }
    );

    ## Assign directories
    my $insample_directory = catdir( $active_parameter_href->{outdata_dir},
        $$sample_id_ref, $$outaligner_dir_ref );
    my $outsample_directory = catdir( $active_parameter_href->{outdata_dir},
        $$sample_id_ref, $$outaligner_dir_ref, lc($program_name) );
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $$sample_id_ref,
        }
    );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$$family_id_ref}{pgatk_combinevariantcallsets}
      {file_tag};
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{pgatk_combinevariantcallsets}
      {file_tag};
    my $infile_prefix       = $$family_id_ref . $infile_tag . $call_type;
    my $file_path_prefix    = catfile( $$temp_directory_ref, $infile_prefix );
    my $outfile_prefix      = $merged_infile_prefix . $outfile_tag;
    my $outfile_path_prefix = catfile( $$temp_directory_ref, $outfile_prefix );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            jobid_chain =>
              $parameter_href->{pgatk_combinevariantcallsets}{chain},
        }
    );
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_eval_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => catfile(
                $infamily_directory, $infile_prefix . $infile_suffix . q{*}
            ),
            outfile_path => $$temp_directory_ref
        }
    );
    say {$FILEHANDLE} q{wait}, "\n";

    ## Select sample id from family id vcf file

    ## GATK SelectVariants
    say {$FILEHANDLE} "## GATK SelectVariants";

    gatk_selectvariants(
        {
            memory_allocation => q{Xmx2g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $$temp_directory_ref,
            java_jar       => catfile(
                $active_parameter_href->{gatk_path},
                "GenomeAnalysisTK.jar"
            ),
            sample_names_ref => [$$sample_id_ref],
            logging_level    => $active_parameter_href->{gatk_logging_level},
            referencefile_path =>
              $active_parameter_href->{human_genome_reference},
            infile_path  => $file_path_prefix . $infile_suffix,
            outfile_path => $outfile_path_prefix . $call_type . $infile_suffix,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

    ## GATK varianteval
    say {$FILEHANDLE} "## GATK varianteval";

    gatk_varianteval(
        {
            memory_allocation => q{Xmx2g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $$temp_directory_ref,
            java_jar       => catfile(
                $active_parameter_href->{gatk_path},
                "GenomeAnalysisTK.jar"
            ),
            infile_paths_ref => [
                catfile(
                    $$temp_directory_ref,
                    $merged_infile_prefix
                      . $infile_tag
                      . $call_type
                      . $infile_suffix
                )
            ],
            outfile_path => $outfile_path_prefix . $call_type . $outfile_suffix,
            logging_level => $active_parameter_href->{gatk_logging_level},
            referencefile_path =>
              $active_parameter_href->{human_genome_reference},
            dbsnp_file_path => $active_parameter_href->{gatk_varianteval_dbsnp},
            indel_gold_standard_file_path =>
              $active_parameter_href->{gatk_varianteval_gold},
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            infile_path  => $outfile_path_prefix . $call_type . $outfile_suffix,
            outfile_path => $outsample_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, "\n";

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                sample_id        => $$sample_id_ref,
                program_name     => 'variantevalall',
                infile           => $merged_infile_prefix,
                outdirectory     => $outsample_directory,
                outfile => $outfile_prefix . $call_type . $outfile_suffix,
            }
        );
    }
    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        slurm_submit_job_sample_id_dependency_family_dead_end(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $$family_id_ref,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
}

sub snpeff {

##snpeff

##Function : snpeff annotates variants from different sources.
##Returns  : "|$xargs_file_counter"
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, $program_info_path, $file_path, $FILEHANDLE, family_id_ref, $temp_directory_ref, $outaligner_dir_ref, $call_type, $xargs_file_counter
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $program_info_path          => The program info path
##         : $file_path                  => File path
##         : $FILEHANDLE                 => Sbatch filehandle to write to
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
    my $call_type;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;
    my $program_info_path;
    my $file_path;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        file_path         => { strict_type => 1, store => \$file_path },
        FILEHANDLE    => { store => \$FILEHANDLE, },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw(get_core_number);
    use MIP::IO::Files qw(migrate_file xargs_migrate_contig_files);
    use MIP::Set::File qw{set_file_suffix};
    use MIP::Get::File qw{get_file_suffix};
    use MIP::Recipes::Xargs qw{ xargs_command };
    use Program::Variantcalling::Snpeff qw(ann);
    use Program::Variantcalling::Snpsift qw(annotate dbnsfp);
    use Program::Variantcalling::Mip qw(vcfparser);
    use MIP::QC::Record qw(add_program_outfile_to_sample_info);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_family);

    my $reduce_io_ref = \$active_parameter_href->{reduce_io};
    my $xargs_file_path_prefix;
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Get core number depending on user supplied input exists or not and max number of cores
    my $core_number = get_core_number(
        {
            module_core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            modifier_core_number => scalar( @{ $file_info_href->{contigs} } ),
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    ## Filehandles
    my $XARGSFILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    unless ( defined($FILEHANDLE) ) {           #Run as individual sbatch script

        $FILEHANDLE = IO::Handle->new();        #Create anonymous filehandle

        use MIP::Script::Setup_script qw(setup_script);

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        ( $file_path, $program_info_path ) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                job_id_href           => $job_id_href,
                FILEHANDLE            => $FILEHANDLE,
                directory_id          => $$family_id_ref,
                program_name          => $program_name,
                program_directory     => catfile( lc($$outaligner_dir_ref) ),
                call_type             => $call_type,
                core_number           => $core_number,
                process_time =>
                  $active_parameter_href->{module_time}{ "p" . $program_name },
                temp_directory => $$temp_directory_ref
            }
        );
    }

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    $parameter_href->{ "p" . $program_name }{indirectory} =
      $outfamily_directory;    #Used downstream

    ## Assign file_tags
    my $infile_tag = $file_info_href->{$$family_id_ref}{pvcfparser}{file_tag};
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{ "p" . $program_name }{file_tag};
    my $infile_prefix       = $$family_id_ref . $infile_tag . $call_type;
    my $file_path_prefix    = catfile( $$temp_directory_ref, $infile_prefix );
    my $outfile_prefix      = $$family_id_ref . $outfile_tag . $call_type;
    my $outfile_path_prefix = catfile( $$temp_directory_ref, $outfile_prefix );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            jobid_chain    => $job_id_chain,
        }
    );
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    my $vcfparser_analysis_type = "";
    my $vcfparser_contigs_ref =
      \@{ $file_info_href->{contigs_size_ordered} };    #Set default

    for (
        my $vcfparser_outfile_counter = 0 ;
        $vcfparser_outfile_counter <
        $active_parameter_href->{vcfparser_outfile_count} ;
        $vcfparser_outfile_counter++
      )
    {

        if ( $vcfparser_outfile_counter == 1 ) {

            $vcfparser_analysis_type = ".selected";     #SelectFile variants
            $vcfparser_contigs_ref =
              \@{ $file_info_href->{sorted_select_file_contigs} }
              ;                                         #Selectfile contigs
        }

        if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

            ## Copy file(s) to temporary directory
            say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
            ($xargs_file_counter) = xargs_migrate_contig_files(
                {
                    FILEHANDLE         => $FILEHANDLE,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    contigs_ref        => $vcfparser_contigs_ref,
                    file_path          => $file_path,
                    program_info_path  => $program_info_path,
                    core_number        => $core_number,
                    xargs_file_counter => $xargs_file_counter,
                    infile             => $infile_prefix,
                    file_ending        => $vcfparser_analysis_type
                      . $infile_suffix . "*",
                    indirectory    => $infamily_directory,
                    temp_directory => $$temp_directory_ref,
                }
            );
        }

        ## SnpSift Annotation
        say {$FILEHANDLE} "## SnpSift Annotation";

        my $annotation_file_counter = 0;

        if ( $active_parameter_href->{snpeff_ann} eq 1 )
        {    #Annotate using snpeff

            ## Create file commands for xargs
            ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
                {
                    FILEHANDLE         => $FILEHANDLE,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    file_path          => $file_path,
                    program_info_path  => $program_info_path,
                    core_number        => $core_number,
                    xargs_file_counter => $xargs_file_counter,
                    first_command      => "java",
                    memory_allocation  => "Xmx4g -XX:-UseConcMarkSweepGC",
                    java_use_large_pages =>
                      $active_parameter_href->{java_use_large_pages},
                    temp_directory => $$temp_directory_ref,
                    java_jar       => catfile(
                        $active_parameter_href->{snpeff_path}, "snpEff.jar"
                    ),
                }
            );

            foreach my $contig (@$vcfparser_contigs_ref) {

                Program::Variantcalling::Snpeff::ann(
                    {
                        verbosity => "v",
                        genome_build_version =>
                          $active_parameter_href->{snpeff_genome_build_version},
                        config_file_path => catfile(
                            $active_parameter_href->{snpeff_path},
                            "snpEff.config"
                        ),
                        infile_path => $file_path_prefix . "_"
                          . $contig
                          . $vcfparser_analysis_type
                          . $infile_suffix,
                        outfile_path => $file_path_prefix . "_"
                          . $contig
                          . $vcfparser_analysis_type
                          . $infile_suffix . "."
                          . $xargs_file_counter,
                        stderrfile_path => $xargs_file_path_prefix . "."
                          . $contig
                          . ".stderr.txt",
                        FILEHANDLE => $XARGSFILEHANDLE,
                    }
                );
                say {$XARGSFILEHANDLE} "\n";
            }
            $annotation_file_counter = $xargs_file_counter;
        }

        while ( my ( $annotation_file, $annotation_info_key ) =
            each( %{ $active_parameter_href->{snpsift_annotation_files} } ) )
        {

            ## Create file commands for xargs
            ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
                {
                    FILEHANDLE         => $FILEHANDLE,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    file_path          => $file_path,
                    program_info_path  => $program_info_path,
                    core_number        => $core_number,
                    xargs_file_counter => $xargs_file_counter,
                    first_command      => "java",
                    memory_allocation  => "Xmx2g -XX:-UseConcMarkSweepGC",
                    java_use_large_pages =>
                      $active_parameter_href->{java_use_large_pages},
                    temp_directory => $$temp_directory_ref,
                    java_jar       => catfile(
                        $active_parameter_href->{snpeff_path},
                        "SnpSift.jar"
                    ),
                }
            );
            ## Get parameters
            my $name_prefix;
            my $info_key;
            if ( defined($annotation_info_key) ) {

                ## Apply specific INFO field output key for easier downstream processing
                if (
                    defined(
                        $active_parameter_href->{snpsift_annotation_outinfo_key}
                          {$annotation_file}
                    )
                  )
                {

                    $name_prefix =
                      $active_parameter_href->{snpsift_annotation_outinfo_key}
                      {$annotation_file};
                }
                $info_key = $annotation_info_key;    #Database
            }

            foreach my $contig (@$vcfparser_contigs_ref) {

                ##Get contig specific parameters
                my $infile_path;
                if ( !$annotation_file_counter ) {    #First file per contig

                    $infile_path =
                        $file_path_prefix . "_"
                      . $contig
                      . $vcfparser_analysis_type
                      . $infile_suffix;
                }
                else {

                    my $annotation_infile_number = $xargs_file_counter - 1;
                    $infile_path =
                        $file_path_prefix . "_"
                      . $contig
                      . $vcfparser_analysis_type
                      . $infile_suffix . "."
                      . $annotation_infile_number;   #Infile from previous round
                }
                Program::Variantcalling::Snpsift::annotate(
                    {
                        verbosity    => "v",
                        infile_path  => $infile_path,
                        outfile_path => $file_path_prefix . "_"
                          . $contig
                          . $vcfparser_analysis_type
                          . $infile_suffix . "."
                          . $xargs_file_counter,
                        config_file_path => catfile(
                            $active_parameter_href->{snpeff_path},
                            "snpEff.config"
                        ),
                        database_path   => $annotation_file,
                        name_prefix     => $name_prefix,
                        info            => $info_key,
                        stderrfile_path => $xargs_file_path_prefix . "."
                          . $contig
                          . ".stderr.txt",
                        append_stderr_info => 1,
                        FILEHANDLE         => $XARGSFILEHANDLE,
                    }
                );
                say {$XARGSFILEHANDLE} "\n";
            }
            $annotation_file_counter++;    #Increment counter
            close $XARGSFILEHANDLE;
        }

        if ( @{ $active_parameter_href->{snpsift_dbnsfp_annotations} } ) {

            ## SnpSiftDbNSFP Annotation
            say {$FILEHANDLE} "## SnpSiftDnNSFP Annotation";

            ## Create file commands for xargs
            ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
                {
                    FILEHANDLE         => $FILEHANDLE,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    file_path          => $file_path,
                    program_info_path  => $program_info_path,
                    core_number        => $core_number,
                    xargs_file_counter => $xargs_file_counter,
                    first_command      => "java",
                    memory_allocation  => "Xmx2g -XX:-UseConcMarkSweepGC",
                    java_use_large_pages =>
                      $active_parameter_href->{java_use_large_pages},
                    temp_directory => $$temp_directory_ref,
                    java_jar       => catfile(
                        $active_parameter_href->{snpeff_path},
                        "SnpSift.jar"
                    ),
                }
            );

            my $annotation_infile_number = $xargs_file_counter - 1;

            foreach my $contig (@$vcfparser_contigs_ref) {

                Program::Variantcalling::Snpsift::dbnsfp(
                    {
                        annotate_fields_ref => \@{
                            $active_parameter_href->{snpsift_dbnsfp_annotations}
                        },
                        infile_path => $file_path_prefix . "_"
                          . $contig
                          . $vcfparser_analysis_type
                          . $infile_suffix . "."
                          . $annotation_infile_number,
                        outfile_path => $file_path_prefix . "_"
                          . $contig
                          . $vcfparser_analysis_type
                          . $infile_suffix . "."
                          . $xargs_file_counter,
                        config_file_path => catfile(
                            $active_parameter_href->{snpeff_path},
                            "snpEff.config"
                        ),
                        database_path =>
                          $active_parameter_href->{snpsift_dbnsfp_file},
                        stderrfile_path => $xargs_file_path_prefix . "."
                          . $contig
                          . ".stderr.txt",
                        append_stderr_info => 1,
                        FILEHANDLE         => $XARGSFILEHANDLE,
                        verbosity          => "v",
                    }
                );
                say {$XARGSFILEHANDLE} "\n";
            }
            close $XARGSFILEHANDLE;
        }

        ## Add INFO headers and FIX_INFO for annotations using vcfparser
        say {$FILEHANDLE}
          "## Add INFO headers and FIX_INFO for annotations using vcfparser";

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                FILEHANDLE         => $FILEHANDLE,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                file_path          => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
                xargs_file_counter => $xargs_file_counter,
            }
        );

        my $annotation_infile_number = $xargs_file_counter - 1;

        foreach my $contig (@$vcfparser_contigs_ref) {

            vcfparser(
                {
                    infile_path => $file_path_prefix . "_"
                      . $contig
                      . $vcfparser_analysis_type
                      . $infile_suffix . "."
                      . $annotation_infile_number,
                    outfile_path => $outfile_path_prefix . "_"
                      . $contig
                      . $vcfparser_analysis_type
                      . $outfile_suffix,
                    stderrfile_path => $xargs_file_path_prefix . "."
                      . $contig
                      . ".stderr.txt ",
                    append_stderr_info => 1,
                    FILEHANDLE         => $XARGSFILEHANDLE,
                }
            );
            say {$XARGSFILEHANDLE} "\n";
        }

        if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

            ## Copies file from temporary directory. Per contig
            say {$FILEHANDLE} q{## Copy file from temporary directory};
            ($xargs_file_counter) = xargs_migrate_contig_files(
                {
                    FILEHANDLE        => $FILEHANDLE,
                    XARGSFILEHANDLE   => $XARGSFILEHANDLE,
                    contigs_ref       => $vcfparser_contigs_ref,
                    file_path         => $file_path,
                    program_info_path => $program_info_path,
                    core_number => $active_parameter_href->{max_cores_per_node},
                    xargs_file_counter => $xargs_file_counter,
                    outfile            => $outfile_prefix,
                    file_ending        => $vcfparser_analysis_type
                      . $outfile_suffix . "*",
                    outdirectory   => $outfamily_directory,
                    temp_directory => $$temp_directory_ref,
                }
            );
        }
        else {

            ## QC Data File(s)
            migrate_file(
                {
                    infile_path => $outfile_path_prefix . "_"
                      . $file_info_href->{contigs_size_ordered}[0]
                      . $vcfparser_analysis_type
                      . $outfile_suffix,
                    outfile_path => $outfamily_directory,
                    FILEHANDLE   => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} q{wait}, "\n";
        }
    }

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        ## Collect QC metadata info for later use
        my $qc_snpeff_outfile =
            $outfile_prefix . q{_}
          . $file_info_href->{contigs_size_ordered}[0]
          . $vcfparser_analysis_type
          . $outfile_suffix;
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => $program_name,
                outdirectory     => $outfamily_directory,
                outfile          => $qc_snpeff_outfile,
            }
        );
    }

    if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

        close $FILEHANDLE;

        if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

            slurm_submit_job_sample_id_dependency_add_to_family(
                {
                    job_id_href             => $job_id_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    sample_ids_ref =>
                      \@{ $active_parameter_href->{sample_ids} },
                    family_id        => $$family_id_ref,
                    path             => $job_id_chain,
                    log              => $log,
                    sbatch_file_name => $file_path,
                }
            );
        }
    }
    if ($$reduce_io_ref) {

        return
          $xargs_file_counter
          ; #Track the number of created xargs scripts per module for Block algorithm
    }
}

sub mvcfparser {

##mvcfparser

##Function : Vcfparser performs parsing of varianteffectpredictor annotated variants
##Returns  : "|$xargs_file_counter"
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, $program_info_path, $file_path, $FILEHANDLE, family_id_ref, $temp_directory_ref, $outaligner_dir_ref, $call_type, $xargs_file_counter
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $program_info_path          => The program info path
##         : $file_path                  => File path
##         : $FILEHANDLE                 => Sbatch filehandle to write to
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type
##         : $xargs_file_counter         => The xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
    my $call_type;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_info_path;
    my $program_name;
    my $file_path;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        file_path         => { strict_type => 1, store => \$file_path },
        FILEHANDLE    => { store => \$FILEHANDLE, },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw(get_core_number);
    use MIP::IO::Files qw(migrate_file xargs_migrate_contig_files);
    use MIP::Set::File qw{set_file_suffix};
    use MIP::Get::File qw{get_file_suffix};
    use MIP::Recipes::Xargs qw{ xargs_command };
    use Program::Variantcalling::Mip qw(vcfparser);
    use MIP::QC::Record qw(add_program_outfile_to_sample_info);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_family);

    my $reduce_io_ref = \$active_parameter_href->{reduce_io};
    my $xargs_file_path_prefix;
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Get core number depending on user supplied input exists or not and max number of cores
    my $core_number = get_core_number(
        {
            module_core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            modifier_core_number => scalar( @{ $file_info_href->{contigs} } ),
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    ## Filehandles
    my $XARGSFILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    unless ( defined($FILEHANDLE) ) {           #Run as individual sbatch script

        $FILEHANDLE = IO::Handle->new();        #Create anonymous filehandle

        use MIP::Script::Setup_script qw(setup_script);

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        ( $file_path, $program_info_path ) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                job_id_href           => $job_id_href,
                FILEHANDLE            => $FILEHANDLE,
                directory_id          => $$family_id_ref,
                program_name          => $program_name,
                program_directory     => catfile( lc($$outaligner_dir_ref) ),
                call_type             => $call_type,
                temp_directory        => $$temp_directory_ref,
                core_number           => $core_number,
                process_time =>
                  $active_parameter_href->{module_time}{ "p" . $program_name },
            }
        );
    }

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    $parameter_href->{ "p" . $program_name }{indirectory} =
      $outfamily_directory;    #Used downstream

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$$family_id_ref}{pvarianteffectpredictor}{file_tag};
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{ "p" . $program_name }{file_tag};
    my $infile_prefix       = $$family_id_ref . $infile_tag . $call_type;
    my $file_path_prefix    = catfile( $$temp_directory_ref, $infile_prefix );
    my $outfile_prefix      = $$family_id_ref . $outfile_tag . $call_type;
    my $outfile_path_prefix = catfile( $$temp_directory_ref, $outfile_prefix );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            jobid_chain    => $job_id_chain,
        }
    );
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        ($xargs_file_counter) = xargs_migrate_contig_files(
            {
                FILEHANDLE      => $FILEHANDLE,
                XARGSFILEHANDLE => $XARGSFILEHANDLE,
                contigs_ref => \@{ $file_info_href->{contigs_size_ordered} },
                file_path   => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
                xargs_file_counter => $xargs_file_counter,
                infile             => $infile_prefix,
                indirectory        => $infamily_directory,
                temp_directory     => $$temp_directory_ref,
            }
        );
    }

    ## vcfparser
    say {$FILEHANDLE} "## vcfparser";

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## Get parameters
        my $padding;
        if ( $contig =~ /MT|M/ ) {

            $padding = 10;    #Special case for mitochondrial contig annotation
        }

        my @select_feature_annotation_columns;
        my $select_file;
        my $select_file_matching_column;
        my $select_outfile;
        if ( $active_parameter_href->{vcfparser_select_file} ) {

            if (
                !check_entry_hash_of_array(
                    {
                        hash_ref => $file_info_href,
                        key      => "select_file_contigs",
                        element  => $contig,
                    }
                )
              )
            {

                $select_file =
                  catfile( $active_parameter_href->{vcfparser_select_file} )
                  ;    #List of genes to analyse separately
                $select_file_matching_column = $active_parameter_href
                  ->{vcfparser_select_file_matching_column}
                  ;    #Column of HGNC Symbol in SelectFile (-sf)

                if (
                    (
                        $active_parameter_href
                        ->{vcfparser_select_feature_annotation_columns}
                    )
                    && (
                        @{
                            $active_parameter_href
                              ->{vcfparser_select_feature_annotation_columns}
                        }
                    )
                  )
                {

                    @select_feature_annotation_columns =
                      @{ $active_parameter_href
                          ->{vcfparser_select_feature_annotation_columns} };
                }
                $select_outfile =
                    $outfile_path_prefix . "_"
                  . $contig
                  . ".selected"
                  . $infile_suffix;
            }
        }

        vcfparser(
            {
                range_feature_annotation_columns_ref => \@{
                    $active_parameter_href
                      ->{vcfparser_range_feature_annotation_columns}
                },
                select_feature_annotation_columns_ref =>
                  \@select_feature_annotation_columns,
                infile_path => $file_path_prefix . "_"
                  . $contig
                  . $infile_suffix,
                outfile_path => $outfile_path_prefix . "_"
                  . $contig
                  . $infile_suffix,
                stderrfile_path => $xargs_file_path_prefix . "."
                  . $contig
                  . ".stderr.txt ",
                range_feature_file_path =>
                  $active_parameter_href->{vcfparser_range_feature_file},
                select_feature_file_path       => $select_file,
                select_feature_matching_column => $select_file_matching_column,
                select_outfile                 => $select_outfile,
                parse_vep  => $active_parameter_href->{pvarianteffectpredictor},
                padding    => $padding,
                FILEHANDLE => $XARGSFILEHANDLE,
            }
        );
        say {$XARGSFILEHANDLE} "\n";
    }

    ## QC Data File(s)
    migrate_file(
        {
            infile_path => $outfile_path_prefix . q{_}
              . $file_info_href->{contigs_size_ordered}[0]
              . $infile_suffix,
            outfile_path => $outfamily_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, "\n";

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        ## Clear old vcfparser entry if present
        if ( defined( $sample_info_href->{$program_name} ) ) {

            delete( $sample_info_href->{$program_name} );
        }

        my %gene_panels = (
            range_file  => "vcfparser_range_feature_file",
            select_file => "vcfparser_select_file",
        );
        while ( my ( $gene_panel_key, $gene_panel_file ) = each(%gene_panels) )
        {

            ## Collect databases(s) from a potentially merged gene panel file and adds them to sample_info
            collect_gene_panels(
                {
                    sample_info_href => $sample_info_href,
                    family_id_ref    => $family_id_ref,
                    program_name_ref => \$program_name,
                    aggregate_gene_panel_file =>
                      $active_parameter_href->{$gene_panel_file},
                    aggregate_gene_panels_key => $gene_panel_key,
                }
            );
        }

        ## Collect QC metadata info for later use
        my $qc_vcfparser_outfile =
            $outfile_prefix . q{_}
          . $file_info_href->{contigs_size_ordered}[0]
          . $infile_suffix;
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => $program_name,
                outdirectory     => $outfamily_directory,
                outfile          => $qc_vcfparser_outfile,
            }
        );
    }

    close $XARGSFILEHANDLE;

    if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

        my $vcfparser_analysis_type = "";
        my @vcfparser_contigs_ref =
          \@{ $file_info_href->{contigs_size_ordered} };

        for (
            my $vcfparser_outfile_counter = 0 ;
            $vcfparser_outfile_counter <
            $active_parameter_href->{vcfparser_outfile_count} ;
            $vcfparser_outfile_counter++
          )
        {

            if ( $vcfparser_outfile_counter == 1 ) {

                $vcfparser_analysis_type = ".selected";    #SelectFile variants
                @vcfparser_contigs_ref =
                  \@{ $file_info_href->{sorted_select_file_contigs} };
            }

            ## Copies file from temporary directory.
            say {$FILEHANDLE} "## Copy file(s) from temporary directory";
            ($xargs_file_counter) = xargs_migrate_contig_files(
                {
                    FILEHANDLE         => $FILEHANDLE,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    contigs_ref        => @vcfparser_contigs_ref,
                    file_path          => $file_path,
                    program_info_path  => $program_info_path,
                    core_number        => $core_number,
                    xargs_file_counter => $xargs_file_counter,
                    outfile            => $outfile_prefix,
                    file_ending        => $vcfparser_analysis_type
                      . $infile_suffix . "*",
                    outdirectory   => $outfamily_directory,
                    temp_directory => $$temp_directory_ref,
                }
            );
        }
        close $FILEHANDLE;
    }

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

            slurm_submit_job_sample_id_dependency_add_to_family(
                {
                    job_id_href             => $job_id_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    sample_ids_ref =>
                      \@{ $active_parameter_href->{sample_ids} },
                    family_id        => $$family_id_ref,
                    path             => $job_id_chain,
                    log              => $log,
                    sbatch_file_name => $file_path,
                }
            );
        }
    }
    if ($$reduce_io_ref) {

        return
          $xargs_file_counter
          ; #Track the number of created xargs scripts per module for Block algorithm
    }
}

sub mpeddy {

##mpeddy

##Function : Compares familial-relationships and sexes.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, family_id_ref, $temp_directory_ref, $outaligner_dir_ref, $call_type
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
    my $call_type;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{get_file_suffix};
    use MIP::Script::Setup_script qw(setup_script);
    use MIP::IO::Files qw(migrate_file);
    use Program::Variantcalling::Peddy qw(peddy);
    use MIP::QC::Record qw(add_program_metafile_to_sample_info);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_family_dead_end);

    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory     => catfile(
                lc($$outaligner_dir_ref),
                "casecheck", lc($program_name)
            ),
            call_type   => $call_type,
            core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
            source_environment_commands_ref => [],
        }
    );

    my ( $volume, $directory, $program_info_file ) =
      File::Spec->splitpath($program_info_path)
      ;    #Split to enable submission to &sample_info_qc later
    my $stderr_file = $program_info_file
      . ".stderr.txt";    #To enable submission to &sample_info_qc later
    my $stdout_file = $program_info_file
      . ".stdout.txt";    #To enable submission to &sample_info_qc later

    ## Assign Directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    my $outfamily_directory = catfile( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref, "casecheck", lc($program_name) );
    my $outfamily_file_directory =
      catfile( $active_parameter_href->{outdata_dir}, $$family_id_ref );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$$family_id_ref}{pgatk_combinevariantcallsets}
      {file_tag};
    my $infile_prefix       = $$family_id_ref . $infile_tag . $call_type;
    my $file_path_prefix    = catfile( $$temp_directory_ref, $infile_prefix );
    my $outfile_path_prefix = catfile( $outfamily_directory, $$family_id_ref );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain_vcf_data
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            jobid_chain =>
              $parameter_href->{pgatk_combinevariantcallsets}{chain},
        }
    );

    my $suffix = ".vcf.gz";

    my $family_file =
      catfile( $outfamily_file_directory, $$family_id_ref . ".fam" );
    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            parameter_href        => $parameter_href,
            active_parameter_href => $active_parameter_href,
            sample_info_href      => $sample_info_href,
            FILEHANDLE            => $FILEHANDLE,
            fam_file_path         => $family_file,
        }
    );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => catfile(
                $infamily_directory, $infile_prefix . $infile_suffix . q{*}
            ),
            outfile_path => $$temp_directory_ref
        }
    );
    say {$FILEHANDLE} q{wait}, "\n";

    ## Reformat variant calling file and index
    view_vcf(
        {
            infile_path         => $file_path_prefix . $infile_suffix,
            outfile_path_prefix => $file_path_prefix,
            output_type         => "z",
            index               => 1,
            index_type          => "tbi",
            FILEHANDLE          => $FILEHANDLE,
        }
    );

    ## peddy
    peddy(
        {
            infile_path         => $file_path_prefix . $suffix,
            outfile_prefix_path => $outfile_path_prefix,
            family_file_path    => $family_file,
            FILEHANDLE          => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        my %peddy_output = (
            ped_check => 'csv',
            sex_check => 'csv',
            peddy     => 'ped',
        );

      PEDDY_OUTPUT_FILES:
        while ( my ( $file_key, $suffix ) = each %peddy_output ) {

            my $outfile_suffix = '.' . $file_key . '.' . $suffix;
            ## Collect QC metadata info for later use
            add_program_metafile_to_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    program_name     => $program_name,
                    metafile_tag     => $file_key,
                    path             => $outfile_path_prefix . $outfile_suffix,
                }
            );
        }
    }

    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        slurm_submit_job_sample_id_dependency_family_dead_end(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $$family_id_ref,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
}

sub mplink {

##mplink

##Function : Tests sample for correct relatives (only performed for samples with relatives defined in pedigree file) performed on sequence data.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, family_id_ref, $temp_directory_ref, $outaligner_dir_ref, $call_type
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
    my $call_type;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use MIP::Get::File qw{get_file_suffix};
    use MIP::Program::Variantcalling::Bcftools
      qw(bcftools_view bcftools_annotate);
    use Program::Variantcalling::Vt qw(vt_uniq);
    use Program::Variantcalling::Plink qw(plink);
    use MIP::QC::Record qw(add_program_outfile_to_sample_info);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_family_dead_end);

    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory     => catfile(
                lc($$outaligner_dir_ref),
                "casecheck", lc($program_name)
            ),
            call_type   => $call_type,
            core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
        }
    );

    my ( $volume, $directory, $program_info_file ) =
      File::Spec->splitpath($program_info_path)
      ;    #Split to enable submission to &sample_info_qc later
    my $stderr_file = $program_info_file
      . ".stderr.txt";    #To enable submission to &sample_info_qc later
    my $stdout_file = $program_info_file
      . ".stdout.txt";    #To enable submission to &sample_info_qc later

    ## Assign Directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    my $outfamily_directory = catfile( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref, "casecheck", lc($program_name) );
    my $outfamily_file_directory =
      catfile( $active_parameter_href->{outdata_dir}, $$family_id_ref );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$$family_id_ref}{pgatk_combinevariantcallsets}
      {file_tag};
    my $infile_prefix = $$family_id_ref . $infile_tag . $call_type;
    my $file_path_prefix = catfile( $$temp_directory_ref, $infile_prefix );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain_vcf_data
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            jobid_chain =>
              $parameter_href->{pgatk_combinevariantcallsets}{chain},
        }
    );

    my $family_file =
      catfile( $outfamily_file_directory, $$family_id_ref . ".fam" );
    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            parameter_href        => $parameter_href,
            active_parameter_href => $active_parameter_href,
            sample_info_href      => $sample_info_href,
            FILEHANDLE            => $FILEHANDLE,
            fam_file_path         => $family_file,
        }
    );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => catfile(
                $infamily_directory, $infile_prefix . $infile_suffix . q{*}
            ),
            outfile_path => $$temp_directory_ref
        }
    );
    say {$FILEHANDLE} q{wait}, "\n";

    ## Prepare input
    say {$FILEHANDLE} "## Remove indels using bcftools ";
    bcftools_view(
        {
            infile_path  => $file_path_prefix . $infile_suffix,
            outfile_path => $file_path_prefix . "_no_indels" . $infile_suffix,
            output_type  => "v",
            exclude_types_ref => ["indels"],
            FILEHANDLE        => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

    say {$FILEHANDLE} "## Create uniq IDs and remove duplicate variants";
    bcftools_annotate(
        {
            remove_ids_ref => ["ID"],
            set_id         => q?+'%CHROM:%POS:%REF:%ALT'?,
            infile_path    => $file_path_prefix . "_no_indels" . $infile_suffix,
            output_type    => "v",
            FILEHANDLE     => $FILEHANDLE,
        }
    );

    print {$FILEHANDLE} "| ";    #Pipe

    ## Drops duplicate variants that appear later in the the VCF file
    Program::Variantcalling::Vt::vt_uniq(
        {
            infile_path  => "-",
            outfile_path => $file_path_prefix
              . "_no_indels_ann_uniq"
              . $infile_suffix,
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

    ### Plink
    say {$FILEHANDLE} "## Create pruning set and uniq IDs";
    plink(
        {
            vcffile_path => $file_path_prefix
              . "_no_indels_ann_uniq"
              . $infile_suffix,
            outfile_prefix =>
              catfile( $$temp_directory_ref, $$family_id_ref . "_data" ),
            vcf_require_gt      => 1,
            vcf_half_call       => "haploid",
            set_missing_var_ids => q?@:#[?
              . $file_info_href->{human_genome_reference_version}
              . q?]\$1,\$2?,
            const_fid           => $$family_id_ref,
            make_bed            => 1,
            indep               => 1,
            indep_window_size   => 50,
            indep_step_size     => 5,
            indep_vif_threshold => 2,
            FILEHANDLE          => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

    say {$FILEHANDLE}
      "## Update Plink fam. Create ped and map file and frequency report";
    plink(
        {
            binary_fileset_prefix =>
              catfile( $$temp_directory_ref, $$family_id_ref . "_data" ),
            fam_file_path => $family_file,
            make_just_fam => 1,
            recode        => 1,
            freqx         => 1,
            outfile_prefix =>
              catfile( $$temp_directory_ref, $$family_id_ref . "_data" ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

    # Only perform if more than 1 sample
    if ( scalar( @{ $active_parameter_href->{sample_ids} } ) > 1 ) {

        say {$FILEHANDLE} "## Calculate inbreeding coefficients per family";
        plink(
            {
                binary_fileset_prefix =>
                  catfile( $$temp_directory_ref, $$family_id_ref . "_data" ),
                outfile_prefix =>
                  catfile( $outfamily_directory, $$family_id_ref ),
                het                     => 1,
                small_sample            => 1,
                inbreeding_coefficients => 1,
                extract_file            => catfile(
                    $$temp_directory_ref, $$family_id_ref . "_data.prune.in"
                ),
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} "\n";

        say {$FILEHANDLE} "## Create Plink .mibs per family";
        plink(
            {
                ped_file_path => catfile(
                    $$temp_directory_ref, $$family_id_ref . "_data.ped"
                ),
                map_file_path => catfile(
                    $$temp_directory_ref, $$family_id_ref . "_data.map"
                ),
                cluster => 1,
                matrix  => 1,
                outfile_prefix =>
                  catfile( $outfamily_directory, $$family_id_ref ),
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} "\n";
    }

    ### Plink sex-check
    ## Get parameters
    my $genome_build;
    if ( $file_info_href->{human_genome_reference_source} eq "GRCh" ) {

        $genome_build = "b" . $file_info_href->{human_genome_reference_version};
    }
    else {

        $genome_build =
          "hg" . $file_info_href->{human_genome_reference_version};
    }
    plink(
        {
            regions_ref => [ 23, 24 ],
            split_x     => $genome_build,
            no_fail     => 1,
            make_bed    => 1,
            binary_fileset_prefix =>
              catfile( $$temp_directory_ref, $$family_id_ref . "_data" ),
            outfile_prefix => catfile(
                $$temp_directory_ref, $$family_id_ref . "_data_unsplit"
            ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

    ##Get parameters
    my $sex_check_min_F;
    if (   ( $consensus_analysis_type eq "wes" )
        || ( $consensus_analysis_type eq "rapid" ) )
    {    #Exome/rapid analysis use combined reference for more power

        $sex_check_min_F = "0.2 0.75";
    }
    my $extract_file;
    if ( scalar( @{ $active_parameter_href->{sample_ids} } ) > 1 ) {

        $extract_file =
          catfile( $$temp_directory_ref, $$family_id_ref . '_data.prune.in' );
    }
    plink(
        {
            check_sex       => 1,
            sex_check_min_F => $sex_check_min_F,
            extract_file    => $extract_file,
            read_freqfile_path =>
              catfile( $$temp_directory_ref, $$family_id_ref . "_data.frqx" ),
            binary_fileset_prefix => catfile(
                $$temp_directory_ref, $$family_id_ref . "_data_unsplit"
            ),
            outfile_prefix => catfile( $outfamily_directory, $$family_id_ref ),
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        if ( scalar( @{ $active_parameter_href->{sample_ids} } ) > 1 )
        {    #Only perform if more than 1 sample

            ## Collect QC metadata info for later use
            add_program_outfile_to_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    program_name     => q{inbreeding_factor},
                    outdirectory     => $outfamily_directory,
                    outfile          => $$family_id_ref . q{.het},
                }
            );

            ## Collect QC metadata info for later use
            add_program_outfile_to_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    program_name     => q{relation_check},
                    outdirectory     => $outfamily_directory,
                    outfile          => $$family_id_ref . q{.mibs},
                }
            );
        }

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => q{plink_sexcheck},
                outdirectory     => $outfamily_directory,
                outfile          => $$family_id_ref . q{.sexcheck},
            }
        );

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => q{plink2},
                outdirectory     => $directory,
                outfile          => $stdout_file,
            }
        );
    }

    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        slurm_submit_job_sample_id_dependency_family_dead_end(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $$family_id_ref,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
}

sub vt {

##vt

##Function : Split multi allelic records into single records and normalize
##Returns  : "|$xargs_file_counter"
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, $program_info_path, $file_path, $stderr_path, $FILEHANDLE, family_id_ref, $temp_directory_ref, $outaligner_dir_ref, $call_type, $xargs_file_counter
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $program_info_path          => The program info path
##         : $file_path                  => File path
##         : $stderr_path                => The stderr path of the block script
##         : $FILEHANDLE                 => Filehandle to write to
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type
##         : $xargs_file_counter         => The xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
    my $call_type;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;
    my $program_info_path;
    my $file_path;
    my $stderr_path;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        file_path         => { strict_type => 1, store => \$file_path },
        stderr_path       => { strict_type => 1, store => \$stderr_path },
        FILEHANDLE    => { store => \$FILEHANDLE, },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Readonly;
    use MIP::Cluster qw(get_core_number);
    use MIP::Set::File qw{set_file_suffix};
    use MIP::Get::File qw{get_file_suffix};
    use MIP::Recipes::Xargs qw{ xargs_command };
    use MIP::Gnu::Coreutils qw(gnu_mv);
    use Program::Variantcalling::Genmod qw(annotate filter);
    use MIP::QC::Record qw(add_program_outfile_to_sample_info);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_family);

    ## Constants
    Readonly my $DOT => q{.};

    my $reduce_io_ref = \$active_parameter_href->{reduce_io};
    my $xargs_file_path_prefix;
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $XARGSFILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    unless ( defined($FILEHANDLE) ) {           #Run as individual sbatch script

        $FILEHANDLE = IO::Handle->new();        #Create anonymous filehandle
    }

    ## Get core number depending on user supplied input exists or not and max number of cores
    my $core_number = get_core_number(
        {
            module_core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            modifier_core_number => scalar( @{ $file_info_href->{contigs} } ),
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

        use MIP::Script::Setup_script qw(setup_script);

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        ( $file_path, $program_info_path ) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                job_id_href           => $job_id_href,
                FILEHANDLE            => $FILEHANDLE,
                directory_id          => $$family_id_ref,
                program_name          => $program_name,
                program_directory     => catfile( lc($$outaligner_dir_ref) ),
                call_type             => $call_type,
                core_number           => $core_number,
                process_time =>
                  $active_parameter_href->{module_time}{ "p" . $program_name },
                temp_directory => $$temp_directory_ref,
            }
        );
        $stderr_path = $program_info_path . ".stderr.txt";
    }
    my ( $volume, $directory, $stderr_file ) =
      File::Spec->splitpath($stderr_path)
      ;    #Split to enable submission to &sample_info_qc later

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    $parameter_href->{ "p" . $program_name }{indirectory} =
      $outfamily_directory;    #Used downstream

    ## Assign file_tags
    my $infile_tag = $file_info_href->{$$family_id_ref}{prhocall}{file_tag};
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{ "p" . $program_name }{file_tag};
    my $infile_prefix       = $$family_id_ref . $infile_tag . $call_type;
    my $file_path_prefix    = catfile( $$temp_directory_ref, $infile_prefix );
    my $outfile_prefix      = $$family_id_ref . $outfile_tag . $call_type;
    my $outfile_path_prefix = catfile( $$temp_directory_ref, $outfile_prefix );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            jobid_chain    => $job_id_chain,
        }
    );
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        ($xargs_file_counter) = xargs_migrate_contig_files(
            {
                FILEHANDLE      => $FILEHANDLE,
                XARGSFILEHANDLE => $XARGSFILEHANDLE,
                contigs_ref => \@{ $file_info_href->{contigs_size_ordered} },
                file_path   => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
                xargs_file_counter => $xargs_file_counter,
                infile             => $infile_prefix,
                indirectory        => $infamily_directory,
                temp_directory     => $$temp_directory_ref,
            }
        );
    }

    say {$FILEHANDLE}
"## vt - Decompose (split multi allelic records into single records) and/or normalize variants";

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    my $remove_star_regexp =
      q?perl -nae \'unless\($F\[4\] eq \"\*\") \{print $_\}\' ?
      ; #VEP does not annotate '*' since the alt allele does not exist, this is captured in the upstream indel and SNV record associated with '*'

    ## Split vcf into contigs
    while ( my ( $contig_index, $contig ) =
        each( @{ $file_info_href->{contigs_size_ordered} } ) )
    {

        ## vt - Split multi allelic records into single records and normalize
        analysis_vt_core_rio(
            {
                active_parameter_href => $active_parameter_href,
                FILEHANDLE            => $XARGSFILEHANDLE,
                infile_path           => $file_path_prefix . "_"
                  . $contig
                  . $infile_suffix,
                outfile_path => $outfile_path_prefix . "_"
                  . $contig
                  . $outfile_suffix,
                decompose => $active_parameter_href->{vt_decompose},
                normalize => $active_parameter_href->{vt_normalize},
                uniq      => $active_parameter_href->{vt_uniq},
                gnu_sed   => 1,
                instream  => 0,
                cmd_break => q{;},
                xargs_file_path_prefix => $xargs_file_path_prefix,
                contig                 => $contig,
            }
        );

        if (   ( $contig_index == 0 )
            && ( $active_parameter_href->{ "p" . $program_name } == 1 ) )
        {

            my ( $volume, $directory, $stderr_file ) =
              File::Spec->splitpath($xargs_file_path_prefix)
              ;    #Split to enable submission to &SampleInfoQC later

            ## Collect QC metadata info for later use
            my $qc_vt_outfile =
              $stderr_file . $DOT . $contig . $DOT . q{stderr.txt};
            add_program_outfile_to_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    program_name     => 'vt',
                    outdirectory     => $directory,
                    outfile          => $qc_vt_outfile,
                }
            );
        }

        my $alt_file_tag = "";

        ## Remove decomposed '*' entries
        if ( $active_parameter_href->{vt_missing_alt_allele} ) {

            $alt_file_tag = "_nostar";
            print $XARGSFILEHANDLE catfile(
                $remove_star_regexp . $$temp_directory_ref,
                $outfile_prefix . "_" . $contig . $outfile_suffix
            ) . " ";
            print $XARGSFILEHANDLE "> "
              . $outfile_path_prefix . "_"
              . $contig
              . $alt_file_tag
              . $outfile_suffix . " ";
            print $XARGSFILEHANDLE "2>> "
              . $xargs_file_path_prefix . "."
              . $contig
              . ".stderr.txt "
              ;    #Redirect xargs output to program specific stderr file
            print $XARGSFILEHANDLE "; ";
        }

        ## Remove common variants
        if ( $active_parameter_href->{vt_genmod_filter} ) {

            Program::Variantcalling::Genmod::annotate(
                {
                    infile_path => $outfile_path_prefix . "_"
                      . $contig
                      . $alt_file_tag
                      . $outfile_suffix,
                    outfile_path => catfile( dirname( devnull() ), "stdout" ),
                    stderrfile_path => $xargs_file_path_prefix . "."
                      . $contig
                      . ".stderr.txt",
                    append_stderr_info  => 1,
                    verbosity           => "v",
                    temp_directory_path => $$temp_directory_ref,
                    thousand_g_file_path =>
                      $active_parameter_href->{vt_genmod_filter_1000g},
                    max_af => $active_parameter_href->{vt_genmod_filter_max_af},
                    FILEHANDLE => $XARGSFILEHANDLE,
                }
            );
            print $XARGSFILEHANDLE "| ";

            $alt_file_tag .= "_genmod_filter";    #Update file tag

            Program::Variantcalling::Genmod::filter(
                {
                    infile_path  => "-",
                    outfile_path => $outfile_path_prefix . "_"
                      . $contig
                      . $alt_file_tag
                      . $outfile_suffix,
                    stderrfile_path => $xargs_file_path_prefix . "."
                      . $contig
                      . ".stderr.txt",
                    append_stderr_info => 1,
                    verbosity          => "v",
                    threshold =>
                      $active_parameter_href->{sv_genmod_filter_threshold},
                    FILEHANDLE => $XARGSFILEHANDLE,
                }
            );
            print $XARGSFILEHANDLE "; ";
        }

        gnu_mv(
            {
                infile_path => $outfile_path_prefix . "_"
                  . $contig
                  . $alt_file_tag
                  . $outfile_suffix,
                outfile_path => $outfile_path_prefix . "_"
                  . $contig
                  . $outfile_suffix,
                FILEHANDLE => $XARGSFILEHANDLE,
            }
        );
        say {$XARGSFILEHANDLE} "\n";
    }

    if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                infile_path => $outfile_path_prefix . q{_*}
                  . $outfile_suffix . q{*},
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} q{wait}, "\n";

        close $FILEHANDLE;
    }
    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

            slurm_submit_job_sample_id_dependency_add_to_family(
                {
                    job_id_href             => $job_id_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    sample_ids_ref =>
                      \@{ $active_parameter_href->{sample_ids} },
                    family_id        => $$family_id_ref,
                    path             => $job_id_chain,
                    log              => $log,
                    sbatch_file_name => $file_path,
                }
            );
        }
    }
    if ($$reduce_io_ref) {

        return
          $xargs_file_counter
          ; #Track the number of created xargs scripts per module for Block algorithm
    }
}

sub rhocall {

##rhocall

##Function : Rhocall performs annotation of autozygosity regions
##Returns  : "|$xargs_file_counter"
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, $program_info_path, $file_path, $stderr_path, $FILEHANDLE, family_id_ref, $temp_directory_ref, $outaligner_dir_ref, $call_type, $xargs_file_counter
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $program_info_path          => The program info path
##         : $file_path                  => File path
##         : $stderr_path                => The stderr path of the block script
##         : $FILEHANDLE                 => Filehandle to write to
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type
##         : $xargs_file_counter         => The xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
    my $call_type;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;
    my $program_info_path;
    my $file_path;
    my $stderr_path;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        file_path         => { strict_type => 1, store => \$file_path },
        stderr_path       => { strict_type => 1, store => \$stderr_path },
        FILEHANDLE    => { store => \$FILEHANDLE, },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw(get_core_number);
    use MIP::IO::Files qw(migrate_file xargs_migrate_contig_files);
    use MIP::Set::File qw{set_file_suffix};
    use MIP::Get::File qw{get_file_suffix};
    use MIP::Recipes::Xargs qw{ xargs_command };
    use MIP::Program::Variantcalling::Bcftools qw(bcftools_roh);
    use Program::Variantcalling::Rhocall qw(aggregate annotate);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_family);

    my $reduce_io_ref = \$active_parameter_href->{reduce_io};
    my $xargs_file_path_prefix;
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $XARGSFILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    unless ( defined($FILEHANDLE) ) {           #Run as individual sbatch script

        $FILEHANDLE = IO::Handle->new();        #Create anonymous filehandle
    }

    ## Get core number depending on user supplied input exists or not and max number of cores
    my $core_number = get_core_number(
        {
            module_core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            modifier_core_number => scalar( @{ $file_info_href->{contigs} } ),
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

        use MIP::Script::Setup_script qw(setup_script);
        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        ( $file_path, $program_info_path ) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                job_id_href           => $job_id_href,
                FILEHANDLE            => $FILEHANDLE,
                directory_id          => $$family_id_ref,
                program_name          => $program_name,
                program_directory     => catfile( lc($$outaligner_dir_ref) ),
                call_type             => $call_type,
                core_number           => $core_number,
                process_time =>
                  $active_parameter_href->{module_time}{ "p" . $program_name },
                temp_directory => $$temp_directory_ref,
            }
        );
        $stderr_path = $program_info_path . ".stderr.txt";
    }
    my ( $volume, $directory, $stderr_file ) =
      File::Spec->splitpath($stderr_path)
      ;    #Split to enable submission to &sample_info_qc later

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    $parameter_href->{ "p" . $program_name }{indirectory} =
      $outfamily_directory;    #Used downstream

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$$family_id_ref}{pgatk_combinevariantcallsets}
      {file_tag};
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{ "p" . $program_name }{file_tag};
    my $infile_prefix       = $$family_id_ref . $infile_tag . $call_type;
    my $file_path_prefix    = catfile( $$temp_directory_ref, $infile_prefix );
    my $outfile_prefix      = $$family_id_ref . $outfile_tag . $call_type;
    my $outfile_path_prefix = catfile( $$temp_directory_ref, $outfile_prefix );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            jobid_chain    => $job_id_chain,
        }
    );
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        ($xargs_file_counter) = xargs_migrate_contig_files(
            {
                FILEHANDLE      => $FILEHANDLE,
                XARGSFILEHANDLE => $XARGSFILEHANDLE,
                contigs_ref => \@{ $file_info_href->{contigs_size_ordered} },
                file_path   => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
                xargs_file_counter => $xargs_file_counter,
                infile             => $infile_prefix,
                indirectory        => $infamily_directory,
                temp_directory     => $$temp_directory_ref,
            }
        );
    }

    say {$FILEHANDLE} "## bcftools rho calculation";
    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## Get parameters
        my @sample_ids;
        if (   ( defined( $parameter_href->{dynamic_parameter}{affected} ) )
            && ( @{ $parameter_href->{dynamic_parameter}{affected} } ) )
        {

            push( @sample_ids,
                $parameter_href->{dynamic_parameter}{affected}[0] );
        }
        else {

            push( @sample_ids, $active_parameter_href->{sample_ids}[0] )
              ;    #No affected - pick any sample_id
        }

        bcftools_roh(
            {
                infile_path => $file_path_prefix . "_"
                  . $contig
                  . $infile_suffix,
                outfile_path => $file_path_prefix . "_" . $contig . ".roh",
                af_file_path =>
                  $active_parameter_href->{rhocall_frequency_file},
                sample_ids_ref => \@sample_ids,
                skip_indels =>
                  1,    #Skip indels as their genotypes are enriched for errors
                FILEHANDLE => $XARGSFILEHANDLE,
            }
        );
        print $XARGSFILEHANDLE "; ";

        Program::Variantcalling::Rhocall::annotate(
            {
                infile_path => $file_path_prefix . "_"
                  . $contig
                  . $infile_suffix,
                outfile_path => $outfile_path_prefix . "_"
                  . $contig
                  . $outfile_suffix,
                rohfile_path => $file_path_prefix . "_" . $contig . ".roh",
                v14          => 1,
                FILEHANDLE   => $XARGSFILEHANDLE,
            }
        );
        say {$XARGSFILEHANDLE} "\n";
    }

    if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                infile_path => $outfile_path_prefix . q{_*}
                  . $outfile_suffix . q{*},
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} q{wait}, "\n";

        close $FILEHANDLE;
    }
    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

            slurm_submit_job_sample_id_dependency_add_to_family(
                {
                    job_id_href             => $job_id_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    sample_ids_ref =>
                      \@{ $active_parameter_href->{sample_ids} },
                    family_id        => $$family_id_ref,
                    path             => $job_id_chain,
                    log              => $log,
                    sbatch_file_name => $file_path,
                }
            );
        }
    }
    if ($$reduce_io_ref) {

        return
          $xargs_file_counter
          ; #Track the number of created xargs scripts per module for Block algorithm
    }
}

sub prepareforvariantannotationblock {

##prepareforvariantannotationblock

##Function : Copy files for variantannotationblock to enable restart and skip of modules within block
##Returns  : "|$xargs_file_counter"
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, $program_info_path, $file_path, $stderr_path, $FILEHANDLE, family_id_ref, $temp_directory_ref, $outaligner_dir_ref, $call_type, $xargs_file_counter
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $program_info_path          => The program info path
##         : $file_path                  => File path
##         : $stderr_path                => The stderr path of the block script
##         : $FILEHANDLE                 => Filehandle to write to
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type
##         : $xargs_file_counter         => The xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
    my $call_type;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;
    my $program_info_path;
    my $file_path;
    my $stderr_path;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        file_path         => { strict_type => 1, store => \$file_path },
        stderr_path       => { strict_type => 1, store => \$stderr_path },
        FILEHANDLE    => { store => \$FILEHANDLE, },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw(get_core_number);
    use MIP::IO::Files qw(migrate_files);
    use MIP::Set::File qw{set_file_suffix};
    use MIP::Get::File qw{get_file_suffix};
    use MIP::Recipes::Xargs qw{ xargs_command };
    use Program::Htslib qw(bgzip tabix);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_family);

    my $reduce_io_ref = \$active_parameter_href->{reduce_io};
    my $xargs_file_path_prefix;
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $XARGSFILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    unless ( defined($FILEHANDLE) ) {           #Run as individual sbatch script

        $FILEHANDLE = IO::Handle->new();        #Create anonymous filehandle
    }

    ## Get core number depending on user supplied input exists or not and max number of cores
    my $core_number = get_core_number(
        {
            module_core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            modifier_core_number => scalar( @{ $file_info_href->{contigs} } ),
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

        use MIP::Script::Setup_script qw(setup_script);

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        ( $file_path, $program_info_path ) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                job_id_href           => $job_id_href,
                FILEHANDLE            => $FILEHANDLE,
                directory_id          => $$family_id_ref,
                program_name          => $program_name,
                program_directory     => catfile( lc($$outaligner_dir_ref) ),
                call_type             => $call_type,
                core_number           => $core_number,
                process_time =>
                  $active_parameter_href->{module_time}{ "p" . $program_name },
                temp_directory => $$temp_directory_ref,
            }
        );
        $stderr_path = $program_info_path . ".stderr.txt";
    }
    my ( $volume, $directory, $stderr_file ) =
      File::Spec->splitpath($stderr_path)
      ;    #Split to enable submission to &sample_info_qc later

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    $parameter_href->{ "p" . $program_name }{indirectory} =
      $outfamily_directory;    #Used downstream in removal of files

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$$family_id_ref}{pgatk_combinevariantcallsets}
      {file_tag};
    my $infile_prefix = $$family_id_ref . $infile_tag . $call_type;
    my $file_path_prefix = catfile( $$temp_directory_ref, $infile_prefix );

    ## Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            jobid_chain    => $job_id_chain,
        }
    );
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => catfile(
                $infamily_directory, $infile_prefix . $infile_suffix . q{*}
            ),
            outfile_path => $$temp_directory_ref
        }
    );
    say {$FILEHANDLE} q{wait}, "\n";

    ## Compress or decompress original file or stream to outfile (if supplied)
    bgzip(
        {
            FILEHANDLE      => $FILEHANDLE,
            infile_path     => $file_path_prefix . $infile_suffix,
            outfile_path    => $file_path_prefix . $outfile_suffix,
            write_to_stdout => 1,
        }
    );
    say {$FILEHANDLE} "\n";

    ## Index file using tabix
    tabix(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $file_path_prefix . $outfile_suffix,
            force       => 1,
            preset      => "vcf",
        }
    );
    say {$FILEHANDLE} "\n";

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            file_path          => $file_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Split vcf into contigs
    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        tabix(
            {
                regions_ref => [$contig],
                infile_path => $file_path_prefix . $outfile_suffix,
                with_header => 1,
                FILEHANDLE  => $XARGSFILEHANDLE,
            }
        );
        print $XARGSFILEHANDLE "| ";

        ## Compress or decompress original file or stream to outfile (if supplied)
        bgzip(
            {
                FILEHANDLE   => $XARGSFILEHANDLE,
                outfile_path => $file_path_prefix . "_"
                  . $contig
                  . $outfile_suffix,
                write_to_stdout => 1,
            }
        );
        print $XARGSFILEHANDLE "; ";

        ## Index file using tabix
        tabix(
            {
                FILEHANDLE  => $XARGSFILEHANDLE,
                infile_path => $file_path_prefix . "_"
                  . $contig
                  . $outfile_suffix,
                force  => 1,
                preset => "vcf",
            }
        );
        print $XARGSFILEHANDLE "\n";
    }

    if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                infile_path => $file_path_prefix . q{_*}
                  . $infile_suffix . q{*},
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} q{wait}, "\n";

        close $FILEHANDLE;
    }
    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

            slurm_submit_job_sample_id_dependency_add_to_family(
                {
                    job_id_href             => $job_id_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    sample_ids_ref =>
                      \@{ $active_parameter_href->{sample_ids} },
                    family_id        => $$family_id_ref,
                    path             => $job_id_chain,
                    log              => $log,
                    sbatch_file_name => $file_path,
                }
            );
        }
    }
    if ($$reduce_io_ref) {

        return
          $xargs_file_counter
          ; #Track the number of created xargs scripts per module for Block algorithm
    }
}

sub gatk_combinevariantcallsets {

##gatk_combinevariantcallsets

##Function : GATK CombineVariants to combine all variants call from different callers.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, family_id_ref, $temp_directory_ref, $outaligner_dir_ref, $call_type
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
    my $call_type;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use MIP::IO::Files qw(migrate_file);
    use MIP::Set::File qw{set_file_suffix};
    use MIP::Get::File qw{get_file_suffix};
    use MIP::Language::Java qw{java_core};
    use MIP::Program::Variantcalling::Gatk qw(gatk_combinevariants);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_family);

    my @variant_callers;    #Stores callers that have been executed
    my @parallel_chains
      ;    #Stores the parallel chains that jobIds should be inherited from
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory     => catfile( lc($$outaligner_dir_ref) ),
            call_type             => $call_type,
            core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
            temp_directory => $$temp_directory_ref,
        }
    );

    ## Assign directories
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    $parameter_href->{ "p" . $program_name }{indirectory} =
      $outfamily_directory;    #Used downstream

    ## Assign file_tags
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{pgatk_combinevariantcallsets}
      {file_tag};
    my @file_tags_and_paths;
    my $outfile_prefix = $$family_id_ref . $outfile_tag . $call_type;
    my $outfile_path_prefix = catfile( $$temp_directory_ref, $outfile_prefix );

    ### Assign suffix
    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    ## Collect file info and migrate files
    foreach my $variant_caller (
        @{ $parameter_href->{dynamic_parameter}{variant_callers} } )
    {

        if ( $active_parameter_href->{$variant_caller} > 0 ) {    #Expect vcf

            ## Assign directories
            my $program_outdirectory_name =
              $parameter_href->{$variant_caller}{outdir_name};
            my $infamily_directory =
              catfile( $active_parameter_href->{outdata_dir},
                $$family_id_ref, $$outaligner_dir_ref,
                $program_outdirectory_name );

            ## Assign file_tags
            my $infile_tag =
              $file_info_href->{$$family_id_ref}{$variant_caller}{file_tag};
            my $infile_prefix = $$family_id_ref . $infile_tag . $call_type;

            ## Assign suffix
            my $infile_suffix = get_file_suffix(
                {
                    parameter_href => $parameter_href,
                    suffix_key     => "outfile_suffix",
                    program_name   => $variant_caller,
                }
            );

            push(
                @file_tags_and_paths,
                $program_outdirectory_name . " "
                  . catfile(
                    $$temp_directory_ref, $infile_prefix . $infile_suffix
                  )
            );    #Collect both tag and path in the same string

            unshift( @variant_callers, $program_outdirectory_name )
              ; #To prioritize downstream - 1. gatk 2. samtools determined by order_parameters order

            unless ( $parameter_href->{$variant_caller}{chain} eq "MAIN" )
            {    #Do not add MAIN chains

                push( @parallel_chains,
                    $parameter_href->{$variant_caller}{chain} );
            }

            ## Copy file(s) to temporary directory
            say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
            migrate_file(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => catfile(
                        $infamily_directory,
                        $infile_prefix . $infile_suffix . q{*}
                    ),
                    outfile_path => $$temp_directory_ref
                }
            );
        }
    }
    say {$FILEHANDLE} q{wait}, "\n";

    ## GATK CombineVariants
    say {$FILEHANDLE} "## GATK CombineVariants";

    gatk_combinevariants(
        {
            memory_allocation => q{Xmx2g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $$temp_directory_ref,
            java_jar       => catfile(
                $active_parameter_href->{gatk_path},
                "GenomeAnalysisTK.jar"
            ),
            infile_paths_ref => \@file_tags_and_paths,
            outfile_path     => $outfile_path_prefix . $outfile_suffix,
            logging_level    => $active_parameter_href->{gatk_logging_level},
            referencefile_path =>
              $active_parameter_href->{human_genome_reference},
            genotype_merge_option => $active_parameter_href
              ->{gatk_combinevariants_genotype_merge_option},
            prioritize_caller =>
              $active_parameter_href->{gatk_combinevariants_prioritize_caller},
            exclude_nonvariants => 1,
            FILEHANDLE          => $FILEHANDLE
        }
    );
    say {$FILEHANDLE} "\n";

    if ( $active_parameter_href->{gatk_combinevariantcallsets_bcf_file} ) {

        ## Reformat variant calling file and index
        view_vcf(
            {
                infile_path         => $outfile_path_prefix . $outfile_suffix,
                outfile_path_prefix => $outfile_path_prefix,
                output_type         => "b",
                FILEHANDLE          => $FILEHANDLE,
            }
        );

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                infile_path  => $outfile_path_prefix . q{.bcf*},
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
    }

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            infile_path  => $outfile_path_prefix . $outfile_suffix,
            outfile_path => $outfamily_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, "\n";

    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        $sample_info_href->{vcf_file}{ready_vcf}{path} =
          catfile( $outfamily_directory, $outfile_prefix . $outfile_suffix );

        if ( $active_parameter_href->{gatk_combinevariantcallsets_bcf_file} ) {

            $sample_info_href->{bcf_file}{path} =
              catfile( $outfamily_directory, $outfile_prefix . ".bcf" );
        }

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref => \@{ $active_parameter_href->{sample_ids} },
                parallel_chains_ref => \@parallel_chains,
                family_id           => $$family_id_ref,
                path                => $job_id_chain,
                log                 => $log,
                sbatch_file_name    => $file_path,
            }
        );
    }
}

sub gatk_variantrecalibration {

##gatk_variantrecalibration

##Function : GATK VariantRecalibrator/ApplyRecalibration.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, family_id_ref, $temp_directory_ref, $outaligner_dir_ref, $call_type
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
    my $call_type;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use MIP::IO::Files qw(migrate_file);
    use MIP::Set::File qw{set_file_suffix};
    use MIP::Get::File qw{get_file_suffix};
    use MIP::Delete::List qw{ delete_contig_elements };
    use MIP::Gnu::Coreutils qw(gnu_mv);
    use MIP::Language::Java qw{java_core};
    use MIP::Program::Variantcalling::Bcftools qw(bcftools_norm);
    use MIP::Program::Variantcalling::Gatk
      qw(gatk_variantrecalibrator gatk_applyrecalibration gatk_selectvariants gatk_calculategenotypeposteriors);
    use MIP::QC::Record qw(add_program_outfile_to_sample_info);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_family);

    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Assign directories
    my $program_outdirectory_name =
      $parameter_href->{ "p" . $program_name }{outdir_name};
    my $outfamily_file_directory =
      catdir( $active_parameter_href->{outdata_dir}, $$family_id_ref )
      ;                                    #For ".fam" file
    my $infamily_directory = catfile( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref, $program_outdirectory_name );
    my $outfamily_directory = catfile( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref, $program_outdirectory_name );
    $parameter_href->{ "p" . $program_name }{indirectory} =
      $outfamily_directory;                #Used downstream

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$$family_id_ref}{pgatk_genotypegvcfs}{file_tag};
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{ "p" . $program_name }{file_tag};
    my $infile_prefix       = $$family_id_ref . $infile_tag . $call_type;
    my $file_path_prefix    = catfile( $$temp_directory_ref, $infile_prefix );
    my $outfile_prefix      = $$family_id_ref . $outfile_tag . $call_type;
    my $outfile_path_prefix = catfile( $$temp_directory_ref, $outfile_prefix );

    ### Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "outfile_suffix",
            program_name   => "pgatk_genotypegvcfs",
        }
    );

    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory     => catfile(
                lc($$outaligner_dir_ref),
                lc($program_outdirectory_name)
            ),
            call_type   => $call_type,
            core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
            temp_directory => catfile($$temp_directory_ref),
        }
    );
    my ( $volume, $directory, $stderr_file ) =
      File::Spec->splitpath( $program_info_path . ".stderr.txt" )
      ;    #Split to enable submission to &sample_info_qc later

    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            parameter_href        => $parameter_href,
            active_parameter_href => $active_parameter_href,
            sample_info_href      => $sample_info_href,
            FILEHANDLE            => $FILEHANDLE,
            fam_file_path =>
              catfile( $outfamily_file_directory, $$family_id_ref . ".fam" ),
        }
    );

    ## Check if "--pedigree" and "--pedigreeValidationType" should be included in analysis
    my %commands = gatk_pedigree_flag(
        {
            active_parameter_href => $active_parameter_href,
            fam_file_path =>
              catfile( $outfamily_file_directory, $$family_id_ref . ".fam" ),
            program_name => $program_name,
        }
    );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => catfile(
                $infamily_directory, $infile_prefix . $infile_suffix . q{*}
            ),
            outfile_path => $$temp_directory_ref
        }
    );
    say {$FILEHANDLE} q{wait}, "\n";

    ### GATK VariantRecalibrator
    ## Set mode to be used in variant recalibration
    my @modes = ( "SNP", "INDEL" );

    if (   ( $consensus_analysis_type eq "wes" )
        || ( $consensus_analysis_type eq "rapid" ) )
    {    #Exome/rapid analysis

        @modes = (q{BOTH});
    }

    foreach my $mode (@modes)
    { #SNP and INDEL will be recalibrated successively in the same file because when you specify eg SNP mode, the indels are emitted without modification, and vice-versa. Exome and rapid will be processed using mode BOTH since there are to few INDELS to use in the recalibration model even though using 30 exome BAMS in Haplotypecaller step.

        say {$FILEHANDLE} "## GATK VariantRecalibrator";

        ##Get parameters
        my @infiles;
        my $max_gaussian_level;
        if (   ( $consensus_analysis_type eq "wes" )
            || ( $consensus_analysis_type eq "rapid" ) )
        {    #Exome/rapid analysis use combined reference for more power

            push( @infiles, $file_path_prefix . $infile_suffix )
              ; #Infile HaplotypeCaller combined vcf which used reference gVCFs to create combined vcf (30> samples gCVFs)
        }
        else {    #WGS

            if ( $mode eq "SNP" ) {

                push( @infiles, $file_path_prefix . $infile_suffix );

                if ( $active_parameter_href
                    ->{gatk_variantrecalibration_snv_max_gaussians} ne 0 )
                {

                    $max_gaussian_level = 4;    #Use hard filtering
                }
            }
            if ( $mode eq "INDEL" ) { #Use created recalibrated snp vcf as input

                push( @infiles,
                    $outfile_path_prefix . ".SNV" . $infile_suffix );
            }
        }

        my @annotations =
          @{ $active_parameter_href->{gatk_variantrecalibration_annotations} };
        if ( $consensus_analysis_type ne "wgs" ) {

            ### Special case: Not to be used with hybrid capture. NOTE: Disable when analysing wes + wgs in the same run
            ## Removes an element from array and return new array while leaving orginal elements_ref untouched
            @annotations = delete_contig_elements(
                {
                    elements_ref       => \@annotations,
                    remove_contigs_ref => [qw{ DP }],
                }
            );
        }

        my @resources;
        if ( ( $mode eq "SNP" ) || ( $mode eq q{BOTH} ) ) {

            while (
                my ( $resource_file, $resource_string ) = each(
                    %{
                        $active_parameter_href
                          ->{gatk_variantrecalibration_resource_snv}
                    }
                )
              )
            {

                push( @resources, $resource_string . " " . $resource_file );
            }
        }
        if ( ( $mode eq "INDEL" ) || ( $mode eq q{BOTH} ) ) {

            while (
                my ( $resource_file, $resource_string ) = each(
                    %{
                        $active_parameter_href
                          ->{gatk_variantrecalibration_resource_indel}
                    }
                )
              )
            {

                push( @resources, $resource_string . " " . $resource_file );
            }
            if ( $active_parameter_href
                ->{gatk_variantrecalibration_indel_max_gaussians} ne 0 )
            {

                $max_gaussian_level = 4;    #Use hard filtering
            }
        }

        # Create distinct set i.e. no duplicates.
        @resources = uniq(@resources);

        gatk_variantrecalibrator(
            {
                memory_allocation => "Xmx10g",
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                temp_directory => $$temp_directory_ref,
                java_jar       => catfile(
                    $active_parameter_href->{gatk_path},
                    "GenomeAnalysisTK.jar"
                ),
                infile_paths_ref => \@infiles,
                annotations_ref  => \@annotations,
                resources_ref    => \@resources,
                logging_level => $active_parameter_href->{gatk_logging_level},
                referencefile_path =>
                  $active_parameter_href->{human_genome_reference},
                recal_file_path    => $file_path_prefix . ".intervals",
                rscript_file_path  => $file_path_prefix . ".intervals.plots.R",
                tranches_file_path => $file_path_prefix . ".intervals.tranches",
                max_gaussian_level => $max_gaussian_level,
                mode               => $mode,
                pedigree_validation_type => $commands{pedigree_validation_type},
                pedigree                 => $commands{pedigree},
                FILEHANDLE               => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} "\n";

        ## GATK ApplyRecalibration
        say {$FILEHANDLE} "## GATK ApplyRecalibration";

        ## Get parameters
        my $infile_path;
        my $outfile_path;
        my $ts_filter_level;
        if (   ( $consensus_analysis_type eq "wes" )
            || ( $consensus_analysis_type eq "rapid" ) )
        {    #Exome/rapid analysis use combined reference for more power

            $infile_path = $file_path_prefix
              . $infile_suffix
              ; #Infile genotypegvcfs combined vcf which used reference gVCFs to create combined vcf file
            $outfile_path =
              $outfile_path_prefix . "_filtered" . $outfile_suffix;
            $ts_filter_level = $active_parameter_href
              ->{gatk_variantrecalibration_snv_tsfilter_level};
        }
        else {    #WGS

            if ( $mode eq "SNP" ) {

                $infile_path  = $file_path_prefix . $infile_suffix;
                $outfile_path = $outfile_path_prefix . ".SNV" . $outfile_suffix;
                $ts_filter_level = $active_parameter_href
                  ->{gatk_variantrecalibration_snv_tsfilter_level};
            }
            if ( $mode eq "INDEL" ) { #Use created recalibrated snp vcf as input

                $infile_path  = $outfile_path_prefix . ".SNV" . $outfile_suffix;
                $outfile_path = $outfile_path_prefix . $outfile_suffix;
                $ts_filter_level = $active_parameter_href
                  ->{gatk_variantrecalibration_indel_tsfilter_level};
            }
        }

        gatk_applyrecalibration(
            {
                memory_allocation => "Xmx10g",
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                temp_directory => $$temp_directory_ref,
                java_jar       => catfile(
                    $active_parameter_href->{gatk_path},
                    "GenomeAnalysisTK.jar"
                ),
                infile_path   => $infile_path,
                outfile_path  => $outfile_path,
                logging_level => $active_parameter_href->{gatk_logging_level},
                referencefile_path =>
                  $active_parameter_href->{human_genome_reference},
                recal_file_path    => $file_path_prefix . ".intervals",
                tranches_file_path => $file_path_prefix . ".intervals.tranches",
                ts_filter_level    => $ts_filter_level,
                mode               => $mode,
                pedigree_validation_type => $commands{pedigree_validation_type},
                pedigree                 => $commands{pedigree},
                FILEHANDLE               => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} "\n";
    }

    if (   ( $consensus_analysis_type eq "wes" )
        || ( $consensus_analysis_type eq "rapid" ) )
    {

        ## BcfTools norm, Left-align and normalize indels, split multiallelics
        bcftools_norm(
            {
                FILEHANDLE => $FILEHANDLE,
                reference_path =>
                  $active_parameter_href->{human_genome_reference},
                infile_path => $outfile_path_prefix
                  . "_filtered"
                  . $outfile_suffix,
                output_type  => "v",
                outfile_path => $outfile_path_prefix
                  . "_filtered_normalized"
                  . $outfile_suffix,
                multiallelic    => "-",
                stderrfile_path => $outfile_path_prefix
                  . "_filtered_normalized.stderr",
            }
        );
        say {$FILEHANDLE} "\n";
    }

    ### GATK SelectVariants

    ## Removes all genotype information for exome ref and recalulates meta-data info for remaining samples in new file.
    if (   ( $consensus_analysis_type eq "wes" )
        || ( $consensus_analysis_type eq "rapid" ) )
    {    #Exome/rapid analysis

        say {$FILEHANDLE} "## GATK SelectVariants";

        gatk_selectvariants(
            {
                memory_allocation => q{Xmx2g},
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                temp_directory => $$temp_directory_ref,
                java_jar       => catfile(
                    $active_parameter_href->{gatk_path},
                    "GenomeAnalysisTK.jar"
                ),
                sample_names_ref => \@{ $active_parameter_href->{sample_ids} },
                logging_level => $active_parameter_href->{gatk_logging_level},
                referencefile_path =>
                  $active_parameter_href->{human_genome_reference},
                infile_path => $outfile_path_prefix
                  . "_filtered_normalized"
                  . $outfile_suffix,
                outfile_path        => $outfile_path_prefix . $outfile_suffix,
                exclude_nonvariants => 1,
                FILEHANDLE          => $FILEHANDLE,
            }
        );
        print {$FILEHANDLE} " &";

        ## Produces another vcf file containing non-variant loci (useful for example in MAF comparisons), but is not used downstream in MIP
        if ( $active_parameter_href
            ->{gatk_variantrecalibration_exclude_nonvariants_file} eq 1 )
        {

            say {$FILEHANDLE} "\n\n#GATK SelectVariants", "\n";

            gatk_selectvariants(
                {
                    memory_allocation => q{Xmx2g},
                    java_use_large_pages =>
                      $active_parameter_href->{java_use_large_pages},
                    temp_directory => $$temp_directory_ref,
                    java_jar       => catfile(
                        $active_parameter_href->{gatk_path},
                        "GenomeAnalysisTK.jar"
                    ),
                    sample_names_ref =>
                      \@{ $active_parameter_href->{sample_ids} },
                    logging_level =>
                      $active_parameter_href->{gatk_logging_level},
                    referencefile_path =>
                      $active_parameter_href->{human_genome_reference},
                    infile_path => $outfile_path_prefix
                      . "_filtered"
                      . $outfile_suffix,
                    outfile_path => $outfile_path_prefix
                      . "_incnonvariantloci"
                      . $outfile_suffix,
                    FILEHANDLE => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} "\n\nwait\n";

            ## Copies file from temporary directory.
            say {$FILEHANDLE} q{## Copy file from temporary directory};
            migrate_file(
                {
                    infile_path => $outfile_path_prefix
                      . q{_incnonvariantloci}
                      . $outfile_suffix . q{*},
                    outfile_path => $outfamily_directory,
                    FILEHANDLE   => $FILEHANDLE,
                }
            );
        }
        say {$FILEHANDLE} "\n\nwait\n";
    }

    ## GenotypeRefinement
    if ( $parameter_href->{dynamic_parameter}{trio} ) {

        say {$FILEHANDLE} "## GATK CalculateGenotypePosteriors";

        gatk_calculategenotypeposteriors(
            {
                memory_allocation => "Xmx6g",
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                temp_directory => $$temp_directory_ref,
                java_jar       => catfile(
                    $active_parameter_href->{gatk_path},
                    "GenomeAnalysisTK.jar"
                ),
                logging_level => $active_parameter_href->{gatk_logging_level},
                referencefile_path =>
                  $active_parameter_href->{human_genome_reference},
                infile_path  => $outfile_path_prefix . $outfile_suffix,
                outfile_path => $outfile_path_prefix
                  . "_refined"
                  . $outfile_suffix,
                supporting_callset_file_path => $active_parameter_href
                  ->{gatk_calculategenotypeposteriors_support_set},
                pedigree_validation_type => $commands{pedigree_validation_type},
                pedigree                 => $commands{pedigree},
                FILEHANDLE               => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} "\n";

        ## Change name of file to accomodate downstream
        gnu_mv(
            {
                infile_path => $outfile_path_prefix
                  . "_refined"
                  . $outfile_suffix,
                outfile_path => $outfile_path_prefix . $outfile_suffix,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} "\n";
    }

    ## BcfTools norm, Left-align and normalize indels, split multiallelics
    bcftools_norm(
        {
            FILEHANDLE     => $FILEHANDLE,
            reference_path => $active_parameter_href->{human_genome_reference},
            infile_path    => $outfile_path_prefix . $outfile_suffix,
            output_type    => "v",
            outfile_path   => $outfile_path_prefix
              . "_normalized"
              . $outfile_suffix,
            multiallelic => "-",
        }
    );
    say {$FILEHANDLE} "\n";

    ## Change name of file to accomodate downstream
    gnu_mv(
        {
            infile_path => $outfile_path_prefix
              . "_normalized"
              . $outfile_suffix,
            outfile_path => $outfile_path_prefix . $outfile_suffix,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

    ## Produce a bcf compressed and index from vcf
    if ( $active_parameter_href->{gatk_variantrecalibration_bcf_file} ) {

        ## Reformat variant calling file and index
        view_vcf(
            {
                infile_path         => $outfile_path_prefix . $outfile_suffix,
                outfile_path_prefix => $outfile_path_prefix,
                output_type         => "b",
                FILEHANDLE          => $FILEHANDLE,
            }
        );

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                infile_path  => $outfile_path_prefix . q{.bcf*},
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} q{wait}, "\n";
    }

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    my @outfiles = (
        $outfile_path_prefix . $outfile_suffix . "*",
        $file_path_prefix . ".intervals.tranches.pdf"
    );
    foreach my $outfile (@outfiles) {

        migrate_file(
            {
                infile_path  => $outfile,
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
    }
    say {$FILEHANDLE} q{wait}, "\n";

    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        ## Collect QC metadata info for later use
        # Disabled pedigreeCheck to not include relationship test is qccollect
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => q{pedigree_check},
                outdirectory     => $outfamily_directory,
                outfile          => $outfile_prefix . $outfile_suffix,
            }
        );

        $sample_info_href->{vcf_file}{ready_vcf}{path} =
          catfile( $outfamily_directory, $outfile_prefix . $outfile_suffix );

        if ( $active_parameter_href->{gatk_variantrecalibration_bcf_file} eq 1 )
        {

            $sample_info_href->{bcf_file}{path} =
              catfile( $outfamily_directory, $outfile_prefix . ".bcf" );
        }

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $$family_id_ref,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
}

sub gatk_concatenate_genotypegvcfs {

##gatk_concatenate_genotypegvcfs

##Function : Concatenate GVCFs produced after gatk_genotypegvcfs done per contig.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
    my $call_type;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Processes qw(print_wait);
    use MIP::Script::Setup_script qw(setup_script);
    use MIP::IO::Files qw(migrate_file);
    use MIP::Set::File qw{set_file_suffix};
    use MIP::Get::File qw{get_file_suffix};
    use MIP::Language::Java qw{java_core};
    use MIP::Program::Variantcalling::Gatk qw(gatk_selectvariants);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_family);

    my $core_number =
      $active_parameter_href->{module_core_number}{ "p" . $program_name };
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref, "gatk" );
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref, "gatk" );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$$family_id_ref}{ "p" . $program_name }{file_tag};
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{ "p" . $program_name }{file_tag};
    my $infile_prefix       = $$family_id_ref . $infile_tag . $call_type;
    my $file_path_prefix    = catfile( $$temp_directory_ref, $infile_prefix );
    my $outfile_prefix      = $$family_id_ref . $outfile_tag . $call_type;
    my $outfile_path_prefix = catfile( $$temp_directory_ref, $outfile_prefix );

    ### Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "outfile_suffix",
            program_name   => "p" . $program_name,
        }
    );

    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory => catfile( lc($$outaligner_dir_ref), "gatk" ),
            call_type         => $call_type,
            core_number       => $core_number,
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
            temp_directory => $$temp_directory_ref
        }
    );

    my $process_batches_count = 1;
    while ( my ( $contig_index, $contig ) =
        each( @{ $file_info_href->{contigs} } ) )
    {

        $process_batches_count = print_wait(
            {
                process_counter       => $contig_index,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                FILEHANDLE            => $FILEHANDLE,
            }
        );

        ## Copy file(s) to temporary directory
        migrate_file(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => catfile(
                    $infamily_directory,
                    $infile_prefix . q{_} . $contig . $infile_suffix . q{*}
                ),
                outfile_path => $$temp_directory_ref
            }
        );
    }
    say {$FILEHANDLE} q{wait}, "\n";

    ## Writes sbatch code to supplied filehandle to concatenate variants in vcf format. Each array element is combined with the infile prefix and postfix.
    concatenate_variants(
        {
            active_parameter_href => $active_parameter_href,
            FILEHANDLE            => $FILEHANDLE,
            elements_ref          => \@{ $file_info_href->{contigs} },
            infile_prefix         => $file_path_prefix . "_",
            infile_postfix        => $infile_suffix,
            outfile               => $outfile_path_prefix . $outfile_suffix,
        }
    );

    ## Produce a bcf compressed and index from vcf
    if ( $active_parameter_href->{gatk_concatenate_genotypegvcfs_bcf_file} ) {

        if (   ( $consensus_analysis_type eq "wes" )
            || ( $consensus_analysis_type eq "rapid" ) )
        {    #Exome/rapid analysis

            say {$FILEHANDLE} "###Remove extra reference samples", "\n";

            say {$FILEHANDLE} "##GATK SelectVariants", "\n";

            gatk_selectvariants(
                {
                    memory_allocation => q{Xmx2g},
                    java_use_large_pages =>
                      $active_parameter_href->{java_use_large_pages},
                    temp_directory => $$temp_directory_ref,
                    java_jar       => catfile(
                        $active_parameter_href->{gatk_path},
                        "GenomeAnalysisTK.jar"
                    ),
                    sample_names_ref =>
                      \@{ $active_parameter_href->{sample_ids} },
                    logging_level =>
                      $active_parameter_href->{gatk_logging_level},
                    referencefile_path =>
                      $active_parameter_href->{human_genome_reference},
                    infile_path  => $outfile_path_prefix . $outfile_suffix,
                    outfile_path => $outfile_path_prefix
                      . "_incnonvariantloci"
                      . $outfile_suffix,
                    FILEHANDLE => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} "\n";

            ## Move to original filename
            gnu_mv(
                {
                    infile_path => $outfile_path_prefix
                      . "_incnonvariantloci"
                      . $outfile_suffix,
                    outfile_path => $outfile_path_prefix . $outfile_suffix,
                    FILEHANDLE   => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} "\n";
        }

        ## Reformat variant calling file and index
        view_vcf(
            {
                infile_path         => $outfile_path_prefix . $outfile_suffix,
                outfile_path_prefix => $outfile_path_prefix,
                output_type         => "b",
                FILEHANDLE          => $FILEHANDLE,
            }
        );

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                infile_path  => $outfile_path_prefix . q{.bcf*},
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} q{wait}, "\n";
    }

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            infile_path  => $outfile_path_prefix . $outfile_suffix . q{*},
            outfile_path => $outfamily_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, "\n";

    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        if ( $active_parameter_href->{gatk_concatenate_genotypegvcfs_bcf_file}
            eq 1 )
        {

            $sample_info_href->{gbcf_file}{path} =
              catfile( $outfamily_directory, $outfile_prefix . ".bcf" );
        }

        ## Collect QC metadata info for later use
        $sample_info_href->{vcf_file}{ready_vcf}{path} =
          catfile( $outfamily_directory, $outfile_prefix . $outfile_suffix );

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $$family_id_ref,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
}

sub mgatk_genotypegvcfs {

##gatk_genotypegvcfs

##Function : GATK GenoTypeGVCFs.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $family_id_ref, $outaligner_dir_ref, $call_type, $program_name
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $family_id_ref              => Family id {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type
##         : $program_name               => Program name

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
    my $call_type;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use MIP::IO::Files qw(migrate_file);
    use MIP::Get::File qw{get_file_suffix get_merged_infile_prefix };
    use MIP::Set::File qw{set_file_suffix};
    use MIP::Language::Java qw{java_core};
    use MIP::Program::Variantcalling::Gatk qw(gatk_genotypegvcfs);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_step_in_parallel_to_family);

    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $sbatch_script_tracker = 0;
    my $core_number =
      $active_parameter_href->{module_core_number}{ "p" . $program_name }
      ; #gatk genotype is most safely processed in single thread mode, , but we need some java heap allocation
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    my $process_time =
      $active_parameter_href->{module_time}{ "p" . $program_name };
    if ( $active_parameter_href->{gatk_genotypegvcfs_all_sites} eq 1 ) {

        $process_time = 50; #Including all sites requires longer processing time
    }

    ## Filehandles
    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Assign directories
    my $outfamily_file_directory =
      catdir( $active_parameter_href->{outdata_dir}, $$family_id_ref )
      ;                                    #For ".fam" file
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref, "gatk" );
    $parameter_href->{ "p" . $program_name }{indirectory} =
      $outfamily_directory;                #Used downstream

    ## Assign file_tags
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{ "p" . $program_name }{file_tag};

    ### Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "outfile_suffix",
            program_name =>
              "pgatk_haplotypecaller",     #Attach to haplotypecaller suffix
        }
    );

    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            parameter_href        => $parameter_href,
            active_parameter_href => $active_parameter_href,
            sample_info_href      => $sample_info_href,
            FILEHANDLE            => $FILEHANDLE,
            fam_file_path =>
              catfile( $outfamily_file_directory, $$family_id_ref . ".fam" ),
        }
    );

    ## Split per contig
    foreach my $contig ( @{ $file_info_href->{contigs} } ) {

        ## Assign file_tags
        my $outfile_prefix =
          $$family_id_ref . $outfile_tag . $call_type . "_" . $contig;
        my $outfile_path_prefix =
          catfile( $$temp_directory_ref, $outfile_prefix );

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        my ($file_path) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                job_id_href           => $job_id_href,
                FILEHANDLE            => $FILEHANDLE,
                directory_id          => $$family_id_ref,
                program_name          => $program_name,
                program_directory =>
                  catfile( lc($$outaligner_dir_ref), "gatk" ),
                call_type      => $call_type,
                core_number    => $core_number,
                process_time   => $process_time,
                temp_directory => $$temp_directory_ref,
                sleep          => 1
                , #Let process sleep randomly for 0-60 seconds to avoid race condition
            }
        );

        ## Collect infiles for all sample_ids to enable migration to temporary directory
        my @file_paths;
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Add merged infile name prefix after merging all BAM files per sample_id
            my $merged_infile_prefix = get_merged_infile_prefix(
                {
                    file_info_href => $file_info_href,
                    sample_id      => $sample_id,
                }
            );

            ## Assign directories
            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $$outaligner_dir_ref, "gatk" );

            ## Assign file_tags
            my $infile_tag =
              $file_info_href->{$sample_id}{pgatk_haplotypecaller}{file_tag};
            my $infile_prefix =
              $merged_infile_prefix . $infile_tag . "_" . $contig;

            ## Collect for downstream use
            push(
                @file_paths,
                catfile(
                    $$temp_directory_ref, $infile_prefix . $infile_suffix
                )
            );

            ## Copy file(s) to temporary directory
            say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
            migrate_file(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => catfile(
                        $insample_directory,
                        $infile_prefix . $infile_suffix . q{*}
                    ),
                    outfile_path => $$temp_directory_ref
                }
            );
            say {$FILEHANDLE} q{wait}, "\n";
        }

        ## GATK GenoTypeGVCFs
        say {$FILEHANDLE} "## GATK GenoTypeGVCFs";

        ## Get parameters
        ## Check if "--pedigree" and "--pedigreeValidationType" should be included in analysis
        my %commands = gatk_pedigree_flag(
            {
                active_parameter_href => $active_parameter_href,
                fam_file_path         => catfile(
                    $outfamily_file_directory, $$family_id_ref . ".fam"
                ),
                program_name => $program_name,
            }
        );
        ## Files for genotypegvcfs
        if (   ( $consensus_analysis_type eq "wes" )
            || ( $consensus_analysis_type eq "rapid" ) )
        {

            push( @file_paths,
                $active_parameter_href->{gatk_genotypegvcfs_ref_gvcf} );
        }

        gatk_genotypegvcfs(
            {
                memory_allocation => "Xmx24g",
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                temp_directory => $$temp_directory_ref,
                java_jar       => catfile(
                    $active_parameter_href->{gatk_path},
                    "GenomeAnalysisTK.jar"
                ),
                intervals_ref    => [$contig],
                infile_paths_ref => \@file_paths,
                outfile_path     => $outfile_path_prefix . $outfile_suffix,
                logging_level => $active_parameter_href->{gatk_logging_level},
                referencefile_path =>
                  $active_parameter_href->{human_genome_reference},
                dbsnp_file_path =>
                  $active_parameter_href->{gatk_haplotypecaller_snp_known_set},
                pedigree_validation_type => $commands{pedigree_validation_type},
                pedigree                 => $commands{pedigree},
                include_nonvariant_sites =>
                  $active_parameter_href->{gatk_genotypegvcfs_all_sites},
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} "\n";

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                infile_path  => $outfile_path_prefix . $outfile_suffix . q{*},
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} q{wait}, "\n";

        close $FILEHANDLE;

        if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

            slurm_submit_job_sample_id_dependency_step_in_parallel_to_family(
                {
                    job_id_href             => $job_id_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    sample_ids_ref =>
                      \@{ $active_parameter_href->{sample_ids} },
                    family_id             => $$family_id_ref,
                    path                  => $job_id_chain,
                    log                   => $log,
                    sbatch_file_name      => $file_path,
                    sbatch_script_tracker => $sbatch_script_tracker
                }
            );
        }
        $sbatch_script_tracker++;    #Tracks nr of sbatch scripts
    }
}

sub sv_reformat {

##sv_reformat

##Function : Concatenate contig files.
##Returns  : "|$xargs_file_counter"
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, $program_info_path, $file_path, $FILEHANDLE, family_id_ref, $temp_directory_ref, $reference_dir_ref, $outaligner_dir_ref, $call_type, $xargs_file_counter
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $family_id_ref              => The family_id_ref {REF}
##         : $call_type                  => Variant call type
##         : $program_name               => Program name
##         : $FILEHANDLE                 => Sbatch filehandle to write to
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $reference_dir_ref          => MIP reference directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type
##         : $xargs_file_counter         => The xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $reference_dir_ref;
    my $outaligner_dir_ref;
    my $call_type;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        reference_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{reference_dir},
            strict_type => 1,
            store       => \$reference_dir_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => "SV", strict_type => 1, store => \$call_type, },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use MIP::Get::File qw{get_file_suffix};
    use MIP::Delete::List qw{ delete_contig_elements delete_male_contig };
    use MIP::IO::Files qw(migrate_file xargs_migrate_contig_files);
    use Program::Htslib qw(bgzip tabix);
    use MIP::Gnu::Software::Gnu_grep qw( gnu_grep);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_family);

    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $xargs_file_path_prefix;
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $FILEHANDLE      = IO::Handle->new();    #Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Set the number of cores
    my $core_number =
      $active_parameter_href->{module_core_number}{ "p" . $program_name };

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory     => catfile( lc($$outaligner_dir_ref) ),
            core_number           => $core_number,
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
            temp_directory => $$temp_directory_ref
        }
    );

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    my $outfamily_file_directory =
      catfile( $active_parameter_href->{outdata_dir}, $$family_id_ref );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$$family_id_ref}{psv_rankvariant}{file_tag};
    my $infile_prefix = $$family_id_ref . $infile_tag . $call_type;
    my $file_path_prefix = catfile( $$temp_directory_ref, $infile_prefix );
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{ "p" . $program_name }{file_tag};
    my $outfile_prefix      = $$family_id_ref . $outfile_tag . $call_type;
    my $outfile_path_prefix = catfile( $$temp_directory_ref,
        $$family_id_ref . $outfile_tag . $call_type );
    my $final_path_prefix = catfile( $outfamily_directory, $outfile_prefix );

    ## Assign suffix
    my $file_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            jobid_chain    => $job_id_chain,
        }
    );

    my $vcfparser_analysis_type = "";

    ## Removes an element from array and return new array while leaving orginal elements_ref untouched
    my @contigs = delete_contig_elements(
        {
            elements_ref       => \@{ $file_info_href->{contigs} },
            remove_contigs_ref => [qw{ MT M }],
        }
    );

    my @contigs_size_ordered = delete_contig_elements(
        {
            elements_ref       => \@{ $file_info_href->{contigs_size_ordered} },
            remove_contigs_ref => [qw{ MT M }],
        }
    );

    ### If no males or other remove contig Y from all downstream analysis
    my @contig_arrays = ( \@contigs_size_ordered, \@contigs );

    foreach my $array_ref (@contig_arrays) {

        ## Removes contig_names from contigs array if no male or other found
        $array_ref = delete_male_contig(
            {
                contigs_ref => $array_ref,
                found_male  => $active_parameter_href->{found_male},
            }
        );
    }

    for (
        my $vcfparser_outfile_counter = 0 ;
        $vcfparser_outfile_counter <
        $active_parameter_href->{vcfparser_outfile_count} ;
        $vcfparser_outfile_counter++
      )
    {

        if ( $vcfparser_outfile_counter == 1 ) {

            $vcfparser_analysis_type = ".selected";    #SelectFile variants

            @contigs = delete_contig_elements(
                {
                    elements_ref =>
                      \@{ $file_info_href->{select_file_contigs} },
                    remove_contigs_ref => [qw{ MT M }],
                }
            );

            ## Removes contigs from supplied contigs_ref
            remove_array_element(
                {
                    contigs_ref        => \@contigs,
                    remove_contigs_ref => ["Y"],
                }
            );

            ## Removes an element from array and return new array while leaving orginal elements_ref untouched
            @contigs_size_ordered = delete_contig_elements(
                {
                    elements_ref =>
                      \@{ $file_info_href->{sorted_select_file_contigs} },
                    remove_contigs_ref => [qw{ MT M }],
                }
            );

            ## Removes contigs from supplied contigs_ref
            remove_array_element(
                {
                    contigs_ref        => \@contigs_size_ordered,
                    remove_contigs_ref => ["Y"],
                }
            );
        }

        if (   ( $consensus_analysis_type eq "wgs" )
            || ( $consensus_analysis_type eq "mixed" ) )
        {    #Transfer contig files

            ## Copy file(s) to temporary directory
            say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
            $xargs_file_counter = xargs_migrate_contig_files(
                {
                    FILEHANDLE         => $FILEHANDLE,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    contigs_ref        => \@contigs_size_ordered,
                    file_path          => $file_path,
                    program_info_path  => $program_info_path,
                    core_number        => $core_number,
                    xargs_file_counter => $xargs_file_counter,
                    infile             => $infile_prefix,
                    file_ending        => $vcfparser_analysis_type
                      . $file_suffix . "*",
                    indirectory    => $infamily_directory,
                    temp_directory => $active_parameter_href->{temp_directory},
                }
            );
        }
        else {

            ## Copy file(s) to temporary directory
            say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
            migrate_file(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => catfile(
                        $infamily_directory,
                        $infile_prefix
                          . $vcfparser_analysis_type
                          . $file_suffix
                    ),
                    outfile_path => $$temp_directory_ref
                }
            );
            say {$FILEHANDLE} q{wait}, "\n";
        }

        my $concatenate_ending = "";
        if (   ( $consensus_analysis_type eq "wgs" )
            || ( $consensus_analysis_type eq "mixed" ) )
        {

            $concatenate_ending = "_cat";

            ## Writes sbatch code to supplied filehandle to concatenate variants in vcf format. Each array element is combined with the infile prefix and postfix.
            concatenate_variants(
                {
                    active_parameter_href => $active_parameter_href,
                    FILEHANDLE            => $FILEHANDLE,
                    elements_ref          => \@contigs,
                    infile_prefix         => $file_path_prefix . "_",
                    infile_postfix => $vcfparser_analysis_type . $file_suffix,
                    outfile        => $file_path_prefix
                      . $vcfparser_analysis_type
                      . $concatenate_ending
                      . $file_suffix,
                }
            );
        }

        ## Writes sbatch code to supplied filehandle to sort variants in vcf format
        sort_vcf(
            {
                active_parameter_href => $active_parameter_href,
                FILEHANDLE            => $FILEHANDLE,
                sequence_dict_file    => catfile(
                    $$reference_dir_ref,
                    $file_info_href->{human_genome_reference_name_prefix}
                      . ".dict"
                ),
                infile_paths_ref => [
                        $file_path_prefix
                      . $vcfparser_analysis_type
                      . $concatenate_ending
                      . $file_suffix
                ],
                outfile => $outfile_path_prefix
                  . $vcfparser_analysis_type
                  . $file_suffix,
            }
        );

        print {$FILEHANDLE} "\n";

        ## Remove variants in hgnc_id list from vcf
        if ( $active_parameter_href->{sv_reformat_remove_genes_file} ) {

            ## Removes contig_names from contigs array if no male or other found
            gnu_grep(
                {
                    filter_file_path => catfile(
                        $$reference_dir_ref,
                        $active_parameter_href->{sv_reformat_remove_genes_file}
                    ),
                    infile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $file_suffix,
                    outfile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . "_filtered"
                      . $file_suffix,
                    invert_match => 1,
                    FILEHANDLE   => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} "\n";

            if ( $vcfparser_outfile_counter == 1 ) {

                $sample_info_href->{program}{$program_name}
                  {sv_reformat_remove_genes_file}{clinical}{path} =
                    $final_path_prefix
                  . $vcfparser_analysis_type
                  . "_filtered"
                  . $file_suffix;    #Save filtered file
            }
            else {

                $sample_info_href->{program}{$program_name}
                  {sv_reformat_remove_genes_file}{research}{path} =
                    $final_path_prefix
                  . $vcfparser_analysis_type
                  . "_filtered"
                  . $file_suffix;    #Save filtered file
            }

            ## Copies file from temporary directory.
            say {$FILEHANDLE} q{## Copy file from temporary directory};
            migrate_file(
                {
                    infile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . q{_filtered}
                      . $file_suffix,
                    outfile_path => $outfamily_directory,
                    FILEHANDLE   => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} q{wait}, "\n";
        }

        if ( $active_parameter_href->{sv_rankvariant_binary_file} ) {

            ## Compress or decompress original file or stream to outfile (if supplied)
            bgzip(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $file_suffix,
                    outfile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $file_suffix . ".gz",
                    write_to_stdout => 1,
                }
            );
            say {$FILEHANDLE} "\n";

            ## Index file using tabix
            tabix(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $file_suffix . ".gz",
                    force  => 1,
                    preset => substr( $file_suffix, 1 ),
                }
            );
            say {$FILEHANDLE} "\n";
        }

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                infile_path => $outfile_path_prefix
                  . $vcfparser_analysis_type
                  . $file_suffix . q{*},
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} q{wait}, "\n";

        ## Adds the most complete vcf file to sample_info
        add_most_complete_vcf(
            {
                active_parameter_href => $active_parameter_href,
                sample_info_href      => $sample_info_href,
                program_name          => $program_name,
                path                  => $final_path_prefix
                  . $vcfparser_analysis_type
                  . $file_suffix,
                vcfparser_outfile_counter => $vcfparser_outfile_counter,
                vcf_file_key => "sv_" . substr( $file_suffix, 1 ) . "_file",
            }
        );

        if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

            if ( $vcfparser_outfile_counter == 1 ) {

                # Save clinical candidate list path
                my $clinical_candidate_path =
                  $final_path_prefix . $vcfparser_analysis_type . $file_suffix;
                add_program_metafile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        program_name     => $program_name,
                        metafile_tag     => q{clinical},
                        path             => $clinical_candidate_path,
                    }
                );

                if ( $active_parameter_href->{sv_rankvariant_binary_file} ) {

                    $sample_info_href->{sv_vcf_binary_file}{clinical}{path} =
                        $final_path_prefix
                      . $vcfparser_analysis_type
                      . $file_suffix . ".gz";
                }
            }
            else {

                # Save research candidate list path
                my $research_candidate_path =
                  $final_path_prefix . $vcfparser_analysis_type . $file_suffix;
                add_program_metafile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        program_name     => $program_name,
                        metafile_tag     => q{research},
                        path             => $research_candidate_path,
                    }
                );

                if ( $active_parameter_href->{sv_rankvariant_binary_file} ) {

                    $sample_info_href->{sv_vcf_binary_file}{research}{path} =
                        $final_path_prefix
                      . $vcfparser_analysis_type
                      . $file_suffix . ".gz";
                }
            }
        }
    }
    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $$family_id_ref,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
}

sub sv_rankvariant {

##sv_rankvariant

##Function : Annotate and score SV variants depending on mendelian inheritance, frequency and phenotype etc.
##Returns  : "|$xargs_file_counter"
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, $program_info_path, $FILEHANDLE, family_id_ref, $temp_directory_ref, $reference_dir_ref, $outaligner_dir_ref, $call_type, $xargs_file_counter
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $family_id_ref              => The family_id_ref {REF}
##         : $program_name               => Program name
##         : $FILEHANDLE                 => Sbatch filehandle to write to
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $reference_dir_ref          => MIP reference directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type
##         : $xargs_file_counter         => The xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $reference_dir_ref;
    my $outaligner_dir_ref;
    my $call_type;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        reference_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{reference_dir},
            strict_type => 1,
            store       => \$reference_dir_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => "SV", strict_type => 1, store => \$call_type, },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use MIP::Get::File qw{get_file_suffix};
    use MIP::Delete::List qw{ delete_contig_elements delete_male_contig};
    use MIP::Recipes::Xargs qw{ xargs_command };
    use MIP::IO::Files qw(migrate_file xargs_migrate_contig_files);
    use Program::Variantcalling::Genmod qw(annotate models score compound);
    use MIP::QC::Record
      qw(add_program_outfile_to_sample_info add_program_metafile_to_sample_info);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_family);

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $xargs_file_path_prefix;
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $FILEHANDLE      = IO::Handle->new();    #Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Set the number of cores
    my $core_number =
      $active_parameter_href->{module_core_number}{ "p" . $program_name };

    ## Limit number of cores requested to the maximum number of cores available per node
    my $genmod_core_number = check_max_core_number(
        {
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
            core_number_requested => $core_number,
        }
    );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory     => catfile( lc($$outaligner_dir_ref) ),
            core_number           => $core_number,
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
            temp_directory => $$temp_directory_ref
        }
    );

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    my $outfamily_file_directory =
      catfile( $active_parameter_href->{outdata_dir}, $$family_id_ref );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$$family_id_ref}{psv_vcfparser}{file_tag};
    my $infile_prefix = $$family_id_ref . $infile_tag . $call_type;
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{ "p" . $program_name }{file_tag};
    my $outfile_prefix = $$family_id_ref . $outfile_tag . $call_type;
    my $outfile_path_prefix = catfile( $$temp_directory_ref, $outfile_prefix );

    ## Assign suffix
    my $file_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            jobid_chain    => $job_id_chain,
        }
    );

    my $vcfparser_analysis_type = "";

    ## Removes an element from array and return new array while leaving orginal elements_ref untouched
    my @contigs_size_ordered = delete_contig_elements(
        {
            elements_ref       => \@{ $file_info_href->{contigs_size_ordered} },
            remove_contigs_ref => [qw{ MT M }],
        }
    );
    ## Removes an element from array and return new array while leaving orginal elements_ref untouched
    my @contigs = delete_contig_elements(
        {
            elements_ref       => \@{ $file_info_href->{contigs} },
            remove_contigs_ref => [qw{ MT M }],
        }
    );

    ### If no males or other remove contig Y from all downstream analysis
    my @contig_arrays = ( \@contigs_size_ordered, \@contigs );

    foreach my $array_ref (@contig_arrays) {

        ## Removes contig_names from contigs array if no male or other found
        $array_ref = delete_male_contig(
            {
                contigs_ref => $array_ref,
                found_male  => $active_parameter_href->{found_male},
            }
        );
    }

    my $family_file =
      catfile( $outfamily_file_directory, $$family_id_ref . ".fam" );
    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            parameter_href        => $parameter_href,
            active_parameter_href => $active_parameter_href,
            sample_info_href      => $sample_info_href,
            FILEHANDLE            => $FILEHANDLE,
            fam_file_path         => $family_file,
        }
    );

    for (
        my $vcfparser_outfile_counter = 0 ;
        $vcfparser_outfile_counter <
        $active_parameter_href->{vcfparser_outfile_count} ;
        $vcfparser_outfile_counter++
      )
    {

        if ( $vcfparser_outfile_counter == 1 ) {

            $vcfparser_analysis_type = ".selected";    #SelectFile variants

            ### Always skip MT and Y in select files
            ## Removes an element from array and return new array while leaving orginal elements_ref untouched
            @contigs = delete_contig_elements(
                {
                    elements_ref =>
                      \@{ $file_info_href->{select_file_contigs} },
                    remove_contigs_ref => [qw{ MT M }],
                }
            );

            ## Removes contigs from supplied contigs_ref
            remove_array_element(
                {
                    contigs_ref        => \@contigs,
                    remove_contigs_ref => ["Y"],
                }
            );

            ## Removes an element from array and return new array while leaving orginal elements_ref untouched
            @contigs_size_ordered = delete_contig_elements(
                {
                    elements_ref =>
                      \@{ $file_info_href->{sorted_select_file_contigs} },
                    remove_contigs_ref => [qw{ MT M }],
                }
            );

            ## Removes contigs from supplied contigs_ref
            remove_array_element(
                {
                    contigs_ref        => \@contigs_size_ordered,
                    remove_contigs_ref => ["Y"],
                }
            );
        }

        if (   ( $consensus_analysis_type eq "wgs" )
            || ( $consensus_analysis_type eq "mixed" ) )
        {    #Transfer contig files

            ## Copy file(s) to temporary directory
            say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
            $xargs_file_counter = xargs_migrate_contig_files(
                {
                    FILEHANDLE         => $FILEHANDLE,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    contigs_ref        => \@contigs_size_ordered,
                    file_path          => $file_path,
                    program_info_path  => $program_info_path,
                    core_number        => $core_number,
                    xargs_file_counter => $xargs_file_counter,
                    infile             => $infile_prefix,
                    file_ending        => $vcfparser_analysis_type
                      . $file_suffix . "*",
                    indirectory    => $infamily_directory,
                    temp_directory => $active_parameter_href->{temp_directory},
                }
            );
        }
        else {

            ## Copy file(s) to temporary directory
            say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
            migrate_file(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => catfile(
                        $infamily_directory,
                        $infile_prefix
                          . $vcfparser_analysis_type
                          . $file_suffix
                    ),
                    outfile_path => $$temp_directory_ref
                }
            );
            say {$FILEHANDLE} q{wait}, "\n";
        }

        my $genmod_module = "";   #Track which genmod modules has been processed

        ## Check affected/unaffected status
        if (
            ( defined( $parameter_href->{dynamic_parameter}{unaffected} ) )
            && ( @{ $parameter_href->{dynamic_parameter}{unaffected} } eq
                @{ $active_parameter_href->{sample_ids} } )
          )
        {                         #Only unaffected

            if ( !$vcfparser_outfile_counter ) {

                $log->warn(
"Only unaffected sample(s) in pedigree - skipping genmod 'models', 'score' and 'compound'"
                );
            }
        }

        ## Genmod
        say {$FILEHANDLE} "## Genmod";

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                FILEHANDLE         => $FILEHANDLE,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                file_path          => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $genmod_core_number,
                xargs_file_counter => $xargs_file_counter,
            }
        );

      CONTIG:
        foreach my $contig (@contigs_size_ordered) {

            ## Get parameters
            my $genmod_file_ending_stub       = $infile_prefix;
            my $genmod_outfile_path_prefix    = $outfile_path_prefix;
            my $genmod_xargs_file_path_prefix = $xargs_file_path_prefix;
            my $genmod_indata                 = catfile( $$temp_directory_ref,
                    $genmod_file_ending_stub
                  . $vcfparser_analysis_type
                  . $file_suffix );    #InFile

            # Update endings with contig info
            if (   ( $consensus_analysis_type eq "wgs" )
                || ( $consensus_analysis_type eq "mixed" ) )
            {

                $genmod_file_ending_stub = $infile_prefix . "_" . $contig;
                $genmod_outfile_path_prefix =
                  $outfile_path_prefix . "_" . $contig;
                $genmod_xargs_file_path_prefix =
                  $xargs_file_path_prefix . "." . $contig;

                # Infile
                $genmod_indata = catfile( $$temp_directory_ref,
                        $genmod_file_ending_stub
                      . $vcfparser_analysis_type
                      . $file_suffix );
            }
            $genmod_module = "";    #Restart for next contig

            ## Genmod Annotate
            $genmod_module = "_annotate";

            my $genmod_outfile_path;
            if (
                ( defined( $parameter_href->{dynamic_parameter}{unaffected} ) )
                && ( @{ $parameter_href->{dynamic_parameter}{unaffected} } eq
                    @{ $active_parameter_href->{sample_ids} } )
              )
            {                       #Only unaffected

                ## Write to outputFile - last genmod module
                # Outfile
                $genmod_outfile_path =
                    $genmod_outfile_path_prefix
                  . $vcfparser_analysis_type
                  . $file_suffix;
            }
            else {

                ## Write to outputstream
                $genmod_outfile_path =
                  catfile( dirname( devnull() ), "stdout" );    #OutFile
            }
            Program::Variantcalling::Genmod::annotate(
                {
                    infile_path     => $genmod_indata,
                    outfile_path    => $genmod_outfile_path,
                    stderrfile_path => $genmod_xargs_file_path_prefix
                      . $genmod_module
                      . ".stderr.txt",
                    verbosity           => "v",
                    temp_directory_path => $$temp_directory_ref,
                    annotate_region =>
                      $active_parameter_href->{sv_genmod_annotate_regions},
                    FILEHANDLE => $XARGSFILEHANDLE,
                }
            );
            if (
                ( defined( $parameter_href->{dynamic_parameter}{unaffected} ) )
                && ( @{ $parameter_href->{dynamic_parameter}{unaffected} } eq
                    @{ $active_parameter_href->{sample_ids} } )
              )
            {    #Only unaffected - do nothing
            }
            else {

                ## Write to outputstream
                print $XARGSFILEHANDLE "| ";    #Pipe

                ## Get parameters
                $genmod_indata = "-";           #Preparation for next module

                ## Genmod Models
                $genmod_module .= "_models";

                my $use_vep;
                if (
                    (
                        $active_parameter_href->{psv_varianteffectpredictor} > 0
                    )
                    && ( !$active_parameter_href->{sv_genmod_annotate_regions} )
                  )
                {    #Use VEP annotations in compound models

                    $use_vep = 1;
                }
                models(
                    {
                        infile_path => $genmod_indata,
                        outfile_path =>
                          catfile( dirname( devnull() ), "stdout" ),
                        stderrfile_path => $genmod_xargs_file_path_prefix
                          . $genmod_module
                          . ".stderr.txt",
                        family_file => $family_file,
                        family_type => $active_parameter_href
                          ->{sv_genmod_models_family_type},
                        reduced_penetrance_file_path => $active_parameter_href
                          ->{sv_genmod_models_reduced_penetrance_file},
                        thread_number => 4,
                        vep           => $use_vep,
                        whole_gene =>
                          $active_parameter_href->{sv_genmod_models_whole_gene},
                        verbosity           => "v",
                        temp_directory_path => $$temp_directory_ref,
                        FILEHANDLE          => $XARGSFILEHANDLE,
                    }
                );
                print $XARGSFILEHANDLE "| ";    #Pipe

                ## Get parameters
                $genmod_indata = "-";           #Preparation for next module

                ## Genmod Score
                $genmod_module .= "_score";

                score(
                    {
                        infile_path => $genmod_indata,
                        outfile_path =>
                          catfile( dirname( devnull() ), "stdout" ),
                        stderrfile_path => $genmod_xargs_file_path_prefix
                          . $genmod_module
                          . ".stderr.txt",
                        verbosity   => "v",
                        family_file => $family_file,
                        family_type => $active_parameter_href
                          ->{sv_genmod_models_family_type},
                        rank_result => 1,
                        rank_model_file_path =>
                          $active_parameter_href->{sv_rank_model_file},
                        FILEHANDLE => $XARGSFILEHANDLE,
                    }
                );
                print $XARGSFILEHANDLE "| ";    #Pipe

                ##Genmod Compound
                $genmod_module .= "_compound";

                compound(
                    {
                        infile_path  => $genmod_indata,
                        outfile_path => $genmod_outfile_path_prefix
                          . $vcfparser_analysis_type
                          . $file_suffix,
                        stderrfile_path => $genmod_xargs_file_path_prefix
                          . $genmod_module
                          . ".stderr.txt",
                        verbosity           => "v",
                        temp_directory_path => $$temp_directory_ref,
                        vep                 => $use_vep,
                        FILEHANDLE          => $XARGSFILEHANDLE,
                    }
                );
            }
            say {$XARGSFILEHANDLE} "\n";

            # Update endings with contig info
            if (   ( $consensus_analysis_type eq "wes" )
                || ( $consensus_analysis_type eq "rapid" ) )
            {

# Only perform once for exome samples to avoid risking contigs lacking variants throwing errors
                last;
            }
        }

        if (   ( $consensus_analysis_type eq "wgs" )
            || ( $consensus_analysis_type eq "mixed" ) )
        {

            ## Copies file from temporary directory.
            say {$FILEHANDLE} "## Copy file(s) from temporary directory";
            ($xargs_file_counter) = xargs_migrate_contig_files(
                {
                    FILEHANDLE         => $FILEHANDLE,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    contigs_ref        => \@contigs,
                    file_path          => $file_path,
                    program_info_path  => $program_info_path,
                    core_number        => $core_number,
                    xargs_file_counter => $xargs_file_counter,
                    outfile            => $outfile_prefix,
                    file_ending        => $vcfparser_analysis_type
                      . $file_suffix . "*",
                    outdirectory   => $outfamily_directory,
                    temp_directory => $$temp_directory_ref,
                }
            );
        }
        else {

            ## Copies file from temporary directory.
            say {$FILEHANDLE} q{## Copy file from temporary directory};
            migrate_file(
                {
                    infile_path => catfile(
                        $$temp_directory_ref,
                        $outfile_prefix
                          . $vcfparser_analysis_type
                          . $file_suffix . q{*}
                    ),
                    outfile_path => $outfamily_directory,
                    FILEHANDLE   => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} q{wait}, "\n";

            ## Adds the most complete vcf file to sample_info
            add_most_complete_vcf(
                {
                    active_parameter_href => $active_parameter_href,
                    sample_info_href      => $sample_info_href,
                    program_name          => $program_name,
                    path                  => catfile(
                        $outfamily_directory,
                        $outfile_prefix
                          . $vcfparser_analysis_type
                          . $file_suffix
                    ),
                    vcfparser_outfile_counter => $vcfparser_outfile_counter,
                    vcf_file_key => "sv_" . substr( $file_suffix, 1 ) . "_file",
                }
            );
        }
    }
    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        if ( defined( $active_parameter_href->{sv_rank_model_file} ) )
        {    #Add to Sample_info

            my $sv_rank_model_version;
            if ( $active_parameter_href->{sv_rank_model_file} =~
                /v(\d+\.\d+.\d+|\d+\.\d+)/ )
            {

                $sv_rank_model_version = $1;
            }
            add_program_metafile_to_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    program_name     => q{sv_genmod},
                    metafile_tag     => q{sv_rank_model},
                    file =>
                      basename( $active_parameter_href->{sv_rank_model_file} ),
                    path    => $active_parameter_href->{sv_rank_model_file},
                    version => $sv_rank_model_version,
                }
            );

        }
        my $qc_sv_genmod_outfile =
            $$family_id_ref
          . $outfile_tag
          . $call_type
          . $vcfparser_analysis_type
          . $file_suffix;
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => q{sv_genmod},
                outdirectory     => $outfamily_directory,
                outfile          => $qc_sv_genmod_outfile,
            }
        );

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $$family_id_ref,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
}

sub sv_vcfparser {

##sv_vcfparser

##Function : sv_vcfparser performs parsing of varianteffectpredictor annotated variants
##Returns  : "|$xargs_file_counter"
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, $program_info_path, $FILEHANDLE, family_id_ref, $temp_directory_ref, $outaligner_dir_ref, $call_type, $xargs_file_counter
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $FILEHANDLE                 => Sbatch filehandle to write to
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type
##         : $xargs_file_counter         => The xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
    my $call_type;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => "SV", strict_type => 1, store => \$call_type, },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use MIP::Get::File qw{get_file_suffix};
    use MIP::Delete::List qw{ delete_contig_elements delete_male_contig};
    use MIP::Recipes::Xargs qw{ xargs_command };
    use MIP::IO::Files qw(migrate_file xargs_migrate_contig_files);
    use Program::Variantcalling::Mip qw(vcfparser);
    use MIP::QC::Record qw(add_program_outfile_to_sample_info);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_family);

    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $core_number =
      $active_parameter_href->{module_core_number}{ "p" . $program_name };
    my $xargs_file_path_prefix;
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $FILEHANDLE      = IO::Handle->new();    #Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory     => catfile( lc($$outaligner_dir_ref) ),
            call_type             => $call_type,
            core_number           => $core_number,
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
            temp_directory => $$temp_directory_ref,
        }
    );

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    $parameter_href->{ "p" . $program_name }{indirectory} =
      $outfamily_directory;    #Used downstream

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$$family_id_ref}{psv_varianteffectpredictor}{file_tag};
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{ "p" . $program_name }{file_tag};
    my $infile_prefix       = $$family_id_ref . $infile_tag . $call_type;
    my $outfile_prefix      = $$family_id_ref . $outfile_tag . $call_type;
    my $outfile_path_prefix = catfile( $$temp_directory_ref, $outfile_prefix );

    ## Assign suffix
    my $file_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            jobid_chain    => $job_id_chain,
        }
    );

    ## Removes an element from array and return new array while leaving orginal elements_ref untouched
    my @contigs = delete_contig_elements(
        {
            elements_ref       => \@{ $file_info_href->{contigs_size_ordered} },
            remove_contigs_ref => [qw{ MT M }],
        }
    );

    ### If no males or other remove contig Y from all downstream analysis
    ## Removes contig_names from contigs array if no male or other found
    @contigs = delete_male_contig(
        {
            contigs_ref => \@contigs,
            found_male  => $active_parameter_href->{found_male},
        }
    );

    if (   ( $consensus_analysis_type eq "wgs" )
        || ( $consensus_analysis_type eq "mixed" ) )
    {    #Transfer contig files

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        ($xargs_file_counter) = xargs_migrate_contig_files(
            {
                FILEHANDLE         => $FILEHANDLE,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                contigs_ref        => \@contigs,
                file_path          => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
                xargs_file_counter => $xargs_file_counter,
                infile             => $infile_prefix,
                indirectory        => $infamily_directory,
                temp_directory     => $$temp_directory_ref,
            }
        );
    }
    else {

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        migrate_file(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => catfile(
                    $infamily_directory,
                    $$family_id_ref . $infile_tag . $call_type . $file_suffix
                ),
                outfile_path => $$temp_directory_ref
            }
        );
        say {$FILEHANDLE} q{wait}, "\n";
    }

    ## vcfparser
    say {$FILEHANDLE} "## vcfparser";

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    foreach my $contig (@contigs) {

        ## Get parameters
        my $vcfparser_infile_prefix          = $infile_prefix;
        my $vcfparser_outfile_prefix         = $outfile_prefix;
        my $vcfparser_xargs_file_path_prefix = $xargs_file_path_prefix;

        if (   ( $consensus_analysis_type eq "wgs" )
            || ( $consensus_analysis_type eq "mixed" ) )
        {    #Update endings with contig info

            $vcfparser_infile_prefix  = $infile_prefix . "_" . $contig;
            $vcfparser_outfile_prefix = $outfile_prefix . "_" . $contig;
            $vcfparser_xargs_file_path_prefix =
              $xargs_file_path_prefix . "." . $contig;
        }

        my $padding;
        if ( $contig =~ /MT|M/ ) {

            $padding = 10;    #Special case for mitochondrial contig annotation
        }

        my @select_feature_annotation_columns;
        my $select_file;
        my $select_file_matching_column;
        my $select_outfile;
        if ( $active_parameter_href->{sv_vcfparser_select_file} ) {

            if (
                !check_entry_hash_of_array(
                    {
                        hash_ref => $file_info_href,
                        key      => "select_file_contigs",
                        element  => $contig,
                    }
                )
              )
            {

                $select_file =
                  catfile( $active_parameter_href->{sv_vcfparser_select_file} )
                  ;    #List of genes to analyse separately
                $select_file_matching_column = $active_parameter_href
                  ->{sv_vcfparser_select_file_matching_column}
                  ;    #Column of HGNC Symbol in SelectFile (-sf)

                if (
                    (
                        $active_parameter_href
                        ->{sv_vcfparser_select_feature_annotation_columns}
                    )
                    && (
                        @{
                            $active_parameter_href
                              ->{sv_vcfparser_select_feature_annotation_columns}
                        }
                    )
                  )
                {

                    @select_feature_annotation_columns =
                      @{ $active_parameter_href
                          ->{sv_vcfparser_select_feature_annotation_columns} };
                }
                $select_outfile = catfile( $$temp_directory_ref,
                    $vcfparser_outfile_prefix . ".selected" . $file_suffix );
            }
        }

        vcfparser(
            {
                range_feature_annotation_columns_ref => \@{
                    $active_parameter_href
                      ->{sv_vcfparser_range_feature_annotation_columns}
                },
                select_feature_annotation_columns_ref =>
                  \@select_feature_annotation_columns,
                infile_path => catfile(
                    $$temp_directory_ref,
                    $vcfparser_infile_prefix . $file_suffix
                ),
                outfile_path => catfile(
                    $$temp_directory_ref,
                    $vcfparser_outfile_prefix . $file_suffix
                ),
                stderrfile_path => $vcfparser_xargs_file_path_prefix
                  . ".stderr.txt",
                range_feature_file_path =>
                  $active_parameter_href->{sv_vcfparser_range_feature_file},
                select_feature_file_path       => $select_file,
                select_feature_matching_column => $select_file_matching_column,
                select_outfile                 => $select_outfile,
                parse_vep =>
                  $active_parameter_href->{psv_varianteffectpredictor},
                per_gene   => $active_parameter_href->{sv_vcfparser_per_gene},
                padding    => $padding,
                FILEHANDLE => $XARGSFILEHANDLE,
            }
        );
        say {$XARGSFILEHANDLE} "\n";

        if (   ( $consensus_analysis_type eq "wes" )
            || ( $consensus_analysis_type eq "rapid" ) )
        {    #Update endings with contig info

            last
              ; #Only perform once for exome samples to avoid risking contigs lacking variants throwing errors
        }
    }

    if (   ( $consensus_analysis_type eq "wgs" )
        || ( $consensus_analysis_type eq "mixed" ) )
    {

        ## QC Data File(s)
        migrate_file(
            {
                infile_path => catfile(
                    $$temp_directory_ref,
                    $outfile_prefix . q{_} . $contigs[0] . $file_suffix
                ),    #Add contig info
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} q{wait}, "\n";
    }

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        ## Clear old vcfparser entry if present
        if ( exists( $sample_info_href->{$program_name} ) ) {

            delete( $sample_info_href->{$program_name} );
        }

        my $outfile_sample_info_prefix = $outfile_prefix;

        if (   ( $consensus_analysis_type eq "wgs" )
            || ( $consensus_analysis_type eq "mixed" ) )
        {    #Update endings with contig info

            $outfile_sample_info_prefix .= "_" . $contigs[0];
        }

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => $program_name,
                outdirectory     => $outfamily_directory,
                outfile          => $outfile_sample_info_prefix . $file_suffix,
            }
        );

        my %gene_panels = (
            range_file  => "sv_vcfparser_range_feature_file",
            select_file => "sv_vcfparser_select_file",
        );
        while ( my ( $gene_panel_key, $gene_panel_file ) = each(%gene_panels) )
        {

            ## Collect databases(s) from a potentially merged gene panel file and adds them to sample_info
            collect_gene_panels(
                {
                    sample_info_href => $sample_info_href,
                    family_id_ref    => $family_id_ref,
                    program_name_ref => \$program_name,
                    aggregate_gene_panel_file =>
                      $active_parameter_href->{$gene_panel_file},
                    aggregate_gene_panels_key => $gene_panel_key,
                }
            );
        }
    }

    close $XARGSFILEHANDLE;

    my $vcfparser_analysis_type = "";

    for (
        my $vcfparser_outfile_counter = 0 ;
        $vcfparser_outfile_counter <
        $active_parameter_href->{vcfparser_outfile_count} ;
        $vcfparser_outfile_counter++
      )
    {

        if ( $vcfparser_outfile_counter == 1 ) {

            $vcfparser_analysis_type = ".selected";    #Select file variants

            ## Removes an element from array and return new array while leaving orginal elements_ref untouched
            @contigs = delete_contig_elements(
                {
                    elements_ref =>
                      \@{ $file_info_href->{sorted_select_file_contigs} },
                    remove_contigs_ref => [qw{ MT M }],
                }
            );
        }

        if (   ( $consensus_analysis_type eq "wgs" )
            || ( $consensus_analysis_type eq "mixed" ) )
        {

            ## Copies file from temporary directory.
            say {$FILEHANDLE} "## Copy file(s) from temporary directory";
            ($xargs_file_counter) = xargs_migrate_contig_files(
                {
                    FILEHANDLE         => $FILEHANDLE,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    contigs_ref        => \@contigs,
                    file_path          => $file_path,
                    program_info_path  => $program_info_path,
                    core_number        => $core_number,
                    xargs_file_counter => $xargs_file_counter,
                    outfile            => $outfile_prefix,
                    file_ending        => $vcfparser_analysis_type
                      . $file_suffix . "*",
                    outdirectory   => $outfamily_directory,
                    temp_directory => $$temp_directory_ref,
                }
            );
        }
        else {

            ## Copies file from temporary directory.
            say {$FILEHANDLE} q{## Copy file from temporary directory};
            migrate_file(
                {
                    infile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $file_suffix . q{*},
                    outfile_path => $outfamily_directory,
                    FILEHANDLE   => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} q{wait}, "\n";

            ## Adds the most complete vcf file to sample_info
            add_most_complete_vcf(
                {
                    active_parameter_href => $active_parameter_href,
                    sample_info_href      => $sample_info_href,
                    program_name          => $program_name,
                    path                  => catfile(
                        $outfamily_directory,
                        $outfile_prefix
                          . $vcfparser_analysis_type
                          . $file_suffix
                    ),
                    vcfparser_outfile_counter => $vcfparser_outfile_counter,
                    vcf_file_key => "sv_" . substr( $file_suffix, 1 ) . "_file",
                }
            );
        }
    }
    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $$family_id_ref,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
}

sub sv_combinevariantcallsets {

##sv_combinevariantcallsets

##Function : CombineVariants to combine all structural variants call from different callers.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, family_id_ref, $temp_directory_ref, $reference_dir_ref, $outaligner_dir_ref, $call_type
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $reference_dir_ref          => MIP reference directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $reference_dir_ref;
    my $outaligner_dir_ref;
    my $call_type;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        reference_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{reference_dir},
            strict_type => 1,
            store       => \$reference_dir_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => "SV", strict_type => 1, store => \$call_type, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use MIP::IO::Files qw(migrate_file);
    use MIP::Get::File qw{get_file_suffix get_merged_infile_prefix };
    use MIP::Set::File qw{set_file_suffix};
    use Program::Variantcalling::Svdb qw(merge query);
    use MIP::Program::Variantcalling::Bcftools
      qw (bcftools_merge bcftools_view bcftools_annotate);
    use Program::Htslib qw(bgzip tabix);
    use Program::Variantcalling::Vt qw(decompose);
    use Program::Variantcalling::Genmod qw(annotate);
    use Program::Variantcalling::Vcfanno qw(vcfanno);
    use MIP::QC::Record qw(add_program_outfile_to_sample_info);

    my @structural_variant_callers;    #Stores callers that have been executed
    my @parallel_chains
      ;    #Stores the parallel chains that jobIds should be inherited from
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory     => catfile( lc($$outaligner_dir_ref) ),
            call_type             => $call_type,
            core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
            temp_directory => $$temp_directory_ref,
        }
    );
    my ( $volume, $directory, $stderr_file ) =
      File::Spec->splitpath( $program_info_path . ".stderr.txt" )
      ;    #Split to enable submission to &sample_info_qc later

    ## Assign directories
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    $parameter_href->{ "p" . $program_name }{indirectory} =
      $outfamily_directory;    #Used downstream

    ## Assign file_tags
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{ "p" . $program_name }{file_tag};
    my $merged_file_path_prefix =
      catfile( $$temp_directory_ref, $$family_id_ref . "_" . $call_type );
    my %file_path_prefix;
    my $outfile_path_prefix = catfile( $$temp_directory_ref,
        $$family_id_ref . $outfile_tag . $call_type );

    ### Assign suffix
    my %suffix;

    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    ## Collect infiles for all sample_ids for programs that do not do joint calling to enable migration to temporary directory
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Add merged infile name prefix after merging all BAM files per sample_id
        my $merged_infile_prefix = get_merged_infile_prefix(
            {
                file_info_href => $file_info_href,
                sample_id      => $sample_id,
            }
        );

        foreach my $structural_variant_caller (
            @{
                $parameter_href->{dynamic_parameter}{structural_variant_callers}
            }
          )
        {

            if (
                ( $active_parameter_href->{$structural_variant_caller} > 0 )
                && ( $structural_variant_caller !~
                    /pmanta|pdelly_reformat|ptiddit/ )
              )
            { #Expect vcf. Special case: manta, delly and tiddit are processed by joint calling and per family

                ## Assign directories
                my $program_outdirectory_name =
                  $parameter_href->{$structural_variant_caller}{outdir_name};
                my $insample_directory =
                  catdir( $active_parameter_href->{outdata_dir},
                    $sample_id, $$outaligner_dir_ref,
                    $program_outdirectory_name );

                ## Assign file_tags
                my $infile_tag =
                  $file_info_href->{$sample_id}{$structural_variant_caller}
                  {file_tag};
                my $infile_prefix = $merged_infile_prefix . $infile_tag;
                $file_path_prefix{$sample_id}{$structural_variant_caller} =
                  catfile( $$temp_directory_ref, $infile_prefix );

                ## Assign suffix
                $suffix{$structural_variant_caller} = get_file_suffix(
                    {
                        parameter_href => $parameter_href,
                        suffix_key     => "outfile_suffix",
                        program_name   => $structural_variant_caller,
                    }
                );

                if (
                    !(
                        any {
                            $_ eq $parameter_href->{$structural_variant_caller}
                              {chain}
                        }
                        @parallel_chains
                    )
                  )
                {    #If element is not part of array

                    push( @parallel_chains,
                        $parameter_href->{$structural_variant_caller}{chain} );
                }

                ## Copy file(s) to temporary directory
                say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
                migrate_file(
                    {
                        FILEHANDLE  => $FILEHANDLE,
                        infile_path => catfile(
                            $insample_directory,
                            $infile_prefix
                              . $suffix{$structural_variant_caller} . q{*}
                        ),
                        outfile_path => $$temp_directory_ref
                    }
                );

                say {$FILEHANDLE} q{wait}, "\n";

                ## Reformat variant calling file and index
                view_vcf(
                    {
                        infile_path => $file_path_prefix{$sample_id}
                          {$structural_variant_caller}
                          . $suffix{$structural_variant_caller},
                        outfile_path_prefix => $file_path_prefix{$sample_id}
                          {$structural_variant_caller},
                        output_type => "z",
                        FILEHANDLE  => $FILEHANDLE,
                    }
                );
            }
        }
    }

    ## Merged sample files to one family file (samples > 1) else reformat to standardise
    foreach my $structural_variant_caller (
        @{ $parameter_href->{dynamic_parameter}{structural_variant_callers} } )
    {

        if ( $active_parameter_href->{$structural_variant_caller} > 0
            && (
                $structural_variant_caller !~ /pmanta|pdelly_reformat|ptiddit/ )
          )
        { #Expect vcf. Special case: manta is processed by joint calling and per family

            ## Assemble file paths by adding file ending
            my @file_paths = map {
                    $file_path_prefix{$_}{$structural_variant_caller}
                  . $suffix{$structural_variant_caller} . ".gz"
            } @{ $active_parameter_href->{sample_ids} };

            if ( scalar( @{ $active_parameter_href->{sample_ids} } ) > 1 ) {

                ## Merge all structural variant caller's vcf files per sample_id
                say {$FILEHANDLE}
"## Merge all structural variant caller's vcf files per sample_id";

                bcftools_merge(
                    {
                        infile_paths_ref => \@file_paths,
                        outfile_path     => catfile(
                            $$temp_directory_ref,
                            $$family_id_ref . "_"
                              . $structural_variant_caller
                              . $outfile_suffix
                        ),
                        output_type     => "v",
                        stderrfile_path => $program_info_path . "_"
                          . $structural_variant_caller
                          . "_merge.stderr.txt",
                        FILEHANDLE => $FILEHANDLE,
                    }
                );
                say {$FILEHANDLE} "\n";
            }
            else {

                ## Reformat all structural variant caller's vcf files per sample_id
                say {$FILEHANDLE}
"## Reformat all structural variant caller's vcf files per sample_id";

                bcftools_view(
                    {
                        infile_path  => $file_paths[0],    #Can be only one
                        outfile_path => catfile(
                            $$temp_directory_ref,
                            $$family_id_ref . "_"
                              . $structural_variant_caller
                              . $outfile_suffix
                        ),
                        output_type     => "v",
                        stderrfile_path => $program_info_path . "_"
                          . $structural_variant_caller
                          . "_merge.stderr.txt",
                        FILEHANDLE => $FILEHANDLE,
                    }
                );
                say {$FILEHANDLE} "\n";
            }
        }
    }

    ## Migrate joint calling per family callers like Manta and Delly
    foreach my $structural_variant_caller (
        @{ $parameter_href->{dynamic_parameter}{structural_variant_callers} } )
    {

        if ( $active_parameter_href->{$structural_variant_caller} > 0
            && (
                $structural_variant_caller =~ /pmanta|pdelly_reformat|ptiddit/ )
          )
        { #Expect vcf. Special case: manta, delly, tiddit are processed by joint calling and per family

            ## Assign directories
            my $program_outdirectory_name =
              $parameter_href->{$structural_variant_caller}{outdir_name};
            my $infamily_directory =
              catfile( $active_parameter_href->{outdata_dir},
                $$family_id_ref, $$outaligner_dir_ref,
                $program_outdirectory_name );

            ## Assign file_tags
            my $infile_tag =
              $file_info_href->{$$family_id_ref}{$structural_variant_caller}
              {file_tag};
            my $infile_prefix =
              $$family_id_ref . $infile_tag . "_" . $call_type;

            ## Assign suffix
            my $infile_suffix = get_file_suffix(
                {
                    parameter_href => $parameter_href,
                    suffix_key     => "outfile_suffix",
                    program_name   => $structural_variant_caller,
                }
            );

            if (
                !(
                    any {
                        $_ eq
                          $parameter_href->{$structural_variant_caller}{chain}
                    }
                    @parallel_chains
                )
              )
            {    #If element is not part of array

                push( @parallel_chains,
                    $parameter_href->{$structural_variant_caller}{chain} );
            }

            ## Copy file(s) to temporary directory
            say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
            migrate_file(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => catfile(
                        $infamily_directory,
                        $infile_prefix . $infile_suffix . q{*}
                    ),
                    outfile_path => $$temp_directory_ref
                }
            );
            say {$FILEHANDLE} q{wait}, "\n";

            if ( $active_parameter_href->{sv_vt_decompose} > 0 ) {

                ## Split multiallelic variants
                say {$FILEHANDLE} "## Split multiallelic variants";
                decompose(
                    {
                        infile_path => catfile(
                            $$temp_directory_ref,
                            $infile_prefix . $infile_suffix
                        ),
                        outfile_path => catfile(
                            $$temp_directory_ref,
                            $$family_id_ref . "_"
                              . $structural_variant_caller
                              . $infile_suffix
                        ),
                        FILEHANDLE          => $FILEHANDLE,
                        smart_decomposition => 1,
                    }
                );
                say {$FILEHANDLE} "\n";
            }
        }
    }

    ## Merge structural variant caller's family vcf files
    say {$FILEHANDLE} "## Merge structural variant caller's family vcf files";

    ## Get parameters
    my @infile_paths;
    foreach my $structural_variant_caller (
        @{ $parameter_href->{dynamic_parameter}{structural_variant_callers} } )
    {

        if ( $active_parameter_href->{$structural_variant_caller} > 0 )
        {    #Expect vcf

            my $variant_caller_alias =
              $parameter_href->{$structural_variant_caller}{outdir_name};
            push(
                @infile_paths,
                catfile(
                    $$temp_directory_ref,
                    $$family_id_ref . "_"
                      . $structural_variant_caller
                      . $outfile_suffix . ":"
                      . $variant_caller_alias
                )
            );
        }
    }

    Program::Variantcalling::Svdb::merge(
        {
            infile_paths_ref => \@infile_paths,
            outfile_path     => $merged_file_path_prefix . $outfile_suffix,
            priority   => $active_parameter_href->{sv_svdb_merge_prioritize},
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

    my $alt_file_tag = "";    #Alternative file tag

    if ( $active_parameter_href->{sv_vt_decompose} > 0 ) {

        $alt_file_tag = "_vt";    #Update file tag

        ## Split multiallelic variants
        say {$FILEHANDLE} "## Split multiallelic variants";
        decompose(
            {
                infile_path  => $merged_file_path_prefix . $outfile_suffix,
                outfile_path => $merged_file_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix,
                FILEHANDLE          => $FILEHANDLE,
                smart_decomposition => 1,
            }
        );
        say {$FILEHANDLE} "\n";
    }
    if ( $active_parameter_href->{sv_svdb_query} > 0 ) {

        use MIP::Gnu::Coreutils qw(gnu_mv);

        my $infile_path =
          $merged_file_path_prefix . $alt_file_tag . $outfile_suffix;
        $alt_file_tag .= "_svdbq";    #Update alternative ending

        my $annotation_file_counter = 0;    #Ensure correct infile
        my $outfile_tracker         = 0;    #Ensure correct outfiles

        while ( my ( $query_db_file, $query_db_tag ) =
            each( %{ $active_parameter_href->{sv_svdb_query_db_files} } ) )
        {

            if ($annotation_file_counter) {

                $infile_path =
                    $merged_file_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix . "."
                  . $outfile_tracker;
                $outfile_tracker++;    #Increment now that infile has been set
            }
            query(
                {
                    infile_path  => $infile_path,
                    outfile_path => $merged_file_path_prefix
                      . $alt_file_tag
                      . $outfile_suffix . "."
                      . $outfile_tracker,
                    dbfile_path   => $query_db_file,
                    bnd_distance  => 25000,
                    overlap       => 0.8,
                    hit_tag       => $query_db_tag,
                    frequency_tag => $query_db_tag . "AF",
                    FILEHANDLE    => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} "\n";
            $annotation_file_counter++;
        }

        ## Rename to remove outfile_tracker
        gnu_mv(
            {
                infile_path => $merged_file_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix . "."
                  . $outfile_tracker,
                outfile_path => $merged_file_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix,
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} "\n";
    }

    my $outfile_alt_file_tag .= $alt_file_tag . "_sorted"; #Alternative file tag

    ## Writes sbatch code to supplied filehandle to sort variants in vcf format
    sort_vcf(
        {
            active_parameter_href => $active_parameter_href,
            FILEHANDLE            => $FILEHANDLE,
            sequence_dict_file    => catfile(
                $$reference_dir_ref,
                $file_info_href->{human_genome_reference_name_prefix} . ".dict"
            ),
            infile_paths_ref =>
              [ $merged_file_path_prefix . $alt_file_tag . $outfile_suffix ],
            outfile => $outfile_path_prefix
              . $outfile_alt_file_tag
              . $outfile_suffix,
        }
    );
    print {$FILEHANDLE} "\n";

    $alt_file_tag = $outfile_alt_file_tag;

    ## Remove FILTER ne PASS
    if ( $active_parameter_href->{sv_bcftools_view_filter} > 0 ) {

        say {$FILEHANDLE} "## Remove FILTER ne PASS";
        bcftools_view(
            {
                apply_filters_ref => ["PASS"],
                infile_path       => $outfile_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix,
                outfile_path => $outfile_path_prefix
                  . $alt_file_tag . "_filt"
                  . $outfile_suffix,
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} "\n";

        $alt_file_tag .= "_filt";    #Update file tag
    }

    ## Remove common variants
    if ( $active_parameter_href->{sv_genmod_filter} > 0 ) {

        say {$FILEHANDLE} "## Remove common variants";
        Program::Variantcalling::Genmod::annotate(
            {
                infile_path => $outfile_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix,
                outfile_path => catfile( dirname( devnull() ), "stdout" ),
                verbosity    => "v",
                temp_directory_path => $$temp_directory_ref,
                thousand_g_file_path =>
                  $active_parameter_href->{sv_genmod_filter_1000g},
                FILEHANDLE => $FILEHANDLE,
            }
        );
        print {$FILEHANDLE} "| ";

        $alt_file_tag .= "_genmod_filter";    #Update file tag

        Program::Variantcalling::Genmod::filter(
            {
                infile_path  => "-",
                outfile_path => $outfile_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix,
                verbosity => "v",
                threshold =>
                  $active_parameter_href->{sv_genmod_filter_threshold},
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} "\n";
    }

    ## Annotate 1000G structural variants
    if ( $active_parameter_href->{sv_vcfanno} > 0 ) {

        say {$FILEHANDLE} "## Annotate 1000G structural variants";
        vcfanno(
            {
                infile_path => $outfile_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix,
                toml_configfile_path =>
                  $active_parameter_href->{sv_vcfanno_config},
                lua        => $active_parameter_href->{sv_vcfanno_lua},
                ends       => 1,
                FILEHANDLE => $FILEHANDLE,
            }
        );
        print {$FILEHANDLE} "| ";
        print {$FILEHANDLE}
q?perl -nae 'if($_=~/^#/) {print $_} else {$F[7]=~s/\[||\]//g; print join("\t", @F), "\n"}' ?
          ;    #Remove "[" and "]" from INFO as it breaks vcf format

        $alt_file_tag .= "_vcfanno";    #Update file tag

        say {$FILEHANDLE} "> "
          . $outfile_path_prefix
          . $alt_file_tag
          . $outfile_suffix, "\n";

        if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

            add_program_outfile_to_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    program_name     => q{sv_combinevariantcallsets},
                    outdirectory     => $directory,
                    outfile          => $stderr_file,
                }
            );
        }

        say {$FILEHANDLE}
          "## Add header for 1000G annotation of structural variants";
        bcftools_annotate(
            {
                infile_path => $outfile_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix,
                outfile_path => $outfile_path_prefix
                  . $alt_file_tag
                  . "_bcftools_annotate"
                  . $outfile_suffix,
                output_type => "v",
                headerfile_path =>
                  $active_parameter_href->{sv_vcfannotation_header_lines_file},
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} "\n";
        $alt_file_tag .= "_bcftools_annotate";    #Update file tag
    }

    if ( $alt_file_tag ne "" ) {    #Then we have something to rename

        ## Writes sbatch code to supplied filehandle to sort variants in vcf format
        sort_vcf(
            {
                active_parameter_href => $active_parameter_href,
                FILEHANDLE            => $FILEHANDLE,
                sequence_dict_file    => catfile(
                    $$reference_dir_ref,
                    $file_info_href->{human_genome_reference_name_prefix}
                      . ".dict"
                ),
                infile_paths_ref =>
                  [ $outfile_path_prefix . $alt_file_tag . $outfile_suffix ],
                outfile => $outfile_path_prefix . $outfile_suffix,
            }
        );

        print {$FILEHANDLE} "\n";
    }

    if ( $active_parameter_href->{sv_combinevariantcallsets_bcf_file} ) {

        ## Reformat variant calling file and index
        view_vcf(
            {
                infile_path         => $outfile_path_prefix . $outfile_suffix,
                outfile_path_prefix => $outfile_path_prefix,
                output_type         => "b",
                FILEHANDLE          => $FILEHANDLE,
            }
        );

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                infile_path  => $outfile_path_prefix . q{.bcf*},
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
    }

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            infile_path  => $outfile_path_prefix . $outfile_suffix,
            outfile_path => $outfamily_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, "\n";

    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        my $qc_svdb_outfile =
          $$family_id_ref . $outfile_tag . $call_type . $outfile_suffix;
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => 'svdb',
                outdirectory     => $outfamily_directory,
                outfile          => $qc_svdb_outfile,
            }
        );

        $sample_info_href->{sv_vcf_file}{ready_vcf}{path} =
          catfile( $outfamily_directory,
            $$family_id_ref . $outfile_tag . $call_type . $outfile_suffix );

        if ( $active_parameter_href->{sv_combinevariantcallsets_bcf_file} ) {

            $sample_info_href->{sv_bcf_file}{path} =
              catfile( $outfamily_directory,
                $$family_id_ref . $outfile_tag . $call_type . ".bcf" );
        }

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref => \@{ $active_parameter_href->{sample_ids} },
                parallel_chains_ref => \@parallel_chains,
                family_id           => $$family_id_ref,
                path                => $job_id_chain,
                log                 => $log,
                sbatch_file_name    => $file_path,
            }
        );
    }
}

sub delly_reformat {

##delly_reformat

##Function : Merge, regenotype, and filter using delly
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, $program_info_path, family_id_ref, $temp_directory_ref, $reference_dir_ref, $outaligner_dir_ref, $xargs_file_counter, $call_type
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $reference_dir_ref          => MIP reference directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $xargs_file_counter         => The xargs file counter
##         : $call_type                  => Variant call type

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $reference_dir_ref;
    my $outaligner_dir_ref;
    my $xargs_file_counter;
    my $call_type;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        reference_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{reference_dir},
            strict_type => 1,
            store       => \$reference_dir_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
        call_type =>
          { default => "SV", strict_type => 1, store => \$call_type, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use MIP::IO::Files qw(migrate_file xargs_migrate_contig_files);
    use MIP::Delete::List qw{ delete_contig_elements };
    use MIP::Get::File qw{get_file_suffix get_merged_infile_prefix };
    use MIP::Set::File qw{set_file_suffix};
    use MIP::Recipes::Xargs qw{ xargs_command };
    use MIP::Gnu::Coreutils qw(gnu_mv);
    use Program::Variantcalling::Delly qw(call merge filter);
    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_merge bcftools_index bcftools_concat };
    use MIP::QC::Record qw(add_program_outfile_to_sample_info);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_family);

    my $core_number =
      $active_parameter_href->{module_core_number}{ "p" . $program_name };
    my $program_outdirectory_name =
      $parameter_href->{ "p" . $program_name }{outdir_name};
    my $xargs_file_path_prefix;
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory     => catfile(
                lc($$outaligner_dir_ref),
                lc($program_outdirectory_name)
            ),
            core_number => $core_number,
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
            temp_directory => $$temp_directory_ref
        }
    );

    ## Assign directories
    my $outfamily_directory = catfile(
        $active_parameter_href->{outdata_dir},
        $$family_id_ref,
        lc($$outaligner_dir_ref),
        lc($program_outdirectory_name)
    );
    $parameter_href->{ "p" . $program_name }{indirectory} =
      $outfamily_directory;    #Used downstream

    ## Assign file_tags
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{ "p" . $program_name }{file_tag};
    my $outfile_prefix = $$family_id_ref . $outfile_tag . "_" . $call_type;
    my $outfile_path_prefix = catfile( $$temp_directory_ref, $outfile_prefix );

    ### Assign suffix
    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    ## Removes an element from array and return new array while leaving orginal elements_ref untouched
    my @contigs = delete_contig_elements(
        {
            elements_ref       => \@{ $file_info_href->{contigs_size_ordered} },
            remove_contigs_ref => [qw{ MT M }],
        }
    );

    ## Removes contigs from supplied contigs_ref
    remove_array_element(
        {
            contigs_ref        => \@contigs,
            remove_contigs_ref => ["Y"]
            , #Skip contig Y throughout since sometimes there are no variants particularly for INS
        }
    );

    ## Collect files and suffix for all sample_ids
    my %infile_path_prefix;
    my %file_path_prefix;
    my %suffix;
    my @program_tag_keys = ( "pgatk_baserecalibration", "pdelly_call" );
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Assign directories
        my $insample_directory_bam =
          catdir( $active_parameter_href->{outdata_dir},
            $sample_id, $$outaligner_dir_ref );
        my $insample_directory_bcf =
          catdir( $active_parameter_href->{outdata_dir},
            $sample_id, $$outaligner_dir_ref,
            $parameter_href->{pdelly_call}{outdir_name} );

        foreach my $infile_tag_key (@program_tag_keys) {

            ## Add merged infile name prefix after merging all BAM files per sample_id
            my $merged_infile_prefix = get_merged_infile_prefix(
                {
                    file_info_href => $file_info_href,
                    sample_id      => $sample_id,
                }
            );

            ## Assign file_tags
            my $infile_tag =
              $file_info_href->{$sample_id}{$infile_tag_key}{file_tag};
            my $infile_prefix = $merged_infile_prefix . $infile_tag;

            $infile_path_prefix{$sample_id}{$infile_tag_key} =
              catfile( $$temp_directory_ref, $infile_prefix );
            $file_path_prefix{$sample_id} =
              catfile( $$temp_directory_ref,
                $merged_infile_prefix . $outfile_tag );    #Used downstream

            if ( $infile_tag_key eq "pdelly_call" ) {      #BCFs

                ## Assign suffix
                $suffix{$infile_tag_key} = get_file_suffix(
                    {
                        parameter_href => $parameter_href,
                        suffix_key     => "outfile_suffix",
                        program_name   => $infile_tag_key,
                    }
                );

              SV_TYPE:
                foreach
                  my $sv_type ( @{ $active_parameter_href->{delly_types} } )
                {

                    my $file_ending = "_"
                      . $sv_type
                      . substr( $suffix{$infile_tag_key}, 0, 2 ) . "*";

                    if ( $sv_type ne "TRA" ) {

                        ## Copy file(s) to temporary directory
                        say {$FILEHANDLE}
                          q{## Copy file(s) to temporary directory};
                        ($xargs_file_counter) = xargs_migrate_contig_files(
                            {
                                FILEHANDLE         => $FILEHANDLE,
                                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                                contigs_ref        => \@contigs,
                                file_path          => $file_path,
                                program_info_path  => $program_info_path,
                                core_number        => $core_number,
                                xargs_file_counter => $xargs_file_counter,
                                infile => $merged_infile_prefix . $infile_tag,
                                indirectory    => $insample_directory_bcf,
                                file_ending    => $file_ending,
                                temp_directory => $$temp_directory_ref,
                            }
                        );
                    }
                    else {

                        say {$FILEHANDLE}
                          q{## Copy file(s) to temporary directory};
                        migrate_file(
                            {
                                FILEHANDLE  => $FILEHANDLE,
                                infile_path => catfile(
                                    $insample_directory_bcf,
                                    $merged_infile_prefix
                                      . $infile_tag
                                      . $file_ending
                                ),
                                outfile_path =>
                                  $active_parameter_href->{temp_directory}
                            }
                        );
                    }
                }
            }
            else {    #BAMs

                ## Assign suffix
                $suffix{$infile_tag_key} = get_file_suffix(
                    {
                        parameter_href => $parameter_href,
                        suffix_key     => q{alignment_file_suffix},
                        jobid_chain =>
                          $parameter_href->{$infile_tag_key}{chain},
                    }
                );

                my $file_ending =
                  substr( $suffix{$infile_tag_key}, 0, 2 ) . "*";

                ## Copy file(s) to temporary directory
                say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
                migrate_file(
                    {
                        FILEHANDLE  => $FILEHANDLE,
                        infile_path => catfile(
                            $insample_directory_bam,
                            $merged_infile_prefix . $infile_tag . $file_ending
                        ),
                        outfile_path => $active_parameter_href->{temp_directory}
                    }
                );

                ## Copy file(s) to temporary directory
                say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
                ($xargs_file_counter) = xargs_migrate_contig_files(
                    {
                        FILEHANDLE        => $FILEHANDLE,
                        XARGSFILEHANDLE   => $XARGSFILEHANDLE,
                        contigs_ref       => \@contigs,
                        file_path         => $file_path,
                        program_info_path => $program_info_path,
                        core_number       => ( $core_number - 1 )
                        ,    #Compensate for cp of entire BAM TRA, see above
                        xargs_file_counter => $xargs_file_counter,
                        infile         => $merged_infile_prefix . $infile_tag,
                        indirectory    => $insample_directory_bam,
                        file_ending    => $file_ending,
                        temp_directory => $$temp_directory_ref,
                    }
                );
            }
            say {$FILEHANDLE} q{wait}, "\n";
        }
    }

    if ( scalar( @{ $active_parameter_href->{sample_ids} } ) > 1 ) {

        ### Delly merge
        say {$FILEHANDLE} "## delly merge \n";

        say {$FILEHANDLE}
          "## Fix locale bug using old centosOS and Boost library";
        say {$FILEHANDLE} q?LC_ALL="C"; export LC_ALL ?, "\n\n";

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                FILEHANDLE         => $FILEHANDLE,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                file_path          => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
                xargs_file_counter => $xargs_file_counter,
            }
        );
      SV_TYPE:
        foreach my $sv_type ( @{ $active_parameter_href->{delly_types} } ) {

            if ( $sv_type ne "TRA" ) {

              CONTIG:
                foreach my $contig (@contigs) {

                    ## Assemble file paths by adding file ending
                    my @file_paths = map {
                            $infile_path_prefix{$_}{pdelly_call} . "_"
                          . $contig . "_"
                          . $sv_type
                          . $suffix{pdelly_call}
                    } @{ $active_parameter_href->{sample_ids} };

                    Program::Variantcalling::Delly::merge(
                        {
                            infile_paths_ref => \@file_paths,
                            outfile_path     => $outfile_path_prefix . "_"
                              . $contig . "_"
                              . $sv_type . ".bcf",
                            stdoutfile_path => $xargs_file_path_prefix . "."
                              . $contig . "."
                              . $sv_type
                              . ".stdout.txt",
                            stderrfile_path => $xargs_file_path_prefix . "."
                              . $contig . "."
                              . $sv_type
                              . ".stderr.txt",
                            sv_type    => $sv_type,
                            min_size   => 0,
                            max_size   => 100000000,
                            FILEHANDLE => $XARGSFILEHANDLE,
                        }
                    );
                    say {$XARGSFILEHANDLE} "\n";
                }
            }
            else {

                ## Assemble file paths by adding file ending
                my @file_paths = map {
                        $infile_path_prefix{$_}{pdelly_call} . "_"
                      . $sv_type
                      . $suffix{pdelly_call}
                } @{ $active_parameter_href->{sample_ids} };

                Program::Variantcalling::Delly::merge(
                    {
                        infile_paths_ref => \@file_paths,
                        outfile_path     => $outfile_path_prefix . "_"
                          . $sv_type . ".bcf",
                        stdoutfile_path => $xargs_file_path_prefix . "."
                          . $sv_type
                          . ".stdout.txt",
                        stderrfile_path => $xargs_file_path_prefix . "."
                          . $sv_type
                          . ".stderr.txt",
                        sv_type    => $sv_type,
                        min_size   => 0,
                        max_size   => 100000000,
                        FILEHANDLE => $XARGSFILEHANDLE,
                    }
                );
                say {$XARGSFILEHANDLE} "\n";
            }
        }

        ## Delly call regenotype
        say {$FILEHANDLE} "## delly call regenotype";

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Create file commands for xargs
            ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
                {
                    FILEHANDLE         => $FILEHANDLE,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    file_path          => $file_path,
                    program_info_path  => $program_info_path,
                    core_number        => $core_number,
                    xargs_file_counter => $xargs_file_counter,
                }
            );
          SV_TYPE:
            foreach my $sv_type ( @{ $active_parameter_href->{delly_types} } ) {

                if ( $sv_type ne "TRA" ) {

                  CONTIG:
                    foreach my $contig (@contigs) {

                        ## Assemble file path
                        my $alignment_sample_file_path =
                          $infile_path_prefix{$sample_id}
                          {pgatk_baserecalibration} . "_"
                          . $contig
                          . $suffix{pgatk_baserecalibration};

                        Program::Variantcalling::Delly::call(
                            {
                                infile_path => $alignment_sample_file_path,
                                genotypefile_path => $outfile_path_prefix . "_"
                                  . $contig . "_"
                                  . $sv_type
                                  . $suffix{pdelly_call},
                                outfile_path => $file_path_prefix{$sample_id}
                                  . "_"
                                  . $contig . "_"
                                  . $sv_type . "_geno"
                                  . $suffix{pdelly_call},
                                stdoutfile_path => $xargs_file_path_prefix . "."
                                  . $contig . "."
                                  . $sv_type
                                  . ".stdout.txt",
                                stderrfile_path => $xargs_file_path_prefix . "."
                                  . $contig . "."
                                  . $sv_type
                                  . ".stderr.txt",
                                sv_type => $sv_type,
                                exclude_file_path =>
                                  $active_parameter_href->{delly_exclude_file},
                                referencefile_path => $active_parameter_href
                                  ->{human_genome_reference},
                                FILEHANDLE => $XARGSFILEHANDLE,
                            }
                        );
                        say {$XARGSFILEHANDLE} "\n";
                    }
                }
                else {

                    ## Assemble file path
                    my $alignment_sample_file_path =
                        $infile_path_prefix{$sample_id}{pgatk_baserecalibration}
                      . $suffix{pgatk_baserecalibration};

                    Program::Variantcalling::Delly::call(
                        {
                            infile_path       => $alignment_sample_file_path,
                            genotypefile_path => $outfile_path_prefix . "_"
                              . $sv_type
                              . $suffix{pdelly_call},
                            outfile_path => $file_path_prefix{$sample_id} . "_"
                              . $sv_type . "_geno"
                              . $suffix{pdelly_call},
                            stdoutfile_path => $xargs_file_path_prefix . "."
                              . $sv_type
                              . ".stdout.txt",
                            stderrfile_path => $xargs_file_path_prefix . "."
                              . $sv_type
                              . ".stderr.txt",
                            sv_type => $sv_type,
                            exclude_file_path =>
                              $active_parameter_href->{delly_exclude_file},
                            referencefile_path =>
                              $active_parameter_href->{human_genome_reference},
                            FILEHANDLE => $XARGSFILEHANDLE,
                        }
                    );
                    say {$XARGSFILEHANDLE} "\n";
                }
            }
        }

        ### Merge calls
        say {$FILEHANDLE} "## bcftools merge";

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                FILEHANDLE         => $FILEHANDLE,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                file_path          => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
                xargs_file_counter => $xargs_file_counter,
            }
        );
      SV_TYPE:
        foreach my $sv_type ( @{ $active_parameter_href->{delly_types} } ) {

            if ( $sv_type ne "TRA" ) {

              CONTIG:
                foreach my $contig (@contigs) {

                    ## Assemble file paths by adding file ending
                    my @file_paths = map {
                            $file_path_prefix{$_} . "_"
                          . $contig . "_"
                          . $sv_type . "_geno"
                          . $suffix{pdelly_call}
                    } @{ $active_parameter_href->{sample_ids} };

                    bcftools_merge(
                        {
                            infile_paths_ref => \@file_paths,
                            outfile_path     => $outfile_path_prefix . "_"
                              . $contig . "_"
                              . $sv_type
                              . $suffix{pdelly_call},
                            output_type     => "b",
                            stdoutfile_path => $xargs_file_path_prefix . "."
                              . $contig . "."
                              . $sv_type
                              . ".stdout.txt",
                            stderrfile_path => $xargs_file_path_prefix . "."
                              . $contig . "."
                              . $sv_type
                              . ".stderr.txt",
                            FILEHANDLE => $XARGSFILEHANDLE,
                        }
                    );
                    print $XARGSFILEHANDLE q{;} . q{ };

                    bcftools_index(
                        {
                            infile_path => $outfile_path_prefix . "_"
                              . $contig . "_"
                              . $sv_type
                              . $suffix{pdelly_call},
                            stderrfile_path => $xargs_file_path_prefix . q{.}
                              . $contig . q{.}
                              . $sv_type . "_"
                              . q{index.stderr.txt},
                            output_type => q{csi},
                            FILEHANDLE  => $XARGSFILEHANDLE,
                        }
                    );
                    say {$XARGSFILEHANDLE} "\n";
                }
            }
            else {

                ## Assemble file paths by adding file ending
                my @file_paths = map {
                        $file_path_prefix{$_} . "_"
                      . $sv_type . "_geno"
                      . $suffix{pdelly_call}
                } @{ $active_parameter_href->{sample_ids} };

                bcftools_merge(
                    {
                        infile_paths_ref => \@file_paths,
                        outfile_path     => $outfile_path_prefix . "_"
                          . $sv_type
                          . $suffix{pdelly_call},
                        output_type     => "b",
                        stdoutfile_path => $xargs_file_path_prefix . "."
                          . $sv_type
                          . ".stdout.txt",
                        stderrfile_path => $xargs_file_path_prefix . "."
                          . $sv_type
                          . ".stderr.txt",
                        FILEHANDLE => $XARGSFILEHANDLE,
                    }
                );
                print $XARGSFILEHANDLE q{;} . q{ };

                bcftools_index(
                    {
                        infile_path => $outfile_path_prefix . "_"
                          . $sv_type
                          . $suffix{pdelly_call},
                        stderrfile_path => $xargs_file_path_prefix . q{.}
                          . $sv_type . "_"
                          . q{index.stderr.txt},
                        output_type => q{csi},
                        FILEHANDLE  => $XARGSFILEHANDLE,
                    }
                );
                say {$XARGSFILEHANDLE} "\n";
            }
        }

        ### Concatenate SV types
        say {$FILEHANDLE}
          "## bcftools concat - concatenate SV type per contigs";

        ## Assemble file paths by adding file ending
        my @file_paths;

      SV_TYPE:
        foreach my $sv_type ( @{ $active_parameter_href->{delly_types} } ) {

            if ( $sv_type ne "TRA" ) {

                push(
                    @file_paths,
                    map {
                            $outfile_path_prefix . "_"
                          . $_ . "_"
                          . $sv_type
                          . $suffix{pdelly_call}
                    } @contigs
                );
            }
            else {

                push( @file_paths,
                        $outfile_path_prefix . "_"
                      . $sv_type
                      . $suffix{pdelly_call} );
            }

            bcftools_concat(
                {
                    infile_paths_ref => \@file_paths,
                    outfile_path     => $outfile_path_prefix . "_"
                      . $sv_type
                      . "_concat"
                      . $suffix{pdelly_call},
                    output_type     => "b",
                    stderrfile_path => $program_info_path . "_"
                      . $sv_type
                      . "_concat.stderr.txt",
                    allow_overlaps => 1,
                    FILEHANDLE     => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} "\n";

            bcftools_index(
                {
                    infile_path => $outfile_path_prefix . "_"
                      . $sv_type
                      . "_concat"
                      . $suffix{pdelly_call},
                    stderrfile_path => $xargs_file_path_prefix . "."
                      . $sv_type
                      . "_index.stderr.txt",
                    output_type => "csi",
                    FILEHANDLE  => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} "\n";
        }

        ### Filter calls
        say {$FILEHANDLE} "## Delly filter";

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                FILEHANDLE         => $FILEHANDLE,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                file_path          => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
                xargs_file_counter => $xargs_file_counter,
            }
        );
      SV_TYPE:
        foreach my $sv_type ( @{ $active_parameter_href->{delly_types} } ) {

            if ( $sv_type ne "TRA" ) {

                Program::Variantcalling::Delly::filter(
                    {
                        infile_path => $outfile_path_prefix . "_"
                          . $sv_type
                          . "_concat"
                          . $suffix{pdelly_call},
                        outfile_path => $outfile_path_prefix . "_"
                          . $sv_type
                          . "_filtered"
                          . $suffix{pdelly_call},
                        stdoutfile_path => $xargs_file_path_prefix . "."
                          . $sv_type
                          . ".stdout.txt",
                        stderrfile_path => $xargs_file_path_prefix . "."
                          . $sv_type
                          . ".stderr.txt",
                        sv_type     => $sv_type,
                        filter_mode => "germline",
                        FILEHANDLE  => $XARGSFILEHANDLE,
                    }
                );
                say {$XARGSFILEHANDLE} "\n";
            }
            else {

                Program::Variantcalling::Delly::filter(
                    {
                        infile_path => $outfile_path_prefix . "_"
                          . $sv_type
                          . $suffix{pdelly_call},
                        outfile_path => $outfile_path_prefix . "_"
                          . $sv_type
                          . "_filtered"
                          . $suffix{pdelly_call},
                        stdoutfile_path => $xargs_file_path_prefix . "."
                          . $sv_type
                          . ".stdout.txt",
                        stderrfile_path => $xargs_file_path_prefix . "."
                          . $sv_type
                          . ".stderr.txt",
                        sv_type     => $sv_type,
                        filter_mode => "germline",
                        FILEHANDLE  => $XARGSFILEHANDLE,
                    }
                );
                say {$XARGSFILEHANDLE} "\n";
            }
        }
    }

    ## Assemble filepaths
    my @file_paths;

    ### Concatenate SV types
    if ( scalar( @{ $active_parameter_href->{sample_ids} } ) > 1 ) {

        say {$FILEHANDLE} "## bcftools concat - merge all SV types";

        @file_paths = map {
                $outfile_path_prefix . "_"
              . $_
              . "_filtered"
              . $suffix{pdelly_call}
        } @{ $active_parameter_href->{delly_types} };
    }
    else {    #Only one sample

        say {$FILEHANDLE} "## Only one sample - skip merging and regenotyping";
        say {$FILEHANDLE} "## bcftools concat - merge all SV types and contigs";

      SV_TYPE:
        foreach my $sv_type ( @{ $active_parameter_href->{delly_types} } ) {

            if ( $sv_type ne "TRA" ) {

              CONTIG:
                foreach my $contig (@contigs) {

                    ## Assemble file paths by adding file ending
                    push(
                        @file_paths,
                        map {
                                $infile_path_prefix{$_}{pdelly_call} . "_"
                              . $contig . "_"
                              . $sv_type
                              . $suffix{pdelly_call}
                        } @{ $active_parameter_href->{sample_ids} }
                    );
                }
            }
            else {

                push(
                    @file_paths,
                    map {
                            $infile_path_prefix{$_}{pdelly_call} . "_"
                          . $sv_type
                          . $suffix{pdelly_call}
                    } @{ $active_parameter_href->{sample_ids} }
                );
            }
        }
    }
    bcftools_concat(
        {
            infile_paths_ref => \@file_paths,
            outfile_path => $outfile_path_prefix . "_concat" . $outfile_suffix,
            output_type  => "v",
            stderrfile_path => $program_info_path . "_concat.stderr.txt",
            allow_overlaps  => 1,
            FILEHANDLE      => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

    ## Writes sbatch code to supplied filehandle to sort variants in vcf format
    sort_vcf(
        {
            active_parameter_href => $active_parameter_href,
            FILEHANDLE            => $FILEHANDLE,
            sequence_dict_file    => catfile(
                $$reference_dir_ref,
                $file_info_href->{human_genome_reference_name_prefix} . ".dict"
            ),
            infile_paths_ref =>
              [ $outfile_path_prefix . "_concat" . $outfile_suffix ],
            outfile => $outfile_path_prefix . $outfile_suffix,
        }
    );

    ## Copies file from temporary directory.
    say {$FILEHANDLE} "\n## Copy file from temporary directory";
    migrate_file(
        {
            infile_path  => $outfile_path_prefix . $outfile_suffix,
            outfile_path => $outfamily_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, "\n";

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => 'delly',
                outdirectory     => $outfamily_directory,
                outfile          => $outfile_prefix . $outfile_suffix,
            }
        );
    }
    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $$family_id_ref,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
}

sub delly_call {

##delly_call

##Function : Call structural variants using delly
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $sample_id_ref, $program_name, $program_info_path, family_id_ref, $temp_directory_ref, $reference_dir_ref, $outaligner_dir_ref, $xargs_file_counter
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $sample_id_ref              => Sample id
##         : $program_name               => Program name
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $reference_dir_ref          => MIP reference directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $xargs_file_counter         => The xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $reference_dir_ref;
    my $outaligner_dir_ref;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $sample_id_ref;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        sample_id_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$sample_id_ref
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        reference_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{reference_dir},
            strict_type => 1,
            store       => \$reference_dir_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use MIP::IO::Files qw(migrate_file xargs_migrate_contig_files);
    use MIP::Get::File qw{get_file_suffix get_merged_infile_prefix };
    use MIP::Set::File qw{set_file_suffix};
    use MIP::Delete::List qw{ delete_contig_elements };
    use Program::Variantcalling::Delly qw(call);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_sample);

    my $core_number =
      $active_parameter_href->{module_core_number}{ "p" . $program_name };
    my $program_outdirectory_name =
      $parameter_href->{ "p" . $program_name }{outdir_name};
    my $xargs_file_path_prefix;
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $FILEHANDLE      = IO::Handle->new();    #Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$sample_id_ref,
            program_name          => $program_name,
            program_directory     => catfile(
                lc($$outaligner_dir_ref),
                lc($program_outdirectory_name)
            ),
            core_number => $core_number,
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
            temp_directory => $$temp_directory_ref
        }
    );

    ## Assign directories
    my $insample_directory = catdir( $active_parameter_href->{outdata_dir},
        $$sample_id_ref, $$outaligner_dir_ref );
    my $outsample_directory = catdir( $active_parameter_href->{outdata_dir},
        $$sample_id_ref, $$outaligner_dir_ref, $program_outdirectory_name );
    $parameter_href->{ "p" . $program_name }{$$sample_id_ref}{indirectory} =
      $outsample_directory;    #Used downstream

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $$sample_id_ref,
        }
    );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$$sample_id_ref}{pgatk_baserecalibration}{file_tag};
    my $outfile_tag =
      $file_info_href->{$$sample_id_ref}{ "p" . $program_name }{file_tag};
    my $infile_prefix       = $merged_infile_prefix . $infile_tag;
    my $file_path_prefix    = catfile( $$temp_directory_ref, $infile_prefix );
    my $outfile_prefix      = $merged_infile_prefix . $outfile_tag;
    my $outfile_path_prefix = catfile( $$temp_directory_ref, $outfile_prefix );

    ## Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
            jobid_chain    => $parameter_href->{pgatk_baserecalibration}{chain}
            ,    #Get infile_suffix from baserecalibration jobid chain
        }
    );
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    ### Update contigs
    ## Removes an element from array and return new array while leaving orginal elements_ref untouched
    my @contigs = delete_contig_elements(
        {
            elements_ref       => \@{ $file_info_href->{contigs_size_ordered} },
            remove_contigs_ref => [qw{ MT M }],
        }
    );

    ## Removes contigs from supplied contigs_ref
    remove_array_element(
        {
            contigs_ref        => \@contigs,
            remove_contigs_ref => ["Y"]
            , #Skip contig Y throughout since sometimes there are no variants particularly for INS
        }
    );

    if ( ( any { $_ eq "TRA" } @{ $active_parameter_href->{delly_types} } ) )
    {         #If element is part of array

        ## Required for processing complete file (TRA)
        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        migrate_file(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => catfile(
                    $insample_directory,
                    $infile_prefix . substr( $infile_suffix, 0, 2 ) . q{*}
                ),
                outfile_path => $active_parameter_href->{temp_directory}
            }
        );
    }
    if (
        (
            any { $_ =~ /DEL|DUP|INS|INV/ }
            @{ $active_parameter_href->{delly_types} }
        )
      )
    {    #If element is part of array

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        ($xargs_file_counter) = xargs_migrate_contig_files(
            {
                FILEHANDLE        => $FILEHANDLE,
                XARGSFILEHANDLE   => $XARGSFILEHANDLE,
                contigs_ref       => \@contigs,
                file_path         => $file_path,
                program_info_path => $program_info_path,
                core_number       => ( $core_number - 1 )
                ,    #Compensate for cp of entire BAM (INS, TRA), see above
                xargs_file_counter => $xargs_file_counter,
                infile             => $infile_prefix,
                indirectory        => $insample_directory,
                file_ending        => substr( $infile_suffix, 0, 2 ) . "*",
                temp_directory     => $$temp_directory_ref,
            }
        );
        say {$FILEHANDLE} q{wait}, "\n";
    }

    ## delly
    say {$FILEHANDLE} "## delly";

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    foreach my $sv_type ( @{ $active_parameter_href->{delly_types} } ) {

        if ( $sv_type ne "TRA" ) {

            ## Process per contig
            foreach my $contig (@contigs) {

                Program::Variantcalling::Delly::call(
                    {
                        infile_path => $file_path_prefix . "_"
                          . $contig
                          . $infile_suffix,
                        outfile_path => $outfile_path_prefix . "_"
                          . $contig . "_"
                          . $sv_type
                          . $outfile_suffix,
                        stdoutfile_path => $xargs_file_path_prefix . "."
                          . $contig . "."
                          . $sv_type
                          . ".stdout.txt",
                        stderrfile_path => $xargs_file_path_prefix . "."
                          . $contig . "."
                          . $sv_type
                          . ".stderr.txt",
                        sv_type => $sv_type,
                        exclude_file_path =>
                          $active_parameter_href->{delly_exclude_file},
                        referencefile_path =>
                          $active_parameter_href->{human_genome_reference},
                        FILEHANDLE => $XARGSFILEHANDLE,
                    }
                );
                say {$XARGSFILEHANDLE} "\n";
            }
        }
        else {

            Program::Variantcalling::Delly::call(
                {
                    infile_path  => $file_path_prefix . $infile_suffix,
                    outfile_path => $outfile_path_prefix . "_"
                      . $sv_type
                      . $outfile_suffix,
                    stdoutfile_path => $xargs_file_path_prefix . "."
                      . $sv_type
                      . ".stdout.txt",
                    stderrfile_path => $xargs_file_path_prefix . "."
                      . $sv_type
                      . ".stderr.txt",
                    sv_type => $sv_type,
                    exclude_file_path =>
                      $active_parameter_href->{delly_exclude_file},
                    referencefile_path =>
                      $active_parameter_href->{human_genome_reference},
                    FILEHANDLE => $XARGSFILEHANDLE,
                }
            );
            say {$XARGSFILEHANDLE} "\n";
        }
    }

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            infile_path => $outfile_path_prefix . q{*} . $outfile_suffix . q{*},
            outfile_path => $outsample_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, "\n";

    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        slurm_submit_job_sample_id_dependency_add_to_sample(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                family_id               => $$family_id_ref,
                sample_id               => $$sample_id_ref,
                path                    => $job_id_chain,
                log                     => $log,
                sbatch_file_name        => $file_path
            }
        );
    }
}

sub msamtools_mpileup {

## msamtools_mpileup

##Function : samtools_mpileup
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $sample_id_ref, $program_name, family_id_ref, $temp_directory_ref, $outaligner_dir_ref, $call_type, $xargs_file_counter
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $sample_id_ref              => Sample id {REF}
##         : $program_name               => Program name
##         : $family_id_ref              => Family id {REF}
##         : $temp_directory_ref         => Temporary directory {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type
##         : $xargs_file_counter          => The xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
    my $call_type;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $sample_id_ref;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use MIP::IO::Files qw(migrate_file xargs_migrate_contig_files);
    use MIP::Get::File qw{get_file_suffix get_merged_infile_prefix };
    use MIP::Set::File qw{set_file_suffix};
    use MIP::Recipes::Xargs qw{ xargs_command };
    use MIP::Program::Alignment::Samtools qw(samtools_mpileup);
    use MIP::Program::Variantcalling::Bcftools
      qw(bcftools_call bcftools_filter bcftools_norm);
    use MIP::QC::Record qw(add_program_outfile_to_sample_info);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_family);

    my $core_number =
      $active_parameter_href->{module_core_number}{ "p" . $program_name };
    my $program_outdirectory_name =
      $parameter_href->{ "p" . $program_name }{outdir_name};
    my $xargs_file_path_prefix;
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $FILEHANDLE      = IO::Handle->new();    #Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory =>
              catfile( lc($$outaligner_dir_ref), $program_outdirectory_name ),
            core_number => $core_number,
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
            temp_directory => $$temp_directory_ref
        }
    );

    ## Assign directories
    my $outfamily_file_directory =
      catfile( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref, $program_outdirectory_name )
      ;    #For ".fam" file
    my $outfamily_directory = catfile( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref, $program_outdirectory_name );
    $parameter_href->{ "p" . $program_name }{indirectory} =
      $outfamily_directory;    #Used downstream

    ## Assign file_tags
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{ "p" . $program_name }{file_tag};
    my $outfile_prefix = $$family_id_ref . $outfile_tag . $call_type;
    my $outfile_path_prefix = catfile( $$temp_directory_ref, $outfile_prefix );

    ## Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
            jobid_chain    => $parameter_href->{pgatk_baserecalibration}{chain}
            ,    #Get infile_suffix from baserecalibration jobid chain
        }
    );

    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            parameter_href        => $parameter_href,
            active_parameter_href => $active_parameter_href,
            sample_info_href      => $sample_info_href,
            FILEHANDLE            => $FILEHANDLE,
            fam_file_path =>
              catfile( $outfamily_file_directory, $$family_id_ref . ".fam" ),
            include_header => 0,
        }
    );

    my %file_path_prefix;

    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } )
    {    #Collect infiles for all sample_ids

        ## Assign directories
        my $insample_directory = catdir( $active_parameter_href->{outdata_dir},
            $sample_id, $$outaligner_dir_ref );

        ## Add merged infile name prefix after merging all BAM files per sample_id
        my $merged_infile_prefix = get_merged_infile_prefix(
            {
                file_info_href => $file_info_href,
                sample_id      => $sample_id,
            }
        );

        ## Assign file_tags
        my $infile_tag =
          $file_info_href->{$sample_id}{pgatk_baserecalibration}{file_tag};
        my $infile_prefix = $merged_infile_prefix . $infile_tag;

        $file_path_prefix{$sample_id} =
          catfile( $$temp_directory_ref, $infile_prefix );

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        ($xargs_file_counter) = xargs_migrate_contig_files(
            {
                FILEHANDLE      => $FILEHANDLE,
                XARGSFILEHANDLE => $XARGSFILEHANDLE,
                contigs_ref => \@{ $file_info_href->{contigs_size_ordered} },
                file_path   => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
                xargs_file_counter => $xargs_file_counter,
                infile             => $infile_prefix,
                indirectory        => $insample_directory,
                file_ending        => substr( $infile_suffix, 0, 2 )
                  . "*",    #q{.bam} -> ".b*" for getting index as well
                temp_directory => $$temp_directory_ref,
            }
        );
    }

    ## SamTools mpileup
    say {$FILEHANDLE} "## SamTools mpileup";

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Split per contig
    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## Assemble file paths by adding file ending
        my @file_paths =
          map { $file_path_prefix{$_} . "_" . $contig . $infile_suffix }
          @{ $active_parameter_href->{sample_ids} };

        samtools_mpileup(
            {
                infile_paths_ref => \@file_paths,
                output_tags_ref  => [ "DV", "AD" ],
                referencefile_path =>
                  $active_parameter_href->{human_genome_reference},
                stderrfile_path => $xargs_file_path_prefix . "."
                  . $contig
                  . ".stderr.txt",
                output_bcf                       => 1,
                adjust_mq                        => 50,
                per_sample_increased_sensitivity => 1,
                FILEHANDLE                       => $XARGSFILEHANDLE,
            }
        );
        print $XARGSFILEHANDLE "| ";    #Pipe

        ## Get parameter
        my $samples_file;
        my $constrain;
        if ( $parameter_href->{dynamic_parameter}{trio} ) {

            $samples_file =
              catfile( $outfamily_file_directory, $$family_id_ref . ".fam" )
              . " ";
            $constrain = "trio";
        }
        bcftools_call(
            {
                form_fields_ref     => ["GQ"],
                variants_only       => 1,
                multiallelic_caller => 1,
                samples_file        => $samples_file,
                constrain           => $constrain,
                stderrfile_path     => $xargs_file_path_prefix . "."
                  . $contig
                  . "_call.stderr.txt",
                FILEHANDLE => $XARGSFILEHANDLE,
            }
        );
        print $XARGSFILEHANDLE "| ";    #Pipe

        bcftools_filter(
            {
                stderrfile_path => $xargs_file_path_prefix . "."
                  . $contig
                  . "_filter.stderr.txt",
                soft_filter => "LowQual",
                snp_gap     => 3,
                indel_gap   => 10,
                exclude =>
q?\'%QUAL<10 || (RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || %MAX(DV)<=3 || %MAX(DV)/%MAX(DP)<=0.25\'?,
                FILEHANDLE => $XARGSFILEHANDLE,
            }
        );

        if ( $active_parameter_href->{replace_iupac} ) {

            ## Replace the IUPAC code in alternative allels with N for input stream and writes to stream
            replace_iupac(
                {
                    stderr_path => $xargs_file_path_prefix . "."
                      . $contig
                      . "_replace_iupac.stderr.txt",
                    FILEHANDLE => $XARGSFILEHANDLE
                }
            );
        }

        print $XARGSFILEHANDLE "| ";    #Pipe

        ## BcfTools norm, Left-align and normalize indels, split multiallelics
        bcftools_norm(
            {
                FILEHANDLE => $XARGSFILEHANDLE,
                reference_path =>
                  $active_parameter_href->{human_genome_reference},
                output_type  => "v",
                outfile_path => $outfile_path_prefix . "_"
                  . $contig
                  . $outfile_suffix,
                multiallelic    => "-",
                stderrfile_path => catfile(
                        $xargs_file_path_prefix . "."
                      . $contig
                      . "_norm.stderr.txt"
                ),
            }
        );
        say {$XARGSFILEHANDLE} "\n";
    }

    ## Writes sbatch code to supplied filehandle to concatenate variants in vcf format. Each array element is combined with the infile prefix and postfix.
    concatenate_variants(
        {
            active_parameter_href => $active_parameter_href,
            FILEHANDLE            => $FILEHANDLE,
            elements_ref          => \@{ $file_info_href->{contigs} },
            infile_prefix         => $outfile_path_prefix . "_",
            infile_postfix        => $outfile_suffix,
            outfile               => $outfile_path_prefix . $outfile_suffix,
        }
    );

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            infile_path  => $outfile_path_prefix . $outfile_suffix . q{*},
            outfile_path => $outfamily_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, "\n";

    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        ## Collect samtools version in qccollect
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => 'samtools',
                outdirectory     => $outfamily_directory,
                outfile          => $outfile_prefix . $outfile_suffix,
            }
        );
        ## Locating samtools_mpileup file
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => 'samtools_mpileup',
                outdirectory     => $outfamily_directory,
                outfile          => $outfile_prefix . $outfile_suffix,
            }
        );
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => 'bcftools',
                outdirectory     => $outfamily_directory,
                outfile          => $outfile_prefix . $outfile_suffix,
            }
        );

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $$family_id_ref,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
}

sub freebayes {

##freebayes

##Function : Call snv/small indels usig freebayes
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $sample_id_ref, $program_name
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $sample_id_ref              => Sample id {REF}
##         : $program_name               => Program name

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
    my $call_type;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $sample_id_ref;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use MIP::IO::Files qw(migrate_file xargs_migrate_contig_files);
    use MIP::Get::File qw{get_file_suffix get_merged_infile_prefix };
    use MIP::Set::File qw{set_file_suffix};
    use MIP::Recipes::Xargs qw{ xargs_command };
    use Program::Variantcalling::Freebayes qw(calling);
    use MIP::Program::Variantcalling::Bcftools
      qw(bcftools_filter bcftools_norm);
    use MIP::QC::Record qw(add_program_outfile_to_sample_info);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_family);

    my $core_number =
      $active_parameter_href->{module_core_number}{ "p" . $program_name };
    my $program_outdirectory_name =
      $parameter_href->{ "p" . $program_name }{outdir_name};
    my $xargs_file_path_prefix;
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $FILEHANDLE      = IO::Handle->new();    #Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory =>
              catfile( lc($$outaligner_dir_ref), $program_outdirectory_name ),
            core_number => $core_number,
            process_time =>
              $active_parameter_href->{module_time}{ "p" . $program_name },
            temp_directory => $$temp_directory_ref
        }
    );

    ## Assign directories
    my $outfamily_file_directory =
      catfile( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref, $program_outdirectory_name )
      ;    #For ".fam" file
    my $outfamily_directory = catfile( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref, $program_outdirectory_name );
    $parameter_href->{ "p" . $program_name }{indirectory} =
      $outfamily_directory;    #Used downstream

    ## Assign file_tags
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{ "p" . $program_name }{file_tag};
    my $outfile_prefix = $$family_id_ref . $outfile_tag . $call_type;
    my $outfile_path_prefix = catfile( $$temp_directory_ref, $outfile_prefix );

    ## Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
            jobid_chain    => $parameter_href->{pgatk_baserecalibration}{chain}
            ,    #Get infile_suffix from baserecalibration jobid chain
        }
    );

    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    my %file_path_prefix;

  SAMPLE:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Assign directories
        my $insample_directory = catdir( $active_parameter_href->{outdata_dir},
            $sample_id, $$outaligner_dir_ref );

        ## Add merged infile name prefix after merging all BAM files per sample_id
        my $merged_infile_prefix = get_merged_infile_prefix(
            {
                file_info_href => $file_info_href,
                sample_id      => $sample_id,
            }
        );

        ## Assign file_tags
        my $infile_tag =
          $file_info_href->{$sample_id}{pgatk_baserecalibration}{file_tag};
        my $infile_prefix = $merged_infile_prefix . $infile_tag;

        $file_path_prefix{$sample_id} =
          catfile( $$temp_directory_ref, $infile_prefix );

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        ($xargs_file_counter) = xargs_migrate_contig_files(
            {
                FILEHANDLE      => $FILEHANDLE,
                XARGSFILEHANDLE => $XARGSFILEHANDLE,
                contigs_ref => \@{ $file_info_href->{contigs_size_ordered} },
                file_path   => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
                xargs_file_counter => $xargs_file_counter,
                infile             => $infile_prefix,
                indirectory        => $insample_directory,
                file_ending        => substr( $infile_suffix, 0, 2 )
                  . "*",    #q{.bam} -> ".b*" for getting index as well
                temp_directory => $$temp_directory_ref,
            }
        );
    }

    ## Freebayes
    say {$FILEHANDLE} "## Freebayes";

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Split per contig
    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## Assemble file paths by adding file ending
        my @file_paths =
          map { $file_path_prefix{$_} . "_" . $contig . $infile_suffix }
          @{ $active_parameter_href->{sample_ids} };

        calling(
            {
                infile_paths_ref => \@file_paths,
                referencefile_path =>
                  $active_parameter_href->{human_genome_reference},
                stderrfile_path => $xargs_file_path_prefix . "."
                  . $contig
                  . ".stderr.txt",
                apply_standard_filter      => 1,
                calculate_genotype_quality => 1,
                FILEHANDLE                 => $XARGSFILEHANDLE,
            }
        );
        print $XARGSFILEHANDLE "| ";    #Pipe

        bcftools_filter(
            {
                stderrfile_path => $xargs_file_path_prefix . "."
                  . $contig
                  . "_filter.stderr.txt",
                soft_filter => "LowQual",
                snp_gap     => 3,
                indel_gap   => 10,
                exclude     => q?\'%QUAL<10 || (AC<2 && %QUAL<15)\'?,
                FILEHANDLE  => $XARGSFILEHANDLE,
            }
        );

        if ( $active_parameter_href->{replace_iupac} ) {

            ## Replace the IUPAC code in alternative allels with N for input stream and writes to stream
            replace_iupac(
                {
                    stderr_path => $xargs_file_path_prefix . "."
                      . $contig
                      . "_filter.stderr.txt",
                    FILEHANDLE => $XARGSFILEHANDLE
                }
            );
        }
        print $XARGSFILEHANDLE "| ";    #Pipe

        ## BcfTools norm, Left-align and normalize indels, split multiallelics
        bcftools_norm(
            {
                FILEHANDLE => $XARGSFILEHANDLE,
                reference_path =>
                  $active_parameter_href->{human_genome_reference},
                output_type  => "v",
                outfile_path => $outfile_path_prefix . "_"
                  . $contig
                  . $outfile_suffix,
                multiallelic    => "-",
                stderrfile_path => catfile(
                        $xargs_file_path_prefix . "."
                      . $contig
                      . "_norm.stderr.txt"
                ),
            }
        );
        say {$XARGSFILEHANDLE} "\n";
    }

    ## Writes sbatch code to supplied filehandle to concatenate variants in vcf format. Each array element is combined with the infile prefix and postfix.
    concatenate_variants(
        {
            active_parameter_href => $active_parameter_href,
            FILEHANDLE            => $FILEHANDLE,
            elements_ref          => \@{ $file_info_href->{contigs} },
            infile_prefix         => $outfile_path_prefix . "_",
            infile_postfix        => $outfile_suffix,
            outfile               => $outfile_path_prefix . $outfile_suffix,
        }
    );

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            infile_path  => $outfile_path_prefix . $outfile_suffix . q{*},
            outfile_path => $outfamily_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, "\n";

    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => 'freebayes',
                outdirectory     => $outfamily_directory,
                outfile          => $outfile_prefix . $outfile_suffix,
            }
        );
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => 'bcftools',
                outdirectory     => $outfamily_directory,
                outfile          => $outfile_prefix . $outfile_suffix,
            }
        );

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $$family_id_ref,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
}

sub variantannotationblock {

##variantannotationblock

##Function : Run consecutive module
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, family_id_ref, $outaligner_dir_ref, $call_type, $xargs_file_counter
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $program_name               => Program name
##         : $family_id_ref              => Family id {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type
##         : $xargs_file_counter         => The xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $outaligner_dir_ref;
    my $call_type;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);

    my $core_number = $active_parameter_href->{max_cores_per_node};
    my $xargs_file_path_prefix;
    my $time = 80;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Filehandles
    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    if ( $active_parameter_href->{pprepareforvariantannotationblock} > 0 ) {

        $log->info("\t[Prepareforvariantannotationblock]\n");
    }
    if ( $active_parameter_href->{prhocall} > 0 )
    {                                      #Run rhocall. Done per family

        $log->info("\t[rhocall]\n");
    }
    if ( $active_parameter_href->{pvt} > 0 ) {

        $log->info("\t[Vt]\n");            #Run vt. Done per family
    }

    # Run varianteffectpredictor. Family-level
    if ( $active_parameter_href->{pvarianteffectpredictor} > 0 ) {

        $log->info( $TAB . q{[Varianteffectpredictor]} );
    }
    if ( $active_parameter_href->{pvcfparser} > 0 )
    {                                      #Run pvcfparser. Done per family

        $log->info("\t[Vcfparser]\n");
    }
    if ( $active_parameter_href->{psnpeff} > 0 ) {  #Run snpEff. Done per family

        $log->info("\t[Snpeff]\n");
    }
    if ( $active_parameter_href->{prankvariant} > 0 )
    {    #Run rankvariant. Done per family

        $log->info("\t[Rankvariant]\n");
    }
    if ( $active_parameter{pendvariantannotationblock} > 0 )
    {    #Run endvariantannotationblock. Done per family

        $log->info("\t[Endvariantannotationblock]\n");
    }

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory     => lc($$outaligner_dir_ref),
            core_number           => $core_number,
            process_time          => $time,
        }
    );

    ## Copy files for variantannotationblock to enable restart and skip of modules within block
    if ( $active_parameter_href->{pprepareforvariantannotationblock} > 0 ) {

        ($xargs_file_counter) = prepareforvariantannotationblock(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                call_type               => $call_type,
                program_name            => "prepareforvariantannotationblock",
                file_path               => $file_path,
                program_info_path       => $program_info_path,
                FILEHANDLE              => $FILEHANDLE,
                xargs_file_counter      => $xargs_file_counter,
                stderr_path             => $program_info_path . ".stderr.txt",
            }
        );
    }
    if ( $active_parameter_href->{prhocall} > 0 ) {    #Run vt. Done per family

        ($xargs_file_counter) = rhocall(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                call_type               => $call_type,
                program_name            => "rhocall",
                file_path               => $file_path,
                program_info_path       => $program_info_path,
                FILEHANDLE              => $FILEHANDLE,
                xargs_file_counter      => $xargs_file_counter,
                stderr_path             => $program_info_path . ".stderr.txt",
            }
        );
    }
    if ( $active_parameter_href->{pvt} > 0 ) {    #Run vt. Done per family

        ($xargs_file_counter) = vt(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                call_type               => $call_type,
                program_name            => "vt",
                file_path               => $file_path,
                program_info_path       => $program_info_path,
                FILEHANDLE              => $FILEHANDLE,
                xargs_file_counter      => $xargs_file_counter,
                stderr_path             => $program_info_path . ".stderr.txt",
            }
        );
    }

    # Run varianteffectpredictor. Family-level
    if ( $active_parameter_href->{pvarianteffectpredictor} > 0 ) {

        my $program_name = lc q{varianteffectpredictor};

        ($xargs_file_counter) = analysis_vep_rio(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                call_type               => $call_type,
                program_name            => $program_name,
                file_path               => $file_path,
                program_info_path       => $program_info_path,
                FILEHANDLE              => $FILEHANDLE,
                xargs_file_counter      => $xargs_file_counter,
                stderr_path => $program_info_path . $DOT . q{stderr.txt},
            }
        );
    }
    if ( $active_parameter_href->{pvcfparser} > 0 )
    {    #Run vcfparser. Done per family

        ($xargs_file_counter) = mvcfparser(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                call_type               => $call_type,
                program_name            => "vcfparser",
                file_path               => $file_path,
                program_info_path       => $program_info_path,
                FILEHANDLE              => $FILEHANDLE,
                xargs_file_counter      => $xargs_file_counter,
            }
        );
    }
    if ( $active_parameter_href->{psnpeff} > 0 ) {  #Run snpEff. Done per family

        ($xargs_file_counter) = snpeff(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                call_type               => $call_type,
                program_name            => "snpeff",
                file_path               => $file_path,
                program_info_path       => $program_info_path,
                FILEHANDLE              => $FILEHANDLE,
                xargs_file_counter      => $xargs_file_counter,
            }
        );
    }
    if ( $active_parameter_href->{prankvariant} > 0 )
    {    #Run rankvariant. Done per family

        ($xargs_file_counter) = rankvariant(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                call_type               => $call_type,
                program_name            => "rankvariant",
                file_path               => $file_path,
                program_info_path       => $program_info_path,
                FILEHANDLE              => $FILEHANDLE,
                xargs_file_counter      => $xargs_file_counter,
            }
        );
    }
    if ( $active_parameter{pendvariantannotationblock} > 0 )
    {    #Run endvariantannotationblock. Done per family

        ## Run endvariantannotationblock. Done per family
        ($xargs_file_counter) = endvariantannotationblock(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                call_type               => $call_type,
                program_name            => "endvariantannotationblock",
                file_path               => $file_path,
                program_info_path       => $program_info_path,
                FILEHANDLE              => $FILEHANDLE,
                xargs_file_counter      => $xargs_file_counter,
            }
        );
    }
}

sub build_bwa_prerequisites {

##build_bwa_prerequisites

##Function : Creates the BwaPreRequisites using active_parameters{'human_genome_reference'} as reference.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $bwa_build_reference_file_endings_ref, $program_name, $FILEHANDLE, family_id_ref, $temp_directory_ref, $outaligner_dir_ref, $human_genome_reference_ref
##         : $parameter_href                       => Parameter hash {REF}
##         : $active_parameter_href                => Active parameters for this analysis hash {REF}
##         : $sample_info_href                     => Info on samples and family hash {REF}
##         : $file_info_href                       => File info hash {REF}
##         : $infile_lane_prefix_href           => Infile(s) without the ".ending" {REF}
##         : $job_id_href                          => Job id hash {REF}
##         : $bwa_build_reference_file_endings_ref => The bwa reference associated file endings {REF}
##         : $family_id_ref                        => Family ID {REF}
##         : $outaligner_dir_ref                   => The outaligner_dir used in the analysis
##         : $program_name                         => Program name
##         : $family_id_ref                        => Family id {REF}
##         : $temp_directory_ref                   => Temporary directory {REF}
##         : $outaligner_dir_ref                   => Outaligner_dir used in the analysis {REF}
##         : $human_genome_reference_ref           => Human genome reference {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
    my $human_genome_reference_ref;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $bwa_build_reference_file_endings_ref;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        bwa_build_reference_file_endings_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$bwa_build_reference_file_endings_ref
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        human_genome_reference_ref => {
            default =>
              \$arg_href->{active_parameter_href}{human_genome_reference},
            strict_type => 1,
            store       => \$human_genome_reference_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Shell qw{ check_exist_and_move_file };
    use MIP::Recipes::Build::Human_genome_prerequisites
      qw{ build_human_genome_prerequisites };
    use MIP::Script::Setup_script qw(setup_script);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_no_dependency_add_to_samples);

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle
    my $random_integer =
      int( rand(10000) );    #Generate a random integer between 0-10,000.

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$family_id_ref,
            program_name          => $program_name,
            program_directory     => lc($$outaligner_dir_ref),
            process_time          => 3,
        }
    );

    build_human_genome_prerequisites(
        {
            parameter_href          => $parameter_href,
            active_parameter_href   => $active_parameter_href,
            sample_info_href        => $sample_info_href,
            file_info_href          => $file_info_href,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            program                 => $program_name,
            FILEHANDLE              => $FILEHANDLE,
            log                     => $log,
            random_integer          => $random_integer,
        }
    );

    if ( $parameter_href->{bwa_build_reference}{build_file} eq 1 ) {

        $log->warn( "Will try to create required "
              . $$human_genome_reference_ref
              . " index files before executing "
              . $program_name
              . "\n" );

        say {$FILEHANDLE} "## Building BWA index";
        print {$FILEHANDLE} "bwa index ";   #Index sequences in the FASTA format
        print {$FILEHANDLE} "-p "
          . $$human_genome_reference_ref . "_"
          . $random_integer
          . " ";                            #Prefix of the index
        print {$FILEHANDLE} "-a bwtsw ";    #BWT construction algorithm
        say {$FILEHANDLE} $$human_genome_reference_ref,
          "\n";                             #The FASTA reference sequences file

        foreach my $file (@$bwa_build_reference_file_endings_ref) {

            my $intended_file_path = $$human_genome_reference_ref . $file;
            my $temporary_file_path =
              $$human_genome_reference_ref . "_" . $random_integer . $file;

            ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
            check_exist_and_move_file(
                {
                    FILEHANDLE          => $FILEHANDLE,
                    intended_file_path  => $intended_file_path,
                    temporary_file_path => $temporary_file_path,
                }
            );
        }
        $parameter_href->{bwa_build_reference}{build_file} =
          0;    #Ensure that this subrutine is only executed once
    }
    close $FILEHANDLE;

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        my $slurm_path = $parameter_href->{ "p" . $program_name }{chain};

        slurm_submit_job_no_dependency_add_to_samples(
            {
                job_id_href      => $job_id_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $$family_id_ref,
                path             => $slurm_path,
                sbatch_file_name => $file_path,
                log              => $log,
            }
        );
    }
}

sub read_yaml_pedigree_file {

##read_yaml_pedigree_file

##Function : Reads family_id_pedigree file in YAML format. Checks for pedigree data for allowed entries and correct format. Add data to sample_info depending on user info.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $pedigree_href, $file_path, $family_id_ref
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => The associated reference file endings {REF}
##         : $pedigree_href              => Pedigree hash {REF}
##         : $file_path                  => Pedigree file path
##         : $family_id_ref              => Family_id {RF}

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $pedigree_href;
    my $file_path;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        pedigree_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$pedigree_href
        },
        file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_path
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Defines which values are allowed
    my %allowed_values = (
        sex       => [ "male",     "female",     "unknown" ],
        phenotype => [ "affected", "unaffected", "unknown" ],
    );

    my %exom_target_bed_test_file_tracker
      ;    #Use to collect which sample_ids have used a certain capture_kit
    my @pedigree_sample_ids;
    my $family_id = $pedigree_href->{family};
    my @mandatory_family_keys = ( "family", "samples" );
    my @mandatory_sample_keys =
      ( "sample_id", "father", "mother", "sex", "phenotype" );
    my @user_input_sample_ids;

    ### Check input

    my %user_supply_switch = get_user_supplied_info(
        {
            parameter_href        => $parameter_href,
            active_parameter_href => $active_parameter_href,
        }
    );

    if ( $user_supply_switch{sample_ids} != 0 ) {

        ## Prepare CLI supplied sample_ids if comma sep
        my $values_ref = \@{ $parameter_href->{sample_ids}{value} };
        my $element_separator_ref =
          \$parameter_href->{sample_ids}{element_separator};
        @user_input_sample_ids =
          split( $$element_separator_ref,
            join( $$element_separator_ref, @$values_ref ) );
    }

    ## Check that we find mandatory family keys
    foreach my $key (@mandatory_family_keys) {

        if ( !$pedigree_href->{$key} ) {

            $log->fatal( "File: "
                  . $file_path
                  . " cannot find mandatory key: "
                  . $key
                  . " in file\n" );
            exit 1;
        }
    }

    ## Check that supplied cmd and YAML pedigree family_id match
    if ( $pedigree_href->{family} ne $active_parameter_href->{family_id} ) {

        $log->fatal( "File: "
              . $file_path
              . " for  pedigree family_id: '"
              . $pedigree_href->{family}
              . "' and supplied family: '"
              . $active_parameter_href->{family_id}
              . "' does not match\n" );
        exit 1;
    }

    ## Check sample keys and values
    foreach my $pedigree_sample_href ( @{ $pedigree_href->{samples} } ) {

        ## Check that we find mandatory family keys
        foreach my $key (@mandatory_sample_keys) {

            if ( !defined( $pedigree_sample_href->{$key} ) ) {

                $log->fatal( "File: "
                      . $file_path
                      . " cannot find mandatory key: "
                      . $key
                      . " in file\n" );
                exit 1;
            }
            elsif ( $allowed_values{$key} ) {    #Check allowed values

                if (
                    !(
                        any { $_ eq $pedigree_sample_href->{$key} }
                        @{ $allowed_values{$key} }
                    )
                  )
                {    #If element is not part of array

                    $log->fatal(
                        "File: "
                          . $file_path
                          . " found illegal value: "
                          . $pedigree_sample_href->{$key}
                          . " allowed values are '"
                          . join( "' '", @{ $allowed_values{$key} } ),
                        "'\n"
                    );
                    $log->fatal("Please correct the entry before analysis.\n");
                    $log->fatal("\nMIP: Aborting run.\n\n");
                    exit 1;
                }
            }
        }
    }

    ### Add values family level info
    foreach my $key ( keys %$pedigree_href ) {

        unless ( $key eq "samples" ) {

            $sample_info_href->{$key} = $pedigree_href->{$key};
        }
    }

    ### Add values sample level info
    foreach my $pedigree_sample_href ( @{ $pedigree_href->{samples} } ) {

        ## Sample_id
        my $sample_id = $pedigree_sample_href->{sample_id};    #Alias
        push( @pedigree_sample_ids, $sample_id );  #Save pedigree sample_id info

        if ( $user_supply_switch{sample_ids} == 0 ) {

            push( @{ $active_parameter_href->{sample_ids} }, $sample_id )
              ;                                    #Save sample_id info

            ## Reformat pedigree keys to plink format and collect sample info to various hashes
            get_pedigree_sample_info(
                {
                    parameter_href        => $parameter_href,
                    active_parameter_href => $active_parameter_href,
                    sample_info_href      => $sample_info_href,
                    file_info_href        => $file_info_href,
                    exom_target_bed_test_file_tracker_href =>
                      \%exom_target_bed_test_file_tracker,
                    pedigree_sample_href    => $pedigree_sample_href,
                    user_supply_switch_href => \%user_supply_switch,
                    sample_id               => $sample_id,
                }
            );
        }
        else
        { #Save sample_ids in pedigree to check that user supplied info and sample_id in pedigree match

            if ( any { $_ eq $sample_id } @user_input_sample_ids )
            {    #Update sample_id info

                push( @{ $active_parameter_href->{sample_ids} }, $sample_id )
                  ;    #Save sample_id info

                ## Reformat pedigree keys to plink format and collect sample info to various hashes
                get_pedigree_sample_info(
                    {
                        parameter_href        => $parameter_href,
                        active_parameter_href => $active_parameter_href,
                        sample_info_href      => $sample_info_href,
                        file_info_href        => $file_info_href,
                        exom_target_bed_test_file_tracker_href =>
                          \%exom_target_bed_test_file_tracker,
                        pedigree_sample_href    => $pedigree_sample_href,
                        user_supply_switch_href => \%user_supply_switch,
                        sample_id               => $sample_id,
                    }
                );
            }
        }
    }

    ##Check that founder_ids are included in the pedigree info and the analysis run
    check_founder_id(
        {
            pedigree_href => $pedigree_href,
            pedigree_sample_ids_ref =>
              \@{ $active_parameter_href->{sample_ids} },
        }
    );

    if ( !$user_supply_switch{sample_ids} ) {

        @{ $active_parameter_href->{sample_ids} } =
          sort( @{ $active_parameter_href->{sample_ids} } )
          ;    #Lexiographical sort to determine the correct order of ids indata
    }
    else {     #Check that CLI supplied sample_id exists in pedigree

        foreach my $sample_id (@user_input_sample_ids) {

            if ( !( any { $_ eq $sample_id } @pedigree_sample_ids ) )
            {    #If element is not part of array

                $log->fatal(
                    "File: "
                      . $file_path
                      . " provided sample_id: "
                      . $sample_id
                      . " is not present in file",
                    "\n"
                );
                exit 1;
            }
        }

    }
    if (%exom_target_bed_test_file_tracker)
    { #We have read capture kits from pedigree and need to transfer to active_parameters

        foreach
          my $exome_target_bed_file ( keys %exom_target_bed_test_file_tracker )
        {

            $active_parameter_href->{exome_target_bed}{$exome_target_bed_file}
              = join(
                ",",
                @{
                    $exom_target_bed_test_file_tracker{$exome_target_bed_file}
                }
              );
        }
    }
}

sub collect_infiles {

##collect_infiles

##Function : Collects the ".fastq(.gz)" files from the supplied infiles directory. Checks if any files exist.
##Returns  : ""
##Arguments: $active_parameter_href, $indir_path_href, $infile_href
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $indir_path_href       => Indirectories path(s) hash {REF}
##         : $infile_href           => Infiles hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $indir_path_href;
    my $infile_href;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        indir_path_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$indir_path_href
        },
        infile_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_href
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    $log->info("Reads from platform:\n");

    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } )
    {    #Collects inputfiles govern by sample_ids

        ## Return the key if the hash value and query match
        my $infile_directory_ref = \get_matching_values_key(
            {
                active_parameter_href => $active_parameter_href,
                query_value_ref       => \$sample_id,
                parameter_name        => "infile_dirs",
            }
        );

        my @infiles;

        ## Collect all fastq files in supplied indirectories
        my $rule = Path::Iterator::Rule->new;
        $rule->skip_subdirs("original_fastq_files")
          ;    #Ignore if original fastq files sub directory
        $rule->name("*.fastq*");    #Only look for fastq or fastq.gz files
        my $it = $rule->iter($$infile_directory_ref);

        while ( my $file = $it->() ) {    #Iterate over directory

            my ( $volume, $directory, $fastq_file ) =
              File::Spec->splitpath($file);
            push( @infiles, $fastq_file );
        }
        chomp(@infiles);    #Remove newline from every entry in array

        if ( !@infiles ) {  #No "*.fastq*" infiles

            $log->fatal(
"Could not find any '.fastq' files in supplied infiles directory "
                  . $$infile_directory_ref,
                "\n"
            );
            exit 1;
        }
        foreach my $infile (@infiles)
        {    #Check that inFileDirs/infile contains sample_id in filename

            unless ( $infile =~ /$sample_id/ ) {

                $log->fatal(
                    "Could not detect sample_id: "
                      . $sample_id
                      . " in supplied infile: "
                      . $$infile_directory_ref . "/"
                      . $infile,
                    "\n"
                );
                $log->fatal(
"Check that: '--sample_ids' and '--inFileDirs' contain the same sample_id and that the filename of the infile contains the sample_id.",
                    "\n"
                );
                exit 1;
            }
        }
        $log->info( "Sample id: " . $sample_id . "\n" );
        $log->info("\tInputfiles:\n");

        ## Log each file from platform
        foreach my $file (@infiles) {

            $log->info( "\t\t", $file, "\n" );    #Indent for visability
        }
        $indir_path_href->{$sample_id} =
          $$infile_directory_ref;                 #Catch inputdir path
        $infile_href->{$sample_id} = [@infiles];  #Reload files into hash
    }
}

sub infiles_reformat {

##infiles_reformat

##Function : Reformat files for MIP output, which have not yet been created into, correct format so that a sbatch script can be generated with the correct filenames.
##Returns  : "$uncompressed_file_counter"
##Arguments: $active_parameter_href, $sample_info_href, $file_info_href, $infile_href, $indir_path_href, $infile_lane_prefix_href, $infile_both_strands_prefix_href, $lane_href, $job_id_href, $program_name, $outaligner_dir_ref
##         : $active_parameter_href              => Active parameters for this analysis hash {REF}
##         : $sample_info_href                   => Info on samples and family hash {REF}
##         : $file_info_href                     => File info hash {REF}
##         : $infile_href                        => Infiles hash {REF}
##         : $indir_path_href                    => Indirectories path(s) hash {REF}
##         : $infile_lane_prefix_href         => Infile(s) without the ".ending" {REF}
##         : $infile_both_strands_prefix_href => The infile(s) without the ".ending" and strand info {REF}
##         : $lane_href                          => The lane info hash {REF}
##         : $job_id_href                        => Job id hash {REF}
##         : $program_name                       => Program name {REF}
##         : $outaligner_dir_ref                 => Outaligner_dir used in the analysis {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $outaligner_dir_ref;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_href;
    my $indir_path_href;
    my $infile_lane_prefix_href;
    my $infile_both_strands_prefix_href;
    my $lane_href;
    my $job_id_href;
    my $program_name;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_href
        },
        indir_path_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$indir_path_href
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        infile_both_strands_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_both_strands_prefix_href
        },
        lane_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$lane_href
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

# Used to decide later if any inputfiles needs to be compressed before starting analysis
    my $uncompressed_file_counter = 0;

  SAMPLE_ID:
    for my $sample_id ( keys %{$infile_href} ) {

        # Needed to be able to track when lanes are finished
        my $lane_tracker = 0;

        while ( my ( $file_index, $file_name ) =
            each( @{ $infile_href->{$sample_id} } ) )
        {

            ## Check if a file is gzipped.
            my $compressed_switch =
              check_gzipped( { file_name_ref => \$file_name, } );
            my $read_file_command = "zcat";

            if ( !$compressed_switch ) {    #Not compressed

                $uncompressed_file_counter = "uncompressed"
                  ; #File needs compression before starting analysis. Note: All files are rechecked downstream and uncompressed ones are gzipped automatically
                $read_file_command = "cat";
            }

            if (
                $file_name =~ /(\d+)_(\d+)_([^_]+)_([^_]+)_([^_]+)_(\d).fastq/ )
            { #Parse 'new' no "index" format $1=lane, $2=date, $3=Flow-cell, $4=Sample_id, $5=index,$6=direction

                ## Check that the sample_id provided and sample_id in infile name match.
                check_sample_id_match(
                    {
                        active_parameter_href => $active_parameter_href,
                        infile_href           => $infile_href,
                        sample_id             => $sample_id,
                        infile_sample_id => $4,    #$4 = Sample_id from filename
                        file_index => $file_index,
                    }
                );

                ## Adds information derived from infile name to sample_info hash. Tracks the number of lanes sequenced and checks unique array elementents.
                add_infile_info(
                    {
                        active_parameter_href   => $active_parameter_href,
                        sample_info_href        => $sample_info_href,
                        file_info_href          => $file_info_href,
                        lane_href               => $lane_href,
                        infile_href             => $infile_href,
                        indir_path_href         => $indir_path_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        infile_both_strands_prefix_href =>
                          $infile_both_strands_prefix_href,
                        lane              => $1,
                        date              => $2,
                        flowcell          => $3,
                        sample_id         => $4,
                        index             => $5,
                        direction         => $6,
                        lane_tracker_ref  => \$lane_tracker,
                        file_index        => $file_index,
                        compressed_switch => $compressed_switch,
                    }
                );
            }
            else
            {    #No regexp match i.e. file does not follow filename convention

                $log->warn(
                        "Could not detect MIP file name convention for file: "
                      . $file_name
                      . ". \n" );
                $log->warn(
                    "Will try to find mandatory information in fastq header.",
                    "\n" );

                ##Check that file name at least contains sample_id
                if ( $file_name !~ /$sample_id/ ) {

                    $log->fatal(
"Please check that the file name contains the sample_id.",
                        "\n"
                    );
                }

                ## Get run info from fastq file header
                my @fastq_info_headers = get_run_info(
                    {
                        directory         => $indir_path_href->{$sample_id},
                        read_file_command => $read_file_command,
                        file              => $file_name,
                    }
                );

                ## Adds information derived from infile name to sample_info hash. Tracks the number of lanes sequenced and checks unique array elementents.
                add_infile_info(
                    {
                        active_parameter_href   => $active_parameter_href,
                        sample_info_href        => $sample_info_href,
                        file_info_href          => $file_info_href,
                        lane_href               => $lane_href,
                        infile_href             => $infile_href,
                        indir_path_href         => $indir_path_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        infile_both_strands_prefix_href =>
                          $infile_both_strands_prefix_href,
                        lane => $fastq_info_headers[3],
                        date => "000101"
                        , #fastq format does not contain a date of the run, so fake it with constant impossible date
                        flowcell          => $fastq_info_headers[2],
                        sample_id         => $sample_id,
                        index             => $fastq_info_headers[5],
                        direction         => $fastq_info_headers[4],
                        lane_tracker_ref  => \$lane_tracker,
                        file_index        => $file_index,
                        compressed_switch => $compressed_switch,
                    }
                );

                $log->info(
                    "Found following information from fastq header: lane="
                      . $fastq_info_headers[3]
                      . " flow-cell="
                      . $fastq_info_headers[2]
                      . " index="
                      . $fastq_info_headers[5]
                      . " direction="
                      . $fastq_info_headers[4],
                    "\n"
                );
                $log->warn(
"Will add fake date '20010101' to follow file convention since this is not recorded in fastq header\n"
                );
            }
        }
    }
    return $uncompressed_file_counter;
}

sub check_sample_id_match {

##check_sample_id_match

##Function : Check that the sample_id provided and sample_id in infile name match.
##Returns  : ""
##Arguments: $active_parameter_href, $infile_href, $sample_id, $infile_sample_id, $file_index
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $infile_href           => Infiles hash {REF}
##         : $sample_id             => Sample id from user
##         : $infile_sample_id      => Sample_id collect with regexp from infile
##         : $file_index            => Counts the number of infiles

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $infile_href;
    my $sample_id;
    my $infile_sample_id;
    my $file_index;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        infile_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_href
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id,
        },
        infile_sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_sample_id
        },
        file_index => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_index
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my %seen = ( $infile_sample_id => 1 );    #Add input as first increment

    foreach my $sample_id_supplied ( @{ $active_parameter_href->{sample_ids} } )
    {

        $seen{$sample_id_supplied}++;
    }
    unless ( $seen{$infile_sample_id} > 1 ) {

        $log->fatal( $sample_id
              . " supplied and sample_id "
              . $infile_sample_id
              . " found in file : "
              . $infile_href->{$sample_id}[$file_index]
              . " does not match. Please rename file to match sample_id: "
              . $sample_id
              . "\n" );
        exit 1;
    }
}

sub get_run_info {

##get_run_info

##Function : Get run info from fastq file header
##Returns  : ""
##Arguments: $directory, $read_file, $file
##         : $directory       => Directory of file
##         : $read_file_command => Command used to read file
##         : $file            => File to parse

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $directory;
    my $read_file_command;
    my $file;

    my $tmpl = {
        directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$directory
        },
        read_file_command => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$read_file_command
        },
        file =>
          { required => 1, defined => 1, strict_type => 1, store => \$file },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $fastq_header_regexp =
q?perl -nae 'chomp($_); if($_=~/^(@\w+):(\w+):(\w+):(\w+)\S+\s(\w+):\w+:\w+:(\w+)/) {print $1." ".$2." ".$3." ".$4." ".$5." ".$6."\n";} if($.=1) {last;}' ?;

    my $pwd = cwd();      #Save current direcory
    chdir($directory);    #Move to sample_id infile directory

    my $fastq_info_headers = `$read_file_command $file | $fastq_header_regexp;`
      ;                   #Collect fastq header info
    my @fastq_info_headers = split( " ", $fastq_info_headers );

    chdir($pwd);          #Move back to original directory

    unless ( scalar(@fastq_info_headers) eq 6 ) {

        $log->fatal(
"Could not detect reuired sample sequencing run info from fastq file header - PLease proved MIP file in MIP file convention format to proceed\n"
        );
        exit 1;
    }

    return @fastq_info_headers;
}

sub add_infile_info {

##add_infile_info

##Function : Adds information derived from infile name to sample_info hash. Tracks the number of lanes sequenced and checks unique array elementents.
##Returns  : ""
##Arguments: $active_parameter_href, $sample_info_href, $file_info_href, $infile_href, $infile_lane_prefix_href, $infile_both_strands_prefix_href, $indir_path_href, $lane_href, $lane, $date, $flowcell, $sample_id, $index, $direction, $lane_tracker_ref, $file_index, $compressed_switch
##         : $active_parameter_href              => Active parameters for this analysis hash {REF}
##         : $sample_info_href                   => Info on samples and family hash {REF}
##         : $file_info_href                     => File info hash {REF}
##         : $infile_href                        => Infiles hash {REF}
##         : $infile_lane_prefix_href         => Infile(s) without the ".ending" {REF}
##         : $infile_both_strands_prefix_href => The infile(s) without the ".ending" and strand info {REF}
##         : $indir_path_href                    => Indirectories path(s) hash {REF}
##         : $lane_href                          => The lane info hash {REF}
##         : $lane                               => Flow-cell lane
##         : $date                               => Flow-cell sequencing date
##         : $flowcell                           => Flow-cell id
##         : $sample_id                          => Sample id
##         : $index                              => The DNA library preparation molecular barcode
##         : $direction                          => Sequencing read direction
##         : $lane_tracker_ref                   => Counts the number of lanes sequenced {REF}
##         : $file_index                         => Index of file
##         : $compressed_switch                  => ".fastq.gz" or ".fastq" info governs zcat or cat downstream

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_href;
    my $indir_path_href;
    my $infile_lane_prefix_href;
    my $infile_both_strands_prefix_href;
    my $lane_href;
    my $lane_tracker_ref;
    my $sample_id;
    my $lane;
    my $date;
    my $flowcell;
    my $index;
    my $direction;
    my $file_index;
    my $compressed_switch;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_href
        },
        indir_path_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$indir_path_href
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        infile_both_strands_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_both_strands_prefix_href
        },
        lane_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$lane_href
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id,
        },
        lane => {
            required    => 1,
            defined     => 1,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$lane
        },
        lane_tracker_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$lane_tracker_ref
        },
        date =>
          { required => 1, defined => 1, strict_type => 1, store => \$date },
        flowcell => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$flowcell
        },
        index =>
          { required => 1, defined => 1, strict_type => 1, store => \$index },
        direction => {
            required    => 1,
            defined     => 1,
            allow       => [ 1, 2 ],
            strict_type => 1,
            store       => \$direction
        },
        file_index => {
            required    => 1,
            defined     => 1,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$file_index
        },
        compressed_switch => {
            required    => 1,
            defined     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$compressed_switch
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $read_file;
    my $file_at_lane_level_ref;
    my $file_at_direction_level_ref;

    my $parsed_date = Time::Piece->strptime( $date, "%y%m%d" );
    $parsed_date = $parsed_date->ymd;

    if ($compressed_switch) {

        $read_file = "zcat";    #Read file in compressed format
    }
    else {

        $read_file = "cat";     #Read file in uncompressed format
    }

    if ( $direction == 1 ) {    #Read 1

        push( @{ $lane_href->{$sample_id} }, $lane );    #Lane
        $infile_lane_prefix_href->{$sample_id}[$$lane_tracker_ref] =
            $sample_id . "."
          . $date . "_"
          . $flowcell . "_"
          . $index . ".lane"
          . $lane
          ; #Save new format (sample_id_date_flow-cell_index_lane) in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry into hash for both strands and .ending is removed (.fastq).

        $file_at_lane_level_ref =
          \$infile_lane_prefix_href->{$sample_id}[$$lane_tracker_ref];    #Alias
        $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
          {sequence_run_type} = "single_end"; #Single_end until proven otherwise

        ## Collect read length from an infile
        $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
          {sequence_length} = collect_read_length(
            {
                directory         => $indir_path_href->{$sample_id},
                read_file_command => $read_file,
                file              => $infile_href->{$sample_id}[$file_index],
            }
          );

        ## Check if fastq file is interleaved
        $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
          {interleaved} = detect_interleaved(
            {
                directory         => $indir_path_href->{$sample_id},
                read_file_command => $read_file,
                file              => $infile_href->{$sample_id}[$file_index],
            }
          );

        ## Detect "regexp" in string
        $file_info_href->{undetermined_in_file_name}
          { $infile_lane_prefix_href->{$sample_id}[$$lane_tracker_ref] } =
          check_string(
            {
                string => $flowcell,
                regexp => "Undetermined",
            }
          );
        $$lane_tracker_ref++;
    }
    if ( $direction == 2 ) {    #2nd read direction

        $file_at_lane_level_ref =
          \$infile_lane_prefix_href->{$sample_id}[ $$lane_tracker_ref - 1 ]
          ;                     #Alias
        $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
          {sequence_run_type} = 'paired-end'
          ;    #$lane_tracker -1 since it gets incremented after direction eq 1.
    }

    $infile_both_strands_prefix_href->{$sample_id}[$file_index] =
        $sample_id . "."
      . $date . "_"
      . $flowcell . "_"
      . $index . ".lane"
      . $lane . "_"
      . $direction
      ; #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry per strand and .ending is removed (.fastq).

    $file_at_direction_level_ref =
      \$infile_both_strands_prefix_href->{$sample_id}[$file_index];    #Alias
    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{original_file_name}
      = $infile_href->{$sample_id}[$file_index];    #Original file_name

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}
      {original_file_name_prefix} =
        $lane . "_"
      . $date . "_"
      . $flowcell . "_"
      . $sample_id . "_"
      . $index . "_"
      . $direction;    #Original file_name, but no ending

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{lane} =
      $lane;           #Save sample lane

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{date} =
      $parsed_date;    #Save Sequence run date

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{flowcell} =
      $flowcell;       #Save Sequence flow-cell

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{sample_barcode} =
      $index;          #Save sample barcode

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{run_barcode} =
      $date . "_" . $flowcell . "_" . $lane . "_" . $index;    #Save run barcode

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{read_direction} =
      $direction;
}

sub detect_interleaved {

##detect_interleaved

##Function : Detect if fastq file is interleaved
##Returns  : "1(=interleaved)"
##Arguments: $directory, $read_file, $file
##         : $directory         => Directory of file
##         : $read_file_command => Command used to read file
##         : $file              => File to parse

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $directory;
    my $read_file_command;
    my $file;

    my $tmpl = {
        directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$directory
        },
        read_file_command => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$read_file_command
        },
        file =>
          { required => 1, defined => 1, strict_type => 1, store => \$file },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $interleaved_regexp =
q?perl -nae 'chomp($_); if( ($_=~/^@\S+:\w+:\w+:\w+\S+\s(\w+):\w+:\w+:\w+/) && ($.==5) ) {print $1."\n";last;} elsif ($.==6) {last;}' ?;

    my $pwd = cwd();      #Save current direcory
    chdir($directory);    #Move to sample_id infile directory

    my $fastq_info_headers = `$read_file_command $file | $interleaved_regexp;`
      ;                   #Collect interleaved info

    if ( !$fastq_info_headers ) {

        my $interleaved_regexp =
q?perl -nae 'chomp($_); if( ($_=~/^@\w+-\w+:\w+:\w+:\w+:\w+:\w+:\w+\/(\w+)/) && ($.==5) ) {print $1."\n";last;} elsif ($.==6) {last;}' ?;
        $fastq_info_headers = `$read_file_command $file | $interleaved_regexp;`
          ;               #Collect interleaved info
    }

    chdir($pwd);          #Move back to original directory

    unless ( $fastq_info_headers =~ /[1, 2, 3]/ ) {

        $log->fatal("Malformed fastq file!\n");
        $log->fatal( "Read direction is: "
              . $fastq_info_headers
              . " allowed entries are '1', '2', '3'. Please check fastq file\n"
        );
        exit 1;
    }
    if ( $fastq_info_headers > 1 ) {

        $log->info( "Found interleaved fastq file: " . $file, "\n" );
        return 1;
    }
    return;
}

sub add_to_active_parameter {

##add_to_active_parameter

##Function : Checks and sets user input or default values to active_parameters.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $broadcasts_ref, $parameter_name, $associated_programs
##         : $parameter_href        => Holds all parameters
##         : $active_parameter_href => Holds all set parameter for analysis
##         : $sample_info_href      => Info on samples and family hash {REF}
##         : $file_info_href        => File info hash {REF}
##         : $broadcasts_ref        => Holds the parameters info for broadcasting later {REF}
##         : $parameter_name        => Parameter name
##         : $associated_programs   => The parameters program(s) {array, REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $broadcasts_ref;
    my $associated_programs_ref;
    my $parameter_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        broadcasts_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$broadcasts_ref
        },
        associated_programs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$associated_programs_ref
        },
        parameter_name =>
          { required => 1, defined => 1, store => \$parameter_name },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $element_separator_ref =
      \$parameter_href->{$parameter_name}{element_separator};
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};

    foreach my $associated_program (@$associated_programs_ref)
    {    #Check all programs that use parameter

        my $parameter_set_switch = 0;

        if ( defined( $active_parameter_href->{$associated_program} )
            && ( $active_parameter_href->{$associated_program} > 0 ) )
        {    #Only add active programs parameters

            $parameter_set_switch = 1;

            ## Input from cmd
            if (   ( $parameter_href->{$parameter_name}{data_type} eq "ARRAY" )
                && ( defined( $parameter_href->{$parameter_name}{value}[0] ) ) )
            {    #Array reference

                my $values_ref =
                  \@{ $parameter_href->{$parameter_name}{value} };
                @{ $active_parameter_href->{$parameter_name} } =
                  split( $$element_separator_ref,
                    join( $$element_separator_ref, @$values_ref ) );
            }
            elsif (( $parameter_href->{$parameter_name}{data_type} eq "HASH" )
                && ( keys %{ $parameter_href->{$parameter_name}{value} } ) )
            {    #Hash reference

                $active_parameter_href->{$parameter_name} =
                  $parameter_href->{$parameter_name}{value};
            }
            elsif (
                defined( $parameter_href->{$parameter_name}{value} )
                && (
                    ref( $parameter_href->{$parameter_name}{value} ) !~
                    /ARRAY|HASH/ )
              )
            {    #Scalar input from cmd

                $active_parameter_href->{$parameter_name} =
                  $parameter_href->{$parameter_name}{value};
            }
            else {

                if ( defined( $active_parameter_href->{$parameter_name} ) )
                {    #Input from config file

                }
                elsif ( exists( $parameter_href->{$parameter_name}{default} ) )
                {    #Default exists

                    if ( $parameter_href->{$parameter_name}{data_type} eq
                        "ARRAY" )
                    {    #Array reference

                        push(
                            @{ $active_parameter_href->{$parameter_name} },
                            @{ $parameter_href->{$parameter_name}{default} }
                        );
                    }
                    elsif ( $parameter_href->{$parameter_name}{data_type} eq
                        "HASH" )
                    {

                        ## Build default for analysis_type
                        if ( $parameter_name eq "analysis_type" ) {

                            foreach my $sample_id (
                                @{ $active_parameter_href->{sample_ids} } )
                            {

                                $active_parameter_href->{$parameter_name}
                                  {$sample_id} = "wgs";
                            }
                        }
                        ## Build default for infile_dirs
                        elsif ( $parameter_name eq "infile_dirs" ) {

                            foreach my $sample_id (
                                @{ $active_parameter_href->{sample_ids} } )
                            {

                                my $path = catfile(
                                    $active_parameter_href
                                      ->{cluster_constant_path},
                                    $$family_id_ref,
                                    $active_parameter_href->{analysis_type}
                                      {$sample_id},
                                    $sample_id,
                                    "fastq"
                                );
                                $active_parameter_href->{$parameter_name}{$path}
                                  = $sample_id;
                            }
                        }
                        else {

                            $active_parameter_href->{$parameter_name} =
                              $parameter_href->{$parameter_name}{default};
                        }
                    }
                    else {    #Scalar

                        if ( $parameter_name eq "bwa_build_reference" ) {

                            $active_parameter_href->{$parameter_name} =
                              $active_parameter_href->{human_genome_reference}
                              ; #Now we now what human genome refrence to build from
                        }
                        else {

                            $active_parameter_href->{$parameter_name} =
                              $parameter_href->{$parameter_name}{default};
                        }
                    }
                }
                else {          ## No default

                    if (
                        (
                            exists(
                                $parameter_href->{$parameter_name}{mandatory}
                            )
                        )
                        && ( $parameter_href->{$parameter_name}{mandatory} eq
                            "no" )
                      )
                    {           #Not mandatory
                    }
                    else {

                        ## Special cases where the requirement is depending on other variabels
                        if (   ( $parameter_name eq "bwa_mem_rapid_db" )
                            && ( $consensus_analysis_type ne "rapid" ) )
                        { #Do nothing since file is not required unless rapid mode is enabled
                        }
                        elsif (
                            (
                                $parameter_name eq "gatk_genotypegvcfs_ref_gvcf"
                            )
                            && ( $consensus_analysis_type =~ /wgs/ )
                          )
                        { #Do nothing since file is not required unless exome or rapid mode is enabled
                        }
                        elsif (
                            (
                                $parameter_name eq
                                "vcfparser_range_feature_annotation_columns"
                            )
                            && ( $active_parameter_href
                                ->{vcfparser_range_feature_file} eq
                                "nouser_info" )
                          )
                        {    #Do nothing since no SelectFile was given
                        }
                        elsif (
                            (
                                $parameter_name eq
                                "vcfparser_select_feature_annotation_columns"
                            )
                            && ( $active_parameter_href->{vcfparser_select_file}
                                eq "nouser_info" )
                          )
                        {    #Do nothing since no SelectFile was given
                        }
                        elsif (
                            (
                                $parameter_name eq
                                "vcfparser_select_file_matching_column"
                            )
                            && ( $active_parameter_href->{vcfparser_select_file}
                                eq "nouser_info" )
                          )
                        {    #Do nothing since no SelectFile was given
                        }
                        elsif (
                            ( $parameter_name eq "rank_model_file" )
                            && (
                                !defined(
                                    $active_parameter_href->{rank_model_file}
                                )
                            )
                          )
                        { #Do nothing since no rank model was given i.e. use rank scripts deafult supplied with distribution
                        }
                        elsif (
                            (
                                $parameter_name eq
                                "genmod_models_reduced_penetrance_file"
                            )
                            && (
                                !defined(
                                    $active_parameter_href
                                      ->{genmod_models_reduced_penetrance_file}
                                )
                            )
                          )
                        { #Do nothing since no reduced penetrance should be performed
                        }
                        elsif ( $parameter_name eq "exome_target_bed" ) {

                            ## Return a default capture kit as user supplied no info
                            my $capture_kit = add_capture_kit(
                                {
                                    file_info_href => $file_info_href,
                                    supported_capture_kit_href =>
                                      $active_parameter_href
                                      ->{supported_capture_kit},
                                    capture_kit => "latest",
                                }
                            );
                            $active_parameter_href->{$parameter_name}
                              {$capture_kit} = join( ",",
                                @{ $active_parameter_href->{sample_ids} } );

                            ## Update exome_target_bed files with human_genome_reference_source_ref and human_genome_reference_version_ref
                            update_exome_target_bed(
                                {
                                    exome_target_bed_file_href =>
                                      $active_parameter_href
                                      ->{exome_target_bed},
                                    human_genome_reference_source_ref =>
                                      \$file_info_href
                                      ->{human_genome_reference_source},
                                    human_genome_reference_version_ref =>
                                      \$file_info_href
                                      ->{human_genome_reference_version},
                                }
                            );
                            $log->warn(
"Could not detect a supplied capture kit. Will Try to use 'latest' capture kit: "
                                  . $capture_kit,
                                "\n"
                            );
                        }
                        else {

                            if ( defined($log) )
                            {    #We have a logg object and somewhere to write

                                $log->fatal( $USAGE, "\n" );
                                $log->fatal(
                                    "Supply '-"
                                      . $parameter_name
                                      . "' if you want to run "
                                      . $associated_program,
                                    "\n"
                                );
                            }
                            else {

                                warn( $USAGE, "\n" );
                                warn(
                                    "Supply '-"
                                      . $parameter_name
                                      . "' if you want to run "
                                      . $associated_program,
                                    "\n"
                                );
                            }
                            exit 1;
                        }
                    }
                }
            }
        }
        if ( $parameter_set_switch eq 1 )
        {    #No need to set parameter more than once
            last;
        }
    }

    ## Parse Human Genome Reference
    if ( $parameter_name eq "human_genome_reference" ) {

        ## Update path for supplied reference(s) associated with parameter that should reside in the mip reference directory to full path
        update_mip_reference_path(
            {
                parameter_href        => $parameter_href,
                active_parameter_href => $active_parameter_href,
                parameter_name        => "human_genome_reference"
                ,    #Special case to make sure that complete path is supplied
            }
        );

        ## Detect version and source of the human_genome_reference: Source (hg19 or GRCh).
        parse_human_genome_reference(
            {
                file_info_href => $file_info_href,
                human_genome_reference_ref =>
                  \basename( $active_parameter_href->{human_genome_reference} ),
            }
        );

        ## Update exome_target_bed files with human_genome_reference_source_ref and human_genome_reference_version_ref
        update_exome_target_bed(
            {
                exome_target_bed_file_href =>
                  $active_parameter_href->{exome_target_bed},
                human_genome_reference_source_ref =>
                  \$file_info_href->{human_genome_reference_source},
                human_genome_reference_version_ref =>
                  \$file_info_href->{human_genome_reference_version},
            }
        );
    }

    ## Parse pedigree file
    if ( $parameter_name eq "pedigree_file" ) {

        ## Reads family_id_pedigree file in PLINK|YAML format. Checks for pedigree data for allowed entries and correct format. Add data to sample_info depending on user info.
        if ( $active_parameter_href->{pedigree_file} =~ /\.yaml$/ )
        {    #Meta data in YAML format

            ##Loads a YAML file into an arbitrary hash and returns it. Load parameters from previous run from sample_info_file
            my %pedigree = load_yaml(
                { yaml_file => $active_parameter_href->{pedigree_file}, } );
            $log->info( "Loaded: " . $active_parameter_href->{pedigree_file},
                "\n" );

            read_yaml_pedigree_file(
                {
                    parameter_href        => $parameter_href,
                    active_parameter_href => $active_parameter_href,
                    sample_info_href      => $sample_info_href,
                    file_info_href        => $file_info_href,
                    file_path     => $active_parameter_href->{pedigree_file},
                    pedigree_href => \%pedigree,
                }
            );
        }
    }

    ## Parameter set
    if ( defined( $active_parameter_href->{$parameter_name} ) ) {

        my $info = "";    #Hold parameters info

        if ( ref( $active_parameter_href->{$parameter_name} ) eq "ARRAY" )
        {                 #Array reference

            $info = "Set "
              . $parameter_name . " to: "
              . join( $$element_separator_ref,
                @{ $active_parameter_href->{$parameter_name} } );
            push( @$broadcasts_ref, $info );    #Add info to broadcasts
        }
        elsif ( ref( $active_parameter_href->{$parameter_name} ) eq "HASH" ) {

            $info = "Set "
              . $parameter_name . " to: "
              . join( ",",
                map { "$_=$active_parameter_href->{$parameter_name}{$_}" }
                  ( keys %{ $active_parameter_href->{$parameter_name} } ) );
            push( @$broadcasts_ref, $info );    #Add info to broadcasts
        }
        else {

            $info = "Set "
              . $parameter_name . " to: "
              . $active_parameter_href->{$parameter_name};
            push( @$broadcasts_ref, $info );    #Add info to broadcasts
        }
    }
}

sub check_parameter_files {

##check_parameter_files

##Function : Checks that files/directories exists and if file_endings need to be built also updates SampleInfoHash information with information from pedigree
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $broadcasts_ref, $associated_programs_ref, $family_id_ref, $parameter_name $parameter_exists_check
##         : $parameter_href                    => Holds all parameters
##         : $active_parameter_href             => Holds all set parameter for analysis
##         : $sample_info_href                  => Info on samples and family hash {REF}
##         : $file_info_href                    => File info hash {REF}
##         : $broadcasts_ref                    => Holds the parameters info for broadcasting later {REF}
##         : $associated_programs_ref           => The parameters program(s) {REF}
##         : $family_id_ref                     => The family_id_ref {REF}
##         : $parameter_name                    => Parameter name
##         : $parameter_exists_check            => Check if intendent file exists in reference directory

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $broadcasts_ref;
    my $associated_programs_ref;
    my $parameter_name;
    my $parameter_exists_check;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        broadcasts_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$broadcasts_ref
        },
        associated_programs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$associated_programs_ref
        },
        parameter_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_name,
        },
        parameter_exists_check => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_exists_check
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};

    foreach my $associated_program (@$associated_programs_ref)
    {    #Check all programs that use parameter

        my $parameter_set_switch = 0;

        if ( defined( $active_parameter_href->{$associated_program} )
            && ( $active_parameter_href->{$associated_program} > 0 ) )
        {    #Only add active programs parameters

            $parameter_set_switch = 1;

            if ( exists( $parameter_href->{$parameter_name}{reference} ) )
            {    #Expect file to be in reference directory

                ## Update path for supplied reference(s) associated with parameter that should reside in the mip reference directory to full path
                update_mip_reference_path(
                    {
                        parameter_href        => $parameter_href,
                        active_parameter_href => $active_parameter_href,
                        parameter_name        => $parameter_name,
                    }
                );
            }

            if ( defined( $active_parameter_href->{$parameter_name} ) )
            {    #Check parameter existence

                if ( $parameter_href->{$parameter_name}{data_type} eq "SCALAR" )
                {

                    my $path =
                      catfile( $active_parameter_href->{$parameter_name} );

                    if ( $parameter_name eq "bwa_build_reference" ) {

                        ## Checks files to be built by combining filename stub with fileendings
                        check_file_endings_to_build(
                            {
                                parameter_href        => $parameter_href,
                                active_parameter_href => $active_parameter_href,
                                file_endings_ref      => \@{
                                    $file_info_href
                                      ->{bwa_build_reference_file_endings}
                                },
                                parameter_name => "bwa_build_reference",
                                file_name      => $active_parameter_href
                                  ->{human_genome_reference},
                            }
                        );
                    }
                    elsif ( $parameter_name eq "sample_info_file" ) {

                        if (
                            defined(
                                $active_parameter_href->{sample_info_file}
                            )
                          )
                        {

                            if ( -f $active_parameter_href->{sample_info_file} )
                            {

                                if (
                                    defined(
                                        $active_parameter_href->{log_file}
                                    )
                                  )
                                {

                                    $log->info(
                                        "Loaded: "
                                          . $active_parameter_href
                                          ->{sample_info_file},
                                        "\n"
                                    );
                                }

                                ##Loads a YAML file into an arbitrary hash and returns it. Load parameters from previous run from sample_info_file
                                my %temp = load_yaml(
                                    {
                                        yaml_file => $active_parameter_href
                                          ->{sample_info_file},
                                    }
                                );

                                ## Update sample_info with information from pedigree
                                update_sample_info_hash(
                                    {
                                        sample_info_href => $sample_info_href,
                                        temp_href        => \%temp,
                                        family_id_ref    => $family_id_ref,
                                    }
                                );
                            }
                        }
                    }
                    elsif (
                        ( $parameter_name eq "genomic_set" )
                        && ( $active_parameter_href->{genomic_set} eq
                            "nouser_info" )
                      )
                    {    #Do nothing since this is not a required feature
                    }
                    elsif (( $parameter_name eq "bwa_mem_rapid_db" )
                        && ( $consensus_analysis_type ne "rapid" ) )
                    { #Do nothing since file is not required unless rapid mode is enabled
                    }
                    elsif (( $parameter_name eq "gatk_genotypegvcfs_ref_gvcf" )
                        && ( $consensus_analysis_type =~ /wgs/ ) )
                    { #Do nothing since file is not required unless exome mode is enabled
                    }
                    elsif (
                        ( $parameter_name eq "vcfparser_range_feature_file" )
                        && ( $active_parameter_href
                            ->{vcfparser_range_feature_file} eq "nouser_info" )
                      )
                    {    #Do nothing since no RangeFile was given
                    }
                    elsif ($parameter_name eq "vcfparser_select_file"
                        || $parameter_name eq "sv_vcfparser_select_file" )
                    {

                        if ( $active_parameter_href->{vcfparser_select_file} eq
                            "nouser_info" )
                        {    #No select_file was given

                            $active_parameter_href->{vcfparser_outfile_count} =
                              1
                              ; #To track if vcfparser was used with a vcfparser_select_file (=2) or not (=1)
                        }
                        else { #To enable addition of select_file to sample_info

                            check_existance(
                                {
                                    parameter_href => $parameter_href,
                                    active_parameter_href =>
                                      $active_parameter_href,
                                    item_name_ref => \$active_parameter_href
                                      ->{$parameter_name},
                                    parameter_name_ref => \$parameter_name,
                                    item_type_to_check =>
                                      $parameter_exists_check,
                                }
                            );

                            ## Collects sequences contigs used in select file
                            collect_select_file_contigs(
                                {
                                    contigs_ref => \@{
                                        $file_info_href->{select_file_contigs}
                                    },
                                    select_file_path => catfile(
                                        $active_parameter_href
                                          ->{vcfparser_select_file}
                                    ),
                                }
                            );

                            $active_parameter_href->{vcfparser_outfile_count} =
                              2
                              ; #To track if vcfparser was used with a vcfparser_select_file (=2) or not (=1)
                        }
                    }
                    elsif (
                        (
                            $parameter_name eq
                            "genmod_models_reduced_penetrance_file"
                        )
                        && (
                            !defined(
                                $active_parameter_href
                                  ->{genmod_models_reduced_penetrance_file}
                            )
                        )
                      )
                    { #Do nothing since no reduced penetrance should be performed
                    }
                    elsif ( $parameter_name eq "rank_model_file" ) {

                        if (
                            !defined(
                                $active_parameter_href->{rank_model_file}
                            )
                          )
                        { #Do nothing since no rank model config file was given. Use default supplied by ranking script
                        }
                        else
                        { #To enable addition of rankModel file and version to sample_info

                            check_existance(
                                {
                                    parameter_href => $parameter_href,
                                    active_parameter_href =>
                                      $active_parameter_href,
                                    item_name_ref      => \$path,
                                    parameter_name_ref => \$parameter_name,
                                    item_type_to_check =>
                                      $parameter_exists_check,
                                }
                            );
                        }
                    }
                    else {

                        check_existance(
                            {
                                parameter_href        => $parameter_href,
                                active_parameter_href => $active_parameter_href,
                                item_name_ref         => \$path,
                                parameter_name_ref    => \$parameter_name,
                                item_type_to_check => $parameter_exists_check,
                            }
                        );

                        if ( $path =~ /\.gz$/ ) { #Check for tabix index as well

                            $path .= ".tbi";
                            check_existance(
                                {
                                    parameter_href => $parameter_href,
                                    active_parameter_href =>
                                      $active_parameter_href,
                                    item_name_ref      => \$path,
                                    parameter_name_ref => \$parameter_name,
                                    item_type_to_check =>
                                      $parameter_exists_check,
                                }
                            );
                        }

                        if ( $parameter_name eq "human_genome_reference" ) {

                            ## Check the existance of associated Human genome files
                            check_human_genome_file_endings(
                                {
                                    parameter_href => $parameter_href,
                                    active_parameter_href =>
                                      $active_parameter_href,
                                    file_info_href     => $file_info_href,
                                    parameter_name_ref => \$parameter_name
                                }
                            );
                        }
                    }
                }
                if ( $parameter_href->{$parameter_name}{data_type} eq "ARRAY" )
                {

                    foreach my $path (
                        @{ $active_parameter_href->{$parameter_name} } )
                    {

                        check_existance(
                            {
                                parameter_href        => $parameter_href,
                                active_parameter_href => $active_parameter_href,
                                item_name_ref         => \$path,
                                parameter_name_ref    => \$parameter_name,
                                item_type_to_check => $parameter_exists_check,
                            }
                        );

                        if ( $path =~ /\.gz$/ ) { #Check for tabix index as well

                            my $path_index = $path . ".tbi";

                            check_existance(
                                {
                                    parameter_href => $parameter_href,
                                    active_parameter_href =>
                                      $active_parameter_href,
                                    item_name_ref      => \$path_index,
                                    parameter_name_ref => \$path_index,
                                    item_type_to_check =>
                                      $parameter_exists_check,
                                }
                            );
                        }
                    }
                }
                if ( $parameter_href->{$parameter_name}{data_type} eq "HASH" ) {

                    for my $path (
                        keys %{ $active_parameter_href->{$parameter_name} } )
                    {

                        if ( $parameter_name eq "exome_target_bed" ) {

                            ## Check that supplied target file ends with ".bed" and otherwise exists
                            check_target_bed_file_exist(
                                {
                                    active_parameter_href =>
                                      $active_parameter_href,
                                    file           => $path,
                                    parameter_name => $parameter_name,
                                }
                            );

                            ## Checks files to be built by combining filename stub with fileendings
                            check_file_endings_to_build(
                                {
                                    parameter_href => $parameter_href,
                                    active_parameter_href =>
                                      $active_parameter_href,
                                    file_endings_ref =>
                                      \@{ $file_info_href->{exome_target_bed} },
                                    parameter_name => "exome_target_bed",
                                    file_name      => $path,
                                }
                            );
                        }
                        check_existance(
                            {
                                parameter_href        => $parameter_href,
                                active_parameter_href => $active_parameter_href,
                                item_name_ref         => \$path,
                                parameter_name_ref    => \$path,
                                item_type_to_check => $parameter_exists_check,
                            }
                        );

                        if ( $path =~ /\.gz$/ ) { #Check for tabix index as well

                            my $path_index = $path . ".tbi";
                            check_existance(
                                {
                                    parameter_href => $parameter_href,
                                    active_parameter_href =>
                                      $active_parameter_href,
                                    item_name_ref      => \$path,
                                    parameter_name_ref => \$path_index,
                                    item_type_to_check =>
                                      $parameter_exists_check,
                                }
                            );
                        }
                    }
                }
            }
        }
        if ( $parameter_set_switch eq 1 )
        {    #No need to set parameter more than once
            last;
        }
    }
}

sub create_file_endings {

##create_file_endings

##Function : Creates the file_tags depending on which modules are used by the user to relevant chain.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $file_info_href, $infile_lane_prefix_href, $order_parameters_ref, $family_id_ref
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $file_info_href             => Info on files hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $order_parameters_ref       => Order of addition to parameter array {REF}
##         : $family_id_ref              => Family id {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $order_parameters_ref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        order_parameters_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$order_parameters_ref
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my %temp_file_ending
      ;    #Used to enable seqential build-up of file_tags between modules

    foreach my $order_parameter_element (@$order_parameters_ref) {

        if ( defined( $active_parameter_href->{$order_parameter_element} ) )
        {    #Only active parameters

            if (
                (
                    any { $_ eq $order_parameter_element }
                    @{ $parameter_href->{dynamic_parameter}{program} }
                )
              )
            {    #Only process programs

                if ( $parameter_href->{$order_parameter_element}{chain} eq
                    "MAIN" )
                {    #MAIN chain

                    if ( $parameter_href->{$order_parameter_element}{file_tag}
                        ne "nofile_tag" )
                    {    #File_tag exist

                        my $file_ending_ref =
                          \$parameter_href->{$order_parameter_element}
                          {file_tag};    #Alias

###MAIN/Per sample_id
                        foreach my $sample_id (
                            @{ $active_parameter_href->{sample_ids} } )
                        {

                            if ( $active_parameter_href
                                ->{$order_parameter_element} > 0 )
                            {            #File_ending should be added

                                if ( $order_parameter_element eq
                                    q{ppicardtools_mergesamfiles} )
                                {        #Special case

                                    $file_info_href->{$sample_id}
                                      {ppicardtools_mergesamfiles}{file_tag} =
                                      $temp_file_ending{$sample_id} . "";
                                }
                                else {

                                    if (
                                        defined(
                                            $temp_file_ending{$sample_id}
                                        )
                                      )
                                    {

                                        $file_info_href->{$sample_id}
                                          {$order_parameter_element}{file_tag}
                                          = $temp_file_ending{$sample_id}
                                          . $$file_ending_ref;
                                    }
                                    else
                                    {    #First module that should add filending

                                        $file_info_href->{$sample_id}
                                          {$order_parameter_element}{file_tag}
                                          = $$file_ending_ref;
                                    }
                                }
                            }
                            else {       #Do not add new module file_tag

                                $file_info_href->{$sample_id}
                                  {$order_parameter_element}{file_tag} =
                                  $temp_file_ending{$sample_id};
                            }
                            $temp_file_ending{$sample_id} =
                              $file_info_href->{$sample_id}
                              {$order_parameter_element}{file_tag}
                              ;    #To enable sequential build-up of fileending
                        }

###MAIN/Per family_id
                        if ( $active_parameter_href->{$order_parameter_element}
                            > 0 )
                        {          #File_ending should be added

                            if ( $order_parameter_element eq
                                "ppicardtools_mergesamfiles" )
                            {      #Special case - do nothing
                            }
                            else {

                                if (
                                    defined(
                                        $temp_file_ending{$$family_id_ref}
                                    )
                                  )
                                {

                                    $file_info_href->{$$family_id_ref}
                                      {$order_parameter_element}{file_tag} =
                                        $temp_file_ending{$$family_id_ref}
                                      . $$file_ending_ref;
                                }
                                else {   #First module that should add filending

                                    $file_info_href->{$$family_id_ref}
                                      {$order_parameter_element}{file_tag} =
                                      $$file_ending_ref;
                                }
                                $temp_file_ending{$$family_id_ref} =
                                  $file_info_href->{$$family_id_ref}
                                  {$order_parameter_element}{file_tag}
                                  ; #To enable sequential build-up of fileending
                            }
                        }
                        else {      #Do not add new module file_tag

                            $file_info_href->{$$family_id_ref}
                              {$order_parameter_element}{file_tag} =
                              $temp_file_ending{$$family_id_ref};
                        }
                    }
                }
                if ( $parameter_href->{$order_parameter_element}{chain} ne
                    "MAIN" )
                {                   #Other chain(s)

                    my $chain_fork =
                      $parameter_href->{$order_parameter_element}{chain};

                    if ( $parameter_href->{$order_parameter_element}{file_tag}
                        ne "nofile_tag" )
                    {               #File_tag exist

                        my $file_ending_ref =
                          \$parameter_href->{$order_parameter_element}
                          {file_tag};    #Alias

###OTHER/Per sample_id
                        foreach my $sample_id (
                            @{ $active_parameter_href->{sample_ids} } )
                        {

                            if ( $active_parameter_href
                                ->{$order_parameter_element} > 0 )
                            {            #File_ending should be added

                                unless (
                                    defined(
                                        $temp_file_ending{$chain_fork}
                                          {$sample_id}
                                    )
                                  )
                                {

                                    $temp_file_ending{$chain_fork}{$sample_id}
                                      = $temp_file_ending{$sample_id}
                                      ;    #Inherit current MAIN chain.
                                }
                                if (
                                    defined(
                                        $temp_file_ending{$chain_fork}
                                          {$sample_id}
                                    )
                                  )
                                {

                                    $file_info_href->{$sample_id}
                                      {$order_parameter_element}{file_tag} =
                                      $temp_file_ending{$chain_fork}{$sample_id}
                                      . $$file_ending_ref;
                                }
                                else {   #First module that should add filending

                                    $file_info_href->{$sample_id}
                                      {$order_parameter_element}{file_tag} =
                                      $$file_ending_ref;
                                }
                            }
                            else {       #Do not add new module file_tag

                                $file_info_href->{$sample_id}
                                  {$order_parameter_element}{file_tag} =
                                  $temp_file_ending{$chain_fork}{$sample_id};
                            }
                            $temp_file_ending{$chain_fork}{$sample_id} =
                              $file_info_href->{$sample_id}
                              {$order_parameter_element}{file_tag}
                              ;    #To enable sequential build-up of fileending
                        }
###Other/Per family_id

                        if ( $active_parameter_href->{$order_parameter_element}
                            > 0 )
                        {          #File ending should be added

                            unless (
                                defined(
                                    $temp_file_ending{$chain_fork}
                                      {$$family_id_ref}
                                )
                              )
                            {

                                $temp_file_ending{$chain_fork}{$$family_id_ref}
                                  = $temp_file_ending{$$family_id_ref}
                                  ;    #Inherit current MAIN chain.
                            }
                            if (
                                defined(
                                    $temp_file_ending{$chain_fork}
                                      {$$family_id_ref}
                                )
                              )
                            {

                                $file_info_href->{$$family_id_ref}
                                  {$order_parameter_element}{file_tag} =
                                  $temp_file_ending{$chain_fork}
                                  {$$family_id_ref} . $$file_ending_ref;
                            }
                            else {    #First module that should add filending

                                $file_info_href->{$$family_id_ref}
                                  {$order_parameter_element}{file_tag} =
                                  $$file_ending_ref;
                            }
                            $temp_file_ending{$chain_fork}{$$family_id_ref} =
                              $file_info_href->{$$family_id_ref}
                              {$order_parameter_element}{file_tag}
                              ;    #To enable sequential build-up of fileending
                        }
                        else {     #Do not add new module file_tag

                            $file_info_href->{$$family_id_ref}
                              {$order_parameter_element}{file_tag} =
                              $temp_file_ending{$chain_fork}{$$family_id_ref};
                        }
                    }
                }
            }
        }
    }
}

sub gatk_pedigree_flag {

##gatk_pedigree_flag

##Function : Check if "--pedigree" and "--pedigreeValidationType" should be included in analysis
##Returns  : "%command"
##Arguments: $active_parameter_href, $fam_file_path, $program_name, $pedigree_validation_type
##         : $active_parameter_href    => Active parameters for this analysis hash {REF}
##         : $fam_file_path            => The family file path
##         : $program_name             => The program to use the pedigree file
##         : $pedigree_validation_type => The pedigree validation strictness level

    my ($arg_href) = @_;

    ## Default(s)
    my $pedigree_validation_type;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $fam_file_path;
    my $program_name;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        fam_file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$fam_file_path
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        pedigree_validation_type => {
            default     => "SILENT",
            allow       => [ "SILENT", "STRICT" ],
            strict_type => 1,
            store       => \$pedigree_validation_type
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $parent_counter;
    my $pq_parent_counter =
q?perl -ne 'my $parent_counter=0; while (<>) { my @line = split(/\t/, $_); unless ($_=~/^#/) { if ( ($line[2] eq 0) || ($line[3] eq 0) ) { $parent_counter++} } } print $parent_counter; last;'?;
    my $child_counter;
    my $pq_child_counter =
q?perl -ne 'my $child_counter=0; while (<>) { my @line = split(/\t/, $_); unless ($_=~/^#/) { if ( ($line[2] ne 0) || ($line[3] ne 0) ) { $child_counter++} } } print $child_counter; last;'?;
    my %command;

    $parent_counter =
      `$pq_parent_counter $fam_file_path`;    #Count the number of parents
    $child_counter =
      `$pq_child_counter $fam_file_path`;     #Count the number of children

    if ( $parent_counter > 0 ) {              #Parents present

        $command{pedigree_validation_type} = $pedigree_validation_type;
        $command{pedigree} = $fam_file_path;    #Pedigree files for samples
    }
    return %command;
}

sub write_cmd_mip_log {

##write_cmd_mip_log

##Function : Write CMD to MIP log file
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $order_parameters_ref, $script_ref, $log_file_ref
##         : $parameter_href        => Parameter hash {REF}
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $order_parameters_ref  => Order of addition to parameter array {REF}
##         : $script_ref            => The script that is being executed {REF}
##         : $log_file_ref          => The log file {REF}
##         : $mip_version_ref       => The MIP version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $order_parameters_ref;
    my $script_ref;
    my $log_file_ref;
    my $mip_version_ref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        order_parameters_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$order_parameters_ref
        },
        script_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$script_ref
        },
        log_file_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$log_file_ref
        },
        mip_version_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$mip_version_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $cmd_line = $$script_ref . " ";

    my @nowrite = (
        "mip",                  "bwa_build_reference",
        "pbamcalibrationblock", "pvariantannotationblock",
        q{associated_program},
    );

  PARAMETER_KEY:
    foreach my $order_parameter_element ( @{$order_parameters_ref} ) {

        if ( defined $active_parameter_href->{$order_parameter_element} ) {

            ## If no config file do not print
            if (   $order_parameter_element eq q{config_file}
                && $active_parameter_href->{config_file} eq 0 )
            {
            }
            else {

                ## If element is part of array - do nothing
                if ( ( any { $_ eq $order_parameter_element } @nowrite ) ) {
                }
                elsif (
                    ## Array reference
                    (
                        exists $parameter_href->{$order_parameter_element}
                        {data_type}
                    )
                    && ( $parameter_href->{$order_parameter_element}{data_type}
                        eq q{ARRAY} )
                  )
                {

                    my $separator = $parameter_href->{$order_parameter_element}
                      {element_separator};
                    $cmd_line .= "-"
                      . $order_parameter_element . " "
                      . join(
                        $separator,
                        @{
                            $active_parameter_href->{$order_parameter_element}
                        }
                      ) . " ";
                }
                elsif (
                    ## HASH reference
                    (
                        exists $parameter_href->{$order_parameter_element}
                        {data_type}
                    )
                    && ( $parameter_href->{$order_parameter_element}{data_type}
                        eq q{HASH} )
                  )
                {

                    # First key
                    $cmd_line .= "-" . $order_parameter_element . " ";
                    $cmd_line .= join(
                        "-" . $order_parameter_element . " ",
                        map {
"$_=$active_parameter_href->{$order_parameter_element}{$_} "
                          } (
                            keys %{
                                $active_parameter_href
                                  ->{$order_parameter_element}
                            }
                          )
                    );
                }
                else {

                    $cmd_line .= "-"
                      . $order_parameter_element . " "
                      . $active_parameter_href->{$order_parameter_element}
                      . " ";
                }
            }
        }
    }
    $log->info( $cmd_line,                            "\n" );
    $log->info( q{MIP Version: } . $$mip_version_ref, "\n" );
    $log->info(
        q{Script parameters and info from }
          . $$script_ref
          . q{ are saved in file: }
          . $$log_file_ref,
        "\n"
    );
    return;
}

sub check_unique_array_element {

##check_unique_array_element

##Function : Detects if there are items in query_ref that are not present in array_to_check_ref. If unique adds the unique element to array_to_check_ref.
##Returns  : ""
##Arguments: $array_to_check_ref, $array_query_ref
##         : $array_to_check_ref => The arrayref to be queried {REF}
##         : $query_ref          => The query reference can be either array or scalar {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $array_to_check_ref;
    my $query_ref;

    my $tmpl = {
        array_to_check_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$array_to_check_ref
        },
        query_ref => { required => 1, defined => 1, store => \$query_ref },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $array_query_ref;

    if ( ref($query_ref) eq "ARRAY" ) {    #Array reference

        $array_query_ref = $query_ref;
    }
    if ( ref($query_ref) eq "SCALAR" ) {    #Scalar reference

        push( @{$array_query_ref}, $$query_ref );    #Standardize to array
    }

    ##For each array_query_ref element, loop through corresponding array_to_check_ref element(s), add if there are none or an updated/unique entry.
    foreach my $query (@$array_query_ref) {

        if ( !( any { $_ eq $query } @$array_to_check_ref ) )
        {    #If element is not part of array

            push( @$array_to_check_ref, $query );    #Go ahead and add
        }
    }
}

sub determine_nr_of_rapid_nodes {

##determine_nr_of_rapid_nodes

##Function : Determines the number of nodes to allocate depending on the sequence read length, which affects the infile size.
##Returns  : $number_nodes, $read_nr_of_lines
##Arguments: $seq_length, $infile_size
##         : $seq_length  => Length of sequence reads
##         : $infile_size => Size of the infile

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $seq_length;
    my $infile_size;

    my $tmpl = {
        seq_length => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$seq_length
        },
        infile_size => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_size
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $number_nodes         = 0;  #Nodes to allocate
    my $read_position_weight = 1;  #Scales the read_start and read_stop position
    my $read_nr_of_lines;

    if ( $seq_length > 75 && $seq_length <= 101 ) {

        $read_nr_of_lines = 190000000;    #Read batch size
        $number_nodes = floor( $infile_size / ( 12 * $read_nr_of_lines ) )
          ; #Determines the number of nodes to use, 150000000 ~ 37,5 million reads, 13 = 2 sdtdev from sample population - currently poor estimate with compression confunding calculation.
        $log->info( "Number of Nodes: " . $number_nodes, "\n" );
    }
    if ( $seq_length > 50 && $seq_length <= 75 ) {

        $read_nr_of_lines = 190000000;    #Read batch size
        $number_nodes = floor( $infile_size / ( 9.75 * $read_nr_of_lines ) )
          ; #Determines the number of nodes to use, 150000000 ~ 37,5 million reads, 13 = 2 sdtdev from sample population - currently poor estimate with compression confunding calculation.
        $log->info( "Number of Nodes: " . $number_nodes, "\n" );
    }
    if ( $seq_length >= 50 && $seq_length < 75 ) {

        $read_nr_of_lines = 130000000;    #Read batch size
        $number_nodes = floor( $infile_size / ( 7 * $read_nr_of_lines ) )
          ; #Determines the number of nodes to use, 150000000 ~ 37,5 million reads, 13 = 2 sdtdev from sample population - currently poor estimate with compression confunding calculation.
        $log->info( "Number of Nodes: " . $number_nodes, "\n" );
    }
    if ( $seq_length >= 35 && $seq_length < 50 ) {

        $read_nr_of_lines = 95000000;    #Read batch size
        $number_nodes = floor( $infile_size / ( 6 * $read_nr_of_lines ) )
          ; #Determines the number of nodes to use, 150000000 ~ 37,5 million reads, 13 = 2 sdtdev from sample population - currently poor estimate with compression confunding calculation.
        $log->info( "Number of Nodes: " . $number_nodes, "\n" );
    }
    if ( $number_nodes <= 1 ) {

        $number_nodes = 2;    #Ensure that at least 1 readbatch is processed
    }
    return $number_nodes, $read_nr_of_lines;
}

sub check_unique_ids {

##check_unique_ids

##Function : Test that the family_id and the sample_id(s) exists and are unique. Check if id sample_id contains "_".
##Returns  : ""
##Arguments: $active_parameter_href, $sample_ids_ref
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $sample_ids_ref        => Array to loop in for parameter {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_ids_ref;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_ids_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$sample_ids_ref
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my %seen;    #Hash to test duplicate sample_ids later

    if ( !@$sample_ids_ref ) {

        $log->fatal("Please provide sample_id(s)\n");
        exit 1;
    }

    foreach my $sample_id (@$sample_ids_ref) {

        $seen{$sample_id}++;    #Increment instance to check duplicates later

        if ( $$family_id_ref eq $sample_id )
        {                       #Family_id cannot be the same as sample_id

            $log->fatal( "Family_id: "
                  . $$family_id_ref
                  . " equals sample_id: "
                  . $sample_id
                  . ". Please make sure that the family_id and sample_id(s) are unique.\n"
            );
            exit 1;
        }
        if ( $seen{$sample_id} > 1 ) {    #Check sample_id are unique

            $log->fatal( "Sample_id: " . $sample_id . " is not uniqe.\n" );
            exit 1;
        }
        if ( $sample_id =~ /_/ )
        { #Sample_id contains "_", which is not allowed according to filename conventions

            $log->fatal( "Sample_id: "
                  . $sample_id
                  . " contains '_'. Please rename sample_id according to MIP's filename convention, removing the '_'.\n"
            );
            exit 1;
        }
    }
}

sub update_config_file {

##update_config_file

##Function : Updates the config file to particular user/cluster for entries following specifications. Leaves other entries untouched.
##Returns  : ""
##Arguments: $active_parameter_href, $parameter_name_ref, $family_id_ref
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $parameter_name_ref    => MIP Parameter to update {REF}
##         : $family_id_ref         => Sets the family_id {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_name_ref;
    my $family_id_ref;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        parameter_name_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$parameter_name_ref
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( $active_parameter_href->{$$parameter_name_ref} ) {    #Active parameter

        if ( defined( $active_parameter_href->{cluster_constant_path} ) )
        {    #Set the project specific path for this cluster

            $active_parameter_href->{$$parameter_name_ref} =~
s/cluster_constant_path!/$active_parameter_href->{cluster_constant_path}/gi
              ;    #Exchange cluster_constant_path! for current cluster path
        }
        if ( defined( $active_parameter_href->{analysis_constant_path} ) )
        {          #Set the project specific path for this cluster

            $active_parameter_href->{$$parameter_name_ref} =~
s/analysis_constant_path!/$active_parameter_href->{analysis_constant_path}/gi
              ;  #Exchange analysis_constant_path! for the current analysis path
        }
        if ( defined($$family_id_ref) ) {    #Set the family_id

            $active_parameter_href->{$$parameter_name_ref} =~
              s/family_id!/$$family_id_ref/gi
              ;    #Exchange FND! for the current family_id
        }
        if ( defined( $active_parameter_href->{outaligner_dir} ) )
        {          #Set the outaligner_dir used

            $active_parameter_href->{$$parameter_name_ref} =~
              s/aligner_outdir!/$active_parameter_href->{outaligner_dir}/gi
              ;    #Exchange aligner_outdir! for the current outaligner_dir
        }
    }
}

sub check_auto_build {

##check_auto_build

##Function : Checks if autobuild is on and returns "1" if enabled or "0" if not
##Returns  : "0|1"
##Arguments: $parameter_href, $active_parameter_href, $parameter_name_ref, $sample_id_ref
##         : $parameter_href        => Parameter hash {REF}
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $parameter_name_ref    => MIP parameter name {REF}
##         : $sample_id_ref         => SampleId {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $parameter_name_ref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        parameter_name_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$parameter_name_ref
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( $parameter_href->{$$parameter_name_ref}{build_file} eq
        "yes_auto_build" )
    {

        return "1";    #Flag that autobuild is needed
    }
    elsif ( $parameter_href->{$$parameter_name_ref}{build_file} eq 1 )
    {                  #1 for arrays

        return "1";    #Flag that autobuild is needed
    }
    else {

        return "0";    #No autobuild is needed
    }
}

sub parse_human_genome_reference {

##parse_human_genome_reference

##Function : Detect version and source of the human_genome_reference: Source (hg19 or GRCh).
##Returns  : ""
##Arguments: $file_info_href, $human_genome_reference_ref
##         : $file_info_href             => File info hash {REF}
##         : $human_genome_reference_ref => The human genome {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $human_genome_reference_ref;

    my $tmpl = {
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        human_genome_reference_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$human_genome_reference_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    if ( $$human_genome_reference_ref =~ /GRCh(\d+\.\d+|\d+)_homo_sapiens_/ )
    {    #Used to change capture kit genome reference version later

        $file_info_href->{human_genome_reference_version} = $1;
        $file_info_href->{human_genome_reference_source}  = "GRCh";    #Ensembl
    }
    elsif ( $$human_genome_reference_ref =~ /hg(\d+)_homo_sapiens/ )
    {    #Used to change capture kit genome reference version later

        $file_info_href->{human_genome_reference_version} = $1;
        $file_info_href->{human_genome_reference_source}  = "hg";    #Refseq
    }
    else {

        $log->warn(
"MIP cannot detect what kind of human_genome_reference you have supplied. If you want to automatically set the capture kits used please supply the reference on this format: [source]_[species]_[version].",
            "\n"
        );
    }

    ## Removes ".file_ending" in filename.FILENDING(.gz)
    $file_info_href->{human_genome_reference_name_prefix} =
      fileparse( $$human_genome_reference_ref, qr/\.fasta|\.fasta\.gz/ );

    $file_info_href->{human_genome_compressed} =
      check_gzipped( { file_name_ref => $human_genome_reference_ref, } );
}

sub check_file_endings_to_build {

##check_file_endings_to_build

##Function : Checks files to be built by combining filename stub with fileendings.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, file_endings_ref, $parameter_name
##         : $parameter_href        => Parameter hash {REF}
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $file_endings_ref      => Reference to the file_endings to be added to the filename stub {REF}
##         : $parameter_name        => MIP parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $file_endings_ref;
    my $parameter_name;
    my $file_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        file_endings_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$file_endings_ref
        },
        parameter_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_name,
        },
        file_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_name
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    foreach my $file_ending (@$file_endings_ref) {

        check_existance(
            {
                parameter_href        => $parameter_href,
                active_parameter_href => $active_parameter_href,
                item_name_ref         => \catfile( $file_name . $file_ending ),
                parameter_name_ref    => \$parameter_name,
                item_type_to_check    => "file",
            }
        );
    }
}

sub check_existance {

##check_existance

##Function : Checks if a file/directory exists and if auto_build is on or not. If file/directory does not extis and there is no autobuild, croaks and exists.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $item_name_ref, $parameter_name_ref, $item_type_to_check, $temp_directory_ref
##         : $parameter_href        => The parameters hash
##         : $active_parameter_href => The active parameter for this analysis hash
##         : $item_name_ref         => Item to check for existance {REF}
##         : $parameter_name_ref    => MIP parameter name {REF}
##         : $item_type_to_check    => The type of item to check
##         : $temp_directory_ref    => Temporary directory

    my ($arg_href) = @_;

    ## Default(s)
    my $temp_directory_ref;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $item_name_ref;
    my $parameter_name_ref;
    my $item_type_to_check;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        item_name_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$item_name_ref
        },
        parameter_name_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$parameter_name_ref
        },
        item_type_to_check => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$item_type_to_check
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    if ( $item_type_to_check eq "directory" ) {

        unless ( -d $$item_name_ref ) {   #Check existence of supplied directory

            $log->fatal( $USAGE, "\n" );
            $log->fatal(
                "Could not find intended "
                  . $$parameter_name_ref
                  . " directory: "
                  . $$item_name_ref,
                "\n"
            );
            exit 1;
        }
    }
    elsif ( $item_type_to_check eq "file" ) {

        unless ( -f $$item_name_ref ) {    #Check existence of supplied file

            ## Check auto_build or not and return value
            $parameter_href->{$$parameter_name_ref}{build_file} =
              check_auto_build(
                {
                    parameter_href        => $parameter_href,
                    active_parameter_href => $active_parameter_href,
                    parameter_name_ref    => $parameter_name_ref,
                }
              );

            if ( $parameter_href->{$$parameter_name_ref}{build_file} == 0 )
            {    #No autobuild

                $log->fatal( $USAGE, "\n" );
                $log->fatal(
                    "Could not find intended "
                      . $$parameter_name_ref
                      . " file: "
                      . $$item_name_ref,
                    "\n"
                );
                exit 1;
            }
        }
        else {

            if (
                (
                    defined(
                        $parameter_href->{$$parameter_name_ref}{build_file}
                    )
                )
                && ( $parameter_href->{$$parameter_name_ref}{build_file} ne 1 )
              )
            {   #If any of associated files do not exist make sure to build them

                $parameter_href->{$$parameter_name_ref}{build_file} =
                  0;    #File exist in this check
            }
        }
    }
}

sub check_user_supplied_info {

##check_user_supplied_info

##Function : Determine if the user supplied info on parameter either via cmd or config
##Returns  : "0|1"
##Arguments: $active_parameter_href, $data_ref, $parameter_name
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $data_ref              => Data to check for existence {REF}
##         : $parameter_name        => MIP parameter to evaluate

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $data_ref;
    my $parameter_name;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        data_ref => { required => 1, defined => 1, store => \$data_ref },
        parameter_name => { strict_type => 1, store => \$parameter_name },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $user_supplied_info_switch;

    if ( ref($data_ref) eq "ARRAY" ) {    #Array reference

        if ( !@$data_ref ) {              #No user supplied sample info

            if ( defined( $active_parameter_href->{$parameter_name} ) )
            {                             #User supplied info in config file

                $user_supplied_info_switch = 0
                  ; #No user supplied cmd info, but present in config file overwrite using info from pedigree file
            }
            else {    #No sample_ids info in config file

                $user_supplied_info_switch = 0
                  ; #No user supplied cmd info, not defined in config file, ADD it from pedigree file
            }
        }
        else {

            $user_supplied_info_switch = 1
              ; #User supplied cmd info, do NOT overwrite using info from pedigree file
        }
    }
    elsif ( ref($data_ref) eq "HASH" ) {

        if ( !%$data_ref ) {

            if ( defined( $active_parameter_href->{$parameter_name} ) )
            {    #User supplied info in config file

                $user_supplied_info_switch = 0
                  ; #No user supplied cmd info, but present in config file overwrite using info from pedigree file
            }
            else {    #No sample_ids info in config file

                $user_supplied_info_switch = 0
                  ; #No user supplied cmd info, not defined in config file, ADD it from pedigree file
            }
        }
        else {

            $user_supplied_info_switch = 1
              ; #User supplied cmd info, do NOT overwrite using info from pedigree file
        }
    }
    return $user_supplied_info_switch;
}

sub check_gzipped {

##check_gzipped

##Function : Check if a file is gzipped.
##Returns  : "0 (=uncompressed)| 1 (=compressed)"
##Arguments: $file_name_ref
##         : $file_name_ref => File name {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_name_ref;

    my $tmpl = {
        file_name_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$file_name_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $file_compression_status = 0;

    if ( $$file_name_ref =~ /.gz$/ ) {

        $file_compression_status = 1;
    }
    return $file_compression_status;
}

sub check_target_bed_file_exist {

##check_target_bed_file_exist

##Function : Check that supplied target file ends with ".bed" and otherwise exists.
##Returns  : ""
##Arguments: $active_parameter_href, $file, $parameter_name
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $file                  => File to check for existance and file ending
##         : $parameter_name        => MIP parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file;
    my $parameter_name;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        file =>
          { required => 1, defined => 1, strict_type => 1, store => \$file },
        parameter_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_name,
        },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    if ( $file !~ /.bed$/ ) {

        $log->fatal(
            "Could not find intendended '.bed file ending' for target file: "
              . $file
              . " in parameter '-"
              . $parameter_name . "'",
            "\n"
        );
        exit 1;
    }
}

sub compare_array_elements {

##compare_array_elements

##Function : Compares the number of elements in two arrays and exits if the elements are not equal.
##Returns  : ""
##Arguments: $elements_ref, $array_queries_ref, $parameter_name, $parameter_name_query
##         : $elements_ref         => Array to match {REF}
##         : $array_queries_ref    => Array to be compared {REF}
##         : $parameter_name       => MIP reference parameter
##         : $parameter_name_query => MIP query parameter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $elements_ref;
    my $array_queries_ref;
    my $parameter_name;
    my $parameter_name_query;

    my $tmpl = {
        elements_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$elements_ref
        },
        array_queries_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$array_queries_ref
        },
        parameter_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_name,
        },
        parameter_name_query => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_name_query
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    if ( scalar(@$elements_ref) != scalar(@$array_queries_ref) ) {

        $log->fatal(
            "The number of supplied '-"
              . $parameter_name_query . "' (="
              . scalar(@$array_queries_ref)
              . ") do not equal the number of '-"
              . $parameter_name . "' (="
              . scalar(@$elements_ref)
              . "). Please specify a equal number of elements for both parameters",
            "\n"
        );
        exit 1;
    }
}

sub collect_seq_contigs {

##collect_seq_contigs

##Function : Collects sequences contigs used in analysis from human genome sequence dictionnary associated with $human_genome_reference
##Returns  : ""
##Arguments: $contigs_ref, $reference_dir_ref, $human_genome_reference_name_prefix_ref
##         : $contigs_ref                               => Contig array {REF}
##         : $reference_dir_ref                         => The MIP reference directory {REF}
##         : $human_genome_reference_name_prefix_ref => The associated human genome file without file ending

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contigs_ref;
    my $reference_dir_ref;
    my $human_genome_reference_name_prefix_ref;

    my $tmpl = {
        contigs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$contigs_ref
        },
        reference_dir_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$reference_dir_ref,
        },
        human_genome_reference_name_prefix_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$human_genome_reference_name_prefix_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $pqSeqDict =
q?perl -nae 'if($F[0]=~/^\@SQ/) { if($F[1]=~/SN\:(\S+)/) {print $1, ",";} }' ?;
    my $SeqDictLocation = catfile( $$reference_dir_ref,
        $$human_genome_reference_name_prefix_ref . ".dict" );
    @$contigs_ref = `$pqSeqDict $SeqDictLocation `
      ;    #Returns a comma seperated string of sequence contigs from dict file
    @$contigs_ref = split( /,/, join( ',', @$contigs_ref ) );
}

sub collect_select_file_contigs {

##collect_select_file_contigs

##Function : Collects sequences contigs used in select file
##Returns  : ""
##Arguments: $contigs_ref, $select_file_path
##         : $contigs_ref      => Contig array {REF}
##         : $select_file_path => The select file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contigs_ref;
    my $select_file_path;

    my $tmpl = {
        contigs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$contigs_ref
        },
        select_file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$select_file_path
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $pquery_seq_dict =
q?perl -nae 'if ($_=~/contig\=(\w+)/) {print $1, ",";} if($_=~/#CHROM/) {last;}' ?;
    @$contigs_ref = `$pquery_seq_dict $select_file_path `
      ;    #Returns a comma seperated string of sequence contigs from file
    @$contigs_ref = split( /,/, join( ',', @$contigs_ref ) );

    if ( !@$contigs_ref ) {

        $log->fatal(
"Could not detect any '##contig' in meta data header in select file: "
              . $select_file_path
              . "\n" );
        exit 1;
    }
}

sub size_sort_select_file_contigs {

##size_sort_select_file_contigs

##Function : Sorts array depending on reference array. NOTE: Only entries present in reference array will survive in sorted array.
##Returns  : "@sorted_contigs"
##Arguments: $file_info_href, $consensus_analysis_type_ref, $hash_key_to_sort, $hash_key_sort_reference
##         : $file_info_href              => File info hash {REF}
##         : $consensus_analysis_type_ref => Consensus analysis_type {REF}
##         : $hash_key_to_sort            => The keys to sort
##         : $hash_key_sort_reference     => The hash keys sort reference

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $consensus_analysis_type_ref;
    my $hash_key_to_sort;
    my $hash_key_sort_reference;

    my $tmpl = {
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        consensus_analysis_type_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$consensus_analysis_type_ref
        },
        hash_key_to_sort => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$hash_key_to_sort
        },
        hash_key_sort_reference => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$hash_key_sort_reference
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my @sorted_contigs;

    ##Sort the contigs depending on reference array
    if ( $file_info_href->{$hash_key_to_sort} ) {

        foreach my $element ( @{ $file_info_href->{$hash_key_sort_reference} } )
        {

            if (
                !check_entry_hash_of_array(
                    {
                        hash_ref => $file_info_href,
                        key      => $hash_key_to_sort,
                        element  => $element,
                    }
                )
              )
            {

                push( @sorted_contigs, $element );
            }
        }
    }

    ## Test if all contigs collected from select file was sorted by reference contig array
    if (
        (@sorted_contigs)
        && (
            scalar( @{ $file_info_href->{$hash_key_to_sort} } ) !=
            scalar(@sorted_contigs) )
      )
    {

        foreach my $element ( @{ $file_info_href->{$hash_key_to_sort} } ) {

            if ( !( any { $_ eq $element } @sorted_contigs ) )
            {    #If element is not part of array

                unless ( ( $$consensus_analysis_type_ref eq "wes" )
                    && ( $element =~ /MT$|M$/ ) )
                { #Special case when analysing wes since Mitochondrial contigs have no baits in exome capture kits

                    $log->fatal( "Could not detect '##contig'= "
                          . $element
                          . " from meta data header in '-vcfparser_select_file' in reference contigs collected from '-human_genome_reference'\n"
                    );
                    exit 1;
                }
            }
        }
    }
    return @sorted_contigs;
}

sub replace_config_parameters_with_cmd_info {

##replace_config_parameters_with_cmd_info

##Function : Replace config parameter with cmd info for config dynamic parameter
##Returns  :
##Arguments: $parameter_href, $active_parameter_href, $parameter_names_ref
##         : $parameter_href        => Parameter hash {REF}
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $parameter_names_ref   => MIP activate parameter names {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $parameter_names_ref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        parameter_names_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$parameter_names_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    foreach my $parameter_name (@$parameter_names_ref) {

        if ( defined( $parameter_href->{$parameter_name}{value} ) )
        {    #Replace config parameter with cmd info for parameter

            $active_parameter_href->{$parameter_name} =
              $parameter_href->{$parameter_name}{value}
              ;    #Transfer to active parameter
        }
        elsif (( exists( $parameter_href->{$parameter_name}{default} ) )
            && ( !defined( $active_parameter_href->{$parameter_name} ) ) )
        {

            $active_parameter_href->{$parameter_name} =
              $parameter_href->{$parameter_name}{default}
              ;    #Transfer to active parameter
        }
    }
}

sub check_entry_hash_of_array {

##check_entry_hash_of_array

##Function : Test element for being part of hash of array at supplied key.
##Returns  : Return "1" if element is not part of array
##Arguments: $hash_ref, $key, $element
##         : $hash_ref => Hash {REF}
##         : $key      => The key pointing to the array in the $hash_ref
##         : $element  => Element to look for in hash of array

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $hash_ref;
    my $key;
    my $element;

    my $tmpl = {
        hash_ref => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$hash_ref
        },
        key =>
          { required => 1, defined => 1, strict_type => 1, store => \$key },
        element =>
          { required => 1, defined => 1, strict_type => 1, store => \$element },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( defined( $$hash_ref{$key} ) ) {    #Information on entry present

        if ( !( any { $_ eq $element } @{ $hash_ref->{$key} } ) )
        {                                   #If element is not part of array

            return 1;
        }
    }
}

sub check_most_complete_and_remove_file {

##check_most_complete_and_remove_file

##Function  : Checks if the file is recorded as the "most_complete_bam|vcf". If false writes removal of file(s) to supplied filehandle
##Returns   : ""
##Arguments : $FILEHANDLE, $file_path_ref, $file_ending, $most_complete_ref
##          : $FILEHANDLE        => SBATCH script FILEHANDLE to print to
##          : $file_path_ref     => Current file {REF}
##          : $file_ending       => File ending of $file_path_ref
##          : $most_complete_ref => The mostComplete file (BAM|VCF) {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $most_complete_ref;
    my $FILEHANDLE;
    my $file_path_ref;
    my $file_ending;

    my $tmpl = {
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE, },
        file_path_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$file_path_ref
        },
        file_ending => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_ending
        },
        most_complete_ref => { store => \$most_complete_ref },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Coreutils qw(gnu_rm);

    if ( ( defined($$most_complete_ref) ) && ( defined($$file_path_ref) ) )
    {    #Not to disturb first dry_run of analysis

        unless ( $$most_complete_ref eq $$file_path_ref )
        {    #Do not remove mostCompleteBAM|VCF

            ## Modify fileending of file to include e.g. .bai for bams
            my $file_name = modify_file_ending(
                {
                    file_path_ref => $file_path_ref,
                    file_ending   => $file_ending,
                }
            );

            ##Print removal of file to sbatch script
            gnu_rm(
                {
                    infile_path => $file_name,
                    force       => 1,
                    FILEHANDLE  => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} "\n";    #Remove file(s)
        }
    }
    else {

        ## Modify fileending of file to include e.g. .bai for bams
        my $file_name = modify_file_ending(
            {
                file_path_ref => $file_path_ref,
                file_ending   => $file_ending,
            }
        );

        ##Print removal of file to sbatch script
        gnu_rm(
            {
                infile_path => $file_name,
                force       => 1,
                FILEHANDLE  => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} "\n";    #Remove file(s)
    }
}

sub modify_file_ending {

##modify_file_ending

##Function  : Modify fileending of file to include e.g. .bai for bams
##Returns   : ""
##Arguments : $file_path_ref, $file_ending
##          : $file_path_ref => Current file {REF}
##          : $file_ending   => File ending of $file_path_ref

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path_ref;
    my $file_ending;

    my $tmpl = {
        file_path_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$file_path_ref
        },
        file_ending => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_ending
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Removes ".file_ending" in filename.FILENDING
    my ( $file_name, $dir_path ) = fileparse( $$file_path_ref, $file_ending );
    my $file_path_prefix = catfile( $dir_path, $file_name );

    if ( defined($file_path_prefix) ) {    #Successfully removed file ending

        my $end = ".*";    #Remove all files with ending with ".*"

        if ( $file_ending eq q{.bam} ) {    #For BAM files

            $end = ".ba*";                  #Removes both .bam and .bai
        }
        if ( $file_ending eq ".vcf" ) {     #For VCF files

            $end = ".vcf*";                 #Removes both .vcf and .vcf.idx
        }
        return $file_path_prefix . $end;
    }
    else {

        return $$file_path_ref;
    }
}

sub concatenate_vcfs {

##concatenate_vcfs

##Function : Concatenate VCFs
##Returns  : ""
##Arguments: $active_parameter_href, $FILEHANDLE, $arrays_ref, $infile_prefix, $infile_postfix, $outfile, $reorder_swith
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $FILEHANDLE            => SBATCH script FILEHANDLE to print to
##         : $arrays_ref            => Holding the number and part of file names to be combined
##         : $infile_prefix         => Will be combined with the each array element
##         : $infile_postfix        => Will be combined with the each array element
##         : $outfile               => The combined outfile
##         : $reorder_swith         => Reorder header

    my $active_parameter_href = $_[0];
    my $FILEHANDLE            = $_[1];
    my $arrays_ref            = $_[2];
    my $infile_prefix         = $_[3];
    my $infile_postfix        = $_[4];
    my $outfile               = $_[5];
    my $reorder_swith         = $_[6];

    use MIP::Language::Java qw{java_core};

    unless ( defined($infile_postfix) ) {

        $infile_postfix = "";    #No postfix
    }
    unless ( defined($outfile) ) {

        $outfile = $infile_prefix . ".vcf";
    }

    ## Writes java core commands to filehandle.
    java_core(
        {
            FILEHANDLE        => $FILEHANDLE,
            memory_allocation => "Xmx4g",
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $active_parameter_href->{temp_directory},
            java_jar =>
              catfile( $active_parameter_href->{snpeff_path}, "SnpSift.jar" ),
        }
    );

    print {$FILEHANDLE} "split -j ";    #Joinf VCFs together

    foreach my $element (@$arrays_ref) {

        print {$FILEHANDLE} $infile_prefix
          . $element
          . $infile_postfix
          . " ";                        #files to combined
    }
    if ( ( defined( $_[6] ) ) && $reorder_swith eq "reorder_header" ) {

        print {$FILEHANDLE} "| ";            #Pipe
        print {$FILEHANDLE} "vcfparser ";    #Parses the vcf output
    }

    print {$FILEHANDLE} "> " . $outfile;     #OutFile
}

sub concatenate_variants {

##concatenate_variants

##Function : Writes sbatch code to supplied filehandle to concatenate variants in vcf format. Each array element is combined with the infile prefix and postfix.
##Returns  : ""
##Arguments: $active_parameter_href, $FILEHANDLE, $elements_ref, $infile_prefix, $infile_postfix, $outfile, $human_genome_reference_ref
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $FILEHANDLE                 => SBATCH script FILEHANDLE to print to
##         : $elements_ref               => Holding the number and part of file names to be combined
##         : $infile_prefix              => Will be combined with the each array element
##         : $infile_postfix             => Will be combined with the each array element
##         : $outfile                    => The combined outfile
##         : $human_genome_reference_ref => Human genome reference {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $human_genome_reference_ref;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $elements_ref;
    my $FILEHANDLE;
    my $infile_prefix;
    my $infile_postfix;
    my $outfile;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        elements_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$elements_ref
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE, },
        infile_prefix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_prefix
        },
        infile_postfix => { strict_type => 1, store => \$infile_postfix },
        outfile        => { strict_type => 1, store => \$outfile, },
        human_genome_reference_ref => {
            default =>
              \$arg_href->{active_parameter_href}{human_genome_reference},
            strict_type => 1,
            store       => \$human_genome_reference_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Java qw{java_core};
    use MIP::Program::Variantcalling::Gatk qw(gatk_catvariants);

    unless ( defined($infile_postfix) ) {

        $infile_postfix = "";    #No postfix
    }
    unless ( defined($outfile) ) {

        $outfile = $infile_prefix . ".vcf";
    }

    say {$FILEHANDLE} "## GATK CatVariants";

    ## Assemble infile paths
    my @infile_paths =
      map { $infile_prefix . $_ . $infile_postfix } @$elements_ref;

    gatk_catvariants(
        {
            memory_allocation => "Xmx4g",
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $active_parameter_href->{temp_directory},
            gatk_path      => catfile( $active_parameter_href->{gatk_path},
                "GenomeAnalysisTK.jar" )
              . " org.broadinstitute.gatk.tools.CatVariants",
            infile_paths_ref   => \@infile_paths,
            outfile_path       => $outfile,
            referencefile_path => $$human_genome_reference_ref,
            logging_level      => $active_parameter_href->{gatk_logging_level},
            FILEHANDLE         => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";
}

sub sort_vcf {

##sort_vcf

##Function : Writes sbatch code to supplied filehandle to sort variants in vcf format.
##Returns  : ""
##Arguments: $active_parameter_href, $infile_paths_ref, $FILEHANDLE, $sequence_dict_file, $outfile
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $infile_paths_ref      => Infiles to sort {REF}
##         : $FILEHANDLE            => SBATCH script FILEHANDLE to print to
##         : $sequence_dict_file    => Human reference sequence dict file
##         : $outfile               => The sorted outfile

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $infile_paths_ref;
    my $FILEHANDLE;
    my $sequence_dict_file;
    my $outfile;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE, },
        sequence_dict_file => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sequence_dict_file
        },
        outfile => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Variantcalling::Picardtools qw(picardtools_sortvcf);

    say {$FILEHANDLE} "## Picard SortVcf";

    ## Writes java core commands to filehandle.
    picardtools_sortvcf(
        {
            infile_paths_ref    => \@$infile_paths_ref,
            outfile_path        => $outfile,
            sequence_dictionary => $sequence_dict_file,
            memory_allocation   => q{Xmx2g},
            referencefile_path =>
              $active_parameter_href->{human_genome_reference},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $active_parameter_href->{temp_directory},
            java_jar       => catfile(
                $active_parameter_href->{picardtools_path},
                q{picard.jar}
            ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";
}

sub detect_sample_id_gender {

##detect_sample_id_gender

##Function : Detect gender of the current analysis
##Returns  : "$found_male $found_female $found_other"
##Arguments: $active_parameter_href, $sample_info_href
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $sample_info_href      => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $found_male   = 0;
    my $found_female = 0;
    my $found_other  = 0;

    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        if ( $sample_info_href->{sample}{$sample_id}{sex} =~ /1|^male/ ) { #Male

            $found_male = 1;                                               #Male
        }
        elsif ( $sample_info_href->{sample}{$sample_id}{sex} =~ /2|female/ )
        {    #Female

            $found_female = 1;
        }
        else {    #Other

            $found_male =
              1;    #Include since it might be male to enable analysis of Y.
            $found_other = 1;
        }
    }
    return $found_male, $found_female, $found_other;
}

sub remove_pedigree_elements {

##remove_pedigree_elements

##Function : Removes ALL keys at third level except keys in allowed_entries hash.
##Returns  : ""
##Arguments: $hash_ref
##         : $hash_ref => Hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $hash_ref;

    my $tmpl = {
        hash_ref => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$hash_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @allowed_entries = (
        "family",        "default_gene_panels",
        "sample",        "sample_id",
        "sample_name",   "capture_kit",
        "sex",           "mother",
        "father",        "phenotype",
        "sequence_type", "expected_coverage",
    );

  FAMILY_INFO:
    for my $key ( keys %$hash_ref ) {

        if ( !any { $_ eq $key } @allowed_entries )
        {    #If element is not part of array

            delete( $hash_ref->{$key} );
        }
    }

  SAMPLE:
    foreach my $sample_id ( keys %{ $hash_ref->{sample} } ) {

      SAMPLE_INFO:
        for my $pedigree_element ( keys %{ $hash_ref->{sample}{$sample_id} } ) {

            if ( !any { $_ eq $pedigree_element } @allowed_entries )
            {    #If element is not part of array

                delete( $hash_ref->{sample}{$sample_id}{$pedigree_element} );
            }
        }
    }
}

sub deafult_log4perl_file {

##deafult_log4perl_file

##Function : Set the default Log4perl file using supplied dynamic parameters.
##Returns  : "$log_file"
##Arguments: $active_parameter_href, $cmd_input_ref, $script_ref, $date_ref, $date_time_stamp_ref, $family_id_ref, $outdata_dir_ref
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $cmd_input_ref         => User supplied info on cmd for log_file option {REF}
##         : $script_ref            => The script that is executed {REF}
##         : $date_ref              => The date {REF}
##         : $date_time_stamp_ref   => The date and time {REF}
##         : $family_id_ref         => Family id {REF}
##         : $outdata_dir_ref       => Outdata directory {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $outdata_dir_ref;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $cmd_input_ref;
    my $script_ref;
    my $date_ref;
    my $date_time_stamp_ref;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        cmd_input_ref =>
          { default => \$$, strict_type => 1, store => \$cmd_input_ref },
        script_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$script_ref
        },
        date_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$date_ref
        },
        date_time_stamp_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$date_time_stamp_ref
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        outdata_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outdata_dir},
            strict_type => 1,
            store       => \$outdata_dir_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    unless ( defined($$cmd_input_ref) )
    {   #No input from cmd i.e. create default logging directory and set default

        make_path( catfile( $$outdata_dir_ref, "mip_log", $$date_ref ) );
        my $log_file = catfile( $$outdata_dir_ref, "mip_log", $$date_ref,
            $$script_ref . "_" . $$date_time_stamp_ref . ".log" )
          ;    #concatenates log filename
        return $log_file;
    }
}

sub check_human_genome_file_endings {

##check_human_genome_file_endings

##Function : Check the existance of associated Human genome files.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $file_info_href, $parameter_name_ref, $human_genome_reference_name_prefix_ref, $reference_dir_ref,
##         : $parameter_href                            => Parameter hash {REF}
##         : $active_parameter_href                     => Holds all set parameter for analysis {REF}
##         : $file_info_href                            => File info hash {REF}
##         : $parameter_name_ref                        => The parameter under evaluation {REF}
##         : $human_genome_reference_name_prefix_ref => The associated human genome file without file ending {REF}
##         : $reference_dir_ref                         => MIP reference directory {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $reference_dir_ref;
    my $human_genome_reference_name_prefix_ref;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $file_info_href;
    my $parameter_name_ref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        parameter_name_ref =>
          { default => \$$, strict_type => 1, store => \$parameter_name_ref },
        reference_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{reference_dir},
            strict_type => 1,
            store       => \$reference_dir_ref,
        },
        human_genome_reference_name_prefix_ref => {
            default =>
              \$arg_href->{file_info_href}{human_genome_reference_name_prefix},
            strict_type => 1,
            store       => \$human_genome_reference_name_prefix_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    foreach my $file_ending (
        @{ $file_info_href->{human_genome_reference_file_endings} } )
    {

        my $path = catfile( $active_parameter_href->{human_genome_reference} );

        ## Enable auto_build of metafiles
        $parameter_href->{ $$parameter_name_ref . $file_ending }{build_file} =
          "yes_auto_build";

        if ( $file_ending eq ".dict" ) {

            ## Removes ".file_ending" in filename.FILENDING(.gz)
            my ( $file, $dir_path ) =
              fileparse( $active_parameter_href->{human_genome_reference},
                qr/\.fasta|\.fasta\.gz/ );
            $path = catfile( $dir_path, $file );
        }

        $path = $path . $file_ending;    #Add current ending
        my $complete_parameter_name = $$parameter_name_ref . $file_ending;

        check_existance(
            {
                parameter_href        => $parameter_href,
                active_parameter_href => $active_parameter_href,
                item_name_ref         => \$path,
                parameter_name_ref    => \$complete_parameter_name,
                item_type_to_check    => "file",
            }
        );
    }
    if ( $parameter_href->{ $$parameter_name_ref . ".dict" }{build_file} eq 0 )
    {

        ##Collect sequence contigs from human reference ".dict" file since it exists
        collect_seq_contigs(
            {
                contigs_ref       => \@{ $file_info_href->{contigs} },
                reference_dir_ref => $reference_dir_ref,
                human_genome_reference_name_prefix_ref =>
                  $human_genome_reference_name_prefix_ref,
            }
        );
    }
}

sub collect_path_entries {

##collect_path_entries

##Function  : Collects all programs outfile path(s) created by MIP as Path->value located in %sample_info.
##Returns   : ""
##Arguments: $sample_info_href, $paths_ref
##         : $sample_info_href => Info on samples and family hash {REF}
##         : $paths_ref        => Holds the collected paths {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $paths_ref;

    my $tmpl = {
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$paths_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %info =
      %$sample_info_href;    #Copy hash to enable recursive removal of keys

    my @outdirectories
      ;   #Temporary array for collecting outDirectories within the same program
    my @outfileArray
      ;    #Temporary array for collecting outfile within the same program

    while ( my ( $key, $value ) = each %info ) {

        if ( ref($value) eq "HASH" ) {

            collect_path_entries(
                {
                    sample_info_href => $value,
                    paths_ref        => $paths_ref,
                }
            );
        }
        else {

            if ($value) {    #Required for first dry-run

                ## Check if key is "Path" and adds value to @paths_ref if true.
                check_and_add_to_array(
                    {
                        paths_ref => $paths_ref,
                        value     => $value,
                        key       => $key,
                    }
                );

                ## Check if key is "OutDirectory" or "OutFile"  and adds joined value to @paths_ref if true.
                collect_outfile(
                    {
                        paths_ref          => $paths_ref,
                        outdirectories_ref => \@outdirectories,
                        outfiles_ref       => \@outfileArray,
                        value              => $value,
                        key                => $key,
                    }
                );

                delete( $info{$value} );
            }
        }
    }
}

sub check_and_add_to_array {

##check_and_add_to_array

##Function  : Check if KeyName is "Path" and adds to @paths_ref if true.
##Returns   : ""
##Arguments: $paths_ref, $value, $key
##         : $paths_ref => Holds the collected paths {REF}
##         : $value     => Hash value
##         : $keyName   => Hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $paths_ref;
    my $value;
    my $key;

    my $tmpl = {
        paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$paths_ref
        },
        value => { required => 1, store => \$value },
        key =>
          { required => 1, defined => 1, strict_type => 1, store => \$key },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( $key eq "path" ) {

        if ( !( any { $_ eq $value } @$paths_ref ) )
        {    #Do not add same path twice

            push( @$paths_ref, $value );
        }
    }
}

sub collect_outfile {

##collect_outfile

##Function  : Check if KeyName is "OutDirectory" or "OutFile"  and adds to @paths_ref if true.
##Returns   : ""
##Arguments: $paths_ref, $outdirectories_ref, $outfiles_ref, $value, $key
##         : $paths_ref          => Holds the collected paths {REF}
##         : $outdirectories_ref => Holds temporary outdirectory path(s) {Optional, REF}
##         : $outfiles_ref       => Holds temporary outdirectory path(s) {Optional, REF}
##         : $value              => Hash value
##         : $key                => Hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $paths_ref;
    my $outdirectories_ref;
    my $outfiles_ref;
    my $value;
    my $key;

    my $tmpl = {
        paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$paths_ref
        },
        outdirectories_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$outdirectories_ref
        },
        outfiles_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$outfiles_ref
        },
        value => { required => 1, defined => 1, store => \$value },
        key =>
          { required => 1, defined => 1, strict_type => 1, store => \$key },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( $key eq "outdirectory" ) {

        push( @$outdirectories_ref, $value );
    }
    if ( $key eq "outfile" ) {

        push( @$outfiles_ref, $value );
    }
    if ( (@$outdirectories_ref) && (@$outfiles_ref) )
    {    #Both outdirectory and outfile have been collected, time to join

        my $path = catfile( $outdirectories_ref->[0], $outfiles_ref->[0] );

        if ( !( any { $_ eq $path } @$paths_ref ) )
        {    #Do not add same path twice

            push( @$paths_ref,
                catfile( $outdirectories_ref->[0], $outfiles_ref->[0] ) );
            @$outdirectories_ref = ();    #Restart
            @$outfiles_ref       = ();    #Restart
        }
    }
}

sub check_email_address {

##check_email_address

##Function : Check the syntax of the email adress is valid not that it is actually exists.
##Returns  : ""
##Arguments: $email_ref
##         : $email_ref => The email adress

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $email_ref;

    my $tmpl = {
        email_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$email_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    $$email_ref =~
      /[ |\t|\r|\n]*\"?([^\"]+\"?@[^ <>\t]+\.[^ <>\t][^ <>\t]+)[ |\t|\r|\n]*/;

    unless ( defined($1) ) {

        $log->fatal(
            "The supplied email: " . $$email_ref . " seem to be malformed. ",
            "\n" );
        exit 1;
    }
}

sub break_string {

##break_string

##Function : Breaks the string supplied on format key1:value1_value2,keyN:valueN_valueN,..n . Add key to %file_info_href and values as array. This enables association of values to key supplied in config or cmd.
##Returns  : ""
##Arguments: $file_info_href, $parameter_value_ref, $parameter_name, $associated_program
##         : $file_info_href      => File info hash {REF}
##         : $parameter_value_ref => MIP parameter value {REF}
##         : $parameter_name      => MIP parameter name {REF}
##         : $associated_program  => The parameters program {REF}

    my $file_info_href         = $_[0];
    my $parameter_value_ref    = $_[1];
    my $parameter_name_ref     = $_[2];
    my $associated_program_ref = $_[3];

    ## Break string into key value pairs
    my @file_keys = split( ',', join( ',', $$parameter_value_ref ) );

    foreach my $element (@file_keys) {

        my @temps = split( /:/, $element );
        @{ $file_info_href->{$$associated_program_ref}{$$parameter_name_ref}
              { $temps[0] } } =
          split( '_', $temps[1] );    #Save info_key associated with file_name
    }
}

sub add_capture_kit {

##add_capture_kit

##Function : Return a capture kit depending on user info. If arg->{user_supplied_parameter_switchRef} is set, go a head and add capture kit no matter what the switch was.
##Returns  : "Set capture kit or ''"
##Arguments: $file_info_href, $supported_capture_kit_href, $capture_kit, $user_supplied_parameter_switch
##         : $file_info_href                 => File info hash {REF}
##         : $supported_capture_kit_href     => The supported capture kits hash {REF}
##         : $capture_kit                    => The capture kit to add
##         : $user_supplied_parameter_switch => Has user supplied parameter {OPTIONAL}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $supported_capture_kit_href;
    my $capture_kit;
    my $user_supplied_parameter_switch;

    my $tmpl = {
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        supported_capture_kit_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$supported_capture_kit_href
        },
        capture_kit => { strict_type => 1, store => \$capture_kit },
        user_supplied_parameter_switch =>
          { strict_type => 1, store => \$user_supplied_parameter_switch },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    unless ( defined($user_supplied_parameter_switch) )
    {    #No detected supplied capture kit

        if ( defined( $supported_capture_kit_href->{$capture_kit} ) )
        {    #Supported capture kit alias

            return $supported_capture_kit_href->{$capture_kit};
        }
        else {    #Return unchanged capture_kit string

            return $capture_kit;
        }
    }
    if (   ( defined($user_supplied_parameter_switch) )
        && ( !$user_supplied_parameter_switch ) )
    {             #Only add if user supplied no info on parameter

        if ( defined( $supported_capture_kit_href->{$capture_kit} ) )
        {         #Supported capture kit alias

            return $supported_capture_kit_href->{$capture_kit};
        }
        else {    #Return unchanged capture_kit string

            return $capture_kit;
        }
    }
}

sub split_bam {

##split_bam

##Function : Split BAM file per contig and index new BAM. Creates the command line for xargs. Writes to sbatch FILEHANDLE and opens xargs FILEHANDLE
##Returns  : ""
##Arguments: $FILEHANDLE, $XARGSFILEHANDLE, $contigs_ref, $xargs_file_counter, $file_path, $program_info_path, $core_number, $infile, $temp_directory_ref
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $FILEHANDLE            => Sbatch filehandle to write to
##         : $XARGSFILEHANDLE       => XARGS filehandle to write to
##         : $xargs_file_counter    => The xargs file counter
##         : $memory_allocation     => Memory allocation for ja
##         : $contigs_ref           => The contigs to process
##         : $file_path             => File name - ususally sbatch
##         : $program_info_path     => The program info path
##         : $core_number           => The number of cores to use
##         : $infile                => The infile
##         : $temp_directory_ref    => Temporary directory

    my ($arg_href) = @_;

    ## Default(s)
    my $temp_directory_ref;
    my $first_command;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $contigs_ref;
    my $FILEHANDLE;
    my $XARGSFILEHANDLE;
    my $file_path;
    my $program_info_path;
    my $core_number;
    my $infile;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        contigs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$contigs_ref
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE, },
        XARGSFILEHANDLE =>
          { required => 1, defined => 1, store => \$XARGSFILEHANDLE },
        file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_path
        },
        program_info_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_info_path
        },
        core_number => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$core_number
        },
        infile =>
          { required => 1, defined => 1, strict_type => 1, store => \$infile, },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Alignment::Samtools qw(samtools_view);
    use MIP::Language::Java qw{java_core};
    use MIP::Recipes::Xargs qw{ xargs_command };

    my $xargs_file_path_prefix;

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Split by contig
    foreach my $contig (@$contigs_ref) {

        samtools_view(
            {
                infile_path =>
                  catfile( $$temp_directory_ref, $infile . q{.bam} ),
                outfile_path => catfile(
                    $$temp_directory_ref, $infile . "_" . $contig . q{.bam}
                ),
                stderrfile_path => $xargs_file_path_prefix . "."
                  . $contig
                  . ".stderr.txt",
                regions_ref              => $contig,
                FILEHANDLE               => $XARGSFILEHANDLE,
                auto_detect_input_format => 1,
                with_header              => 1,
                uncompressed_bam_output  => 1,
            }
        );
        print $XARGSFILEHANDLE "; ";    #Wait

        ## Writes java core commands to filehandle.
        java_core(
            {
                FILEHANDLE        => $XARGSFILEHANDLE,
                memory_allocation => "Xmx4g",
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                temp_directory => $$temp_directory_ref,
                java_jar       => catfile(
                    $active_parameter_href->{picardtools_path},
                    q{picard.jar}
                ),
            }
        );

        print $XARGSFILEHANDLE "BuildBamIndex ";
        print $XARGSFILEHANDLE "INPUT="
          . catfile( $$temp_directory_ref, $infile . "_" . $contig . q{.bam} )
          . " ";    #InFile
        say {$XARGSFILEHANDLE} "2> "
          . $xargs_file_path_prefix . "."
          . $contig
          . ".stderr.txt "
          ;         #Redirect xargs output to program specific stderr file
    }
    return $xargs_file_counter;
}

sub find_max_seq_length_for_sample_id {

##find_max_seq_length_for_sample_id

##Function : Finds the maximum sequence length of the reads for all sequencing file(s).
##Returns  : $max_sequence_length
##Arguments: $active_parameter_href, $sample_info_href, $infile_lane_prefix_href, $sample_id_ref
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $sample_id_ref              => Sample id {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href   = $arg_href->{active_parameter_href};
    my $sample_info_href        = $arg_href->{sample_info_href};
    my $infile_lane_prefix_href = $arg_href->{infile_lane_prefix_href};
    my $sample_id_ref           = $arg_href->{sample_id_ref};

    my $max_sequence_length = 0;

    foreach my $infile ( @{ $infile_lane_prefix_href->{$$sample_id_ref} } )
    {    #For all infiles per lane

        my $seq_length =
          $sample_info_href->{sample}{$$sample_id_ref}{file}{$infile}
          {sequence_length};

        if ( $seq_length > $max_sequence_length ) {

            $max_sequence_length = $seq_length;
        }
    }
    return $max_sequence_length;
}

sub collect_gene_panels {

##collect_gene_panels

##Function : Collect databases(s) from a database file and adds them to sample_info
##Returns  : ""
##Arguments: $sample_info_href, $family_id_ref, $program_name_ref, $aggregate_gene_panel_file, $aggregate_gene_panels_key
##         : $sample_info_href          => Info on samples and family hash {REF}
##         : $family_id_ref             => The family ID {REF}
##         : $program_name_ref          => Program name {REF}
##         : $aggregate_gene_panel_file => The database file
##         : $aggregate_gene_panels_key => The database key i.e. select or range

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $family_id_ref;
    my $program_name_ref;
    my $aggregate_gene_panel_file;
    my $aggregate_gene_panels_key;

    my $tmpl = {
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        family_id_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$family_id_ref,
        },
        program_name_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$program_name_ref
        },
        aggregate_gene_panels_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$aggregate_gene_panels_key
        },
        aggregate_gene_panel_file =>
          { strict_type => 1, store => \$aggregate_gene_panel_file },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( defined($aggregate_gene_panel_file) ) {

        ## Retrieve logger object
        my $log = Log::Log4perl->get_logger(q{MIP});

        my %gene_panel;    #Collect each gene panel features
        my %header = (
            gene_panel   => "gene_panel",
            version      => "version",
            updated_at   => "updated_at",
            display_name => "display_name",
        );

        my $sub_database_regexp =
q?perl -nae 'if ($_=~/^##gene_panel=/) {chomp($_);my @entries=split(/,/, $_); my $entry = join(",", $_); print $entry.":" } if($_=~/^#\w/) {last;}'?;
        my $ret = `$sub_database_regexp $aggregate_gene_panel_file`
          ;                #Collect header_lines(s) from select_file header
        my @header_lines = split( /:/, $ret )
          ;    #Split each gene panel meta data header line into array element

      LINE:
        foreach my $line (@header_lines) {

            my @features = split( /,/, $line )
              ;    #Split each memember database line into features

          ELEMENT:
            foreach my $feature_element (@features) {

              KEY_VALUE:
                foreach my $gene_panel_header_element ( keys %header )
                {    #Parse the features using defined header keys

                    if ( $feature_element =~ /$gene_panel_header_element=/ ) {

                        my @temps = split( "=", $feature_element );
                        $gene_panel{ $header{$gene_panel_header_element} } =
                          $temps[1];    #Value
                        last;
                    }
                }
            }

            if ( defined( $gene_panel{gene_panel} ) ) {

                my $gene_panel_name =
                  $gene_panel{gene_panel};    #Create unique gene panel ID

                ## Add new entries
                foreach my $feature ( keys %gene_panel ) {

                    $sample_info_href->{$$program_name_ref}
                      {$aggregate_gene_panels_key}{gene_panel}
                      {$gene_panel_name}{$feature} = $gene_panel{$feature};
                }
            }
            else {

                $log->warn( "Unable to write "
                      . $aggregate_gene_panels_key
                      . " aggregate gene panel(s) to qc_sample_info. Lacking ##gene_panel=<ID=[?] or version=[?] in aggregate gene panel(s) header."
                      . "\n" );
            }
            %gene_panel = ();    #Reset hash for next line
        }
        $sample_info_href->{$$program_name_ref}{$aggregate_gene_panels_key}
          {path} = $aggregate_gene_panel_file;
    }
}

sub add_most_complete_vcf {

##add_most_complete_vcf

##Function : Adds the most complete vcf file to sample_info
##Returns  : ""
##Arguments: $active_parameter_href, $sample_info_href, $path, $program_name, $family_id_ref, $vcfparser_outfile_counter
##         : $active_parameter_href     => Active parameters for this analysis hash {REF}
##         : $sample_info_href          => Info on samples and family hash {REF}
##         : $path                      => Path to file
##         : $program_name              => Program name
##         : $family_id_ref             => The family ID {REF}
##         : $vcfparser_outfile_counter => Number of outfile files from in vcfParser (select, range)
##         : $vcf_file_key              => Key for labelling most complete vcf

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $vcfparser_outfile_counter;
    my $vcf_file_key;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;
    my $path;
    my $program_name;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        path =>
          { required => 1, defined => 1, strict_type => 1, store => \$path },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        vcfparser_outfile_counter => {
            default     => 0,
            strict_type => 1,
            store       => \$vcfparser_outfile_counter
        },
        vcf_file_key => {
            default     => "vcf_file",
            allow       => [ "vcf_file", "sv_vcf_file" ],
            strict_type => 1,
            store       => \$vcf_file_key
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        if ( $vcfparser_outfile_counter == 1 ) {

            $sample_info_href->{$vcf_file_key}{clinical}{path} = $path;
        }
        else {

            $sample_info_href->{$vcf_file_key}{research}{path} = $path;
        }
    }
}

sub check_command_in_path {

##check_command_in_path

##Function : Checking commands in your path and executable
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href
##         : $parameter_href        => Parameter hash {REF}
##         : $active_parameter_href => Active parameters for this analysis hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Unix qw{check_binary_in_path};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    # Track program paths that have already been checked
    my %seen;

    foreach my $parameter_name ( keys %$active_parameter_href ) {

        if (   ( exists( $parameter_href->{$parameter_name}{type} ) )
            && ( $parameter_href->{$parameter_name}{type} eq "program" ) )
        {

            my $program_name_paths_ref =
              \@{ $parameter{$parameter_name}{program_name_path} };    #Alias

            if (   (@$program_name_paths_ref)
                && ( $active_parameter_href->{$parameter_name} > 0 ) )
            {    #Only check path(s) for active programs

                foreach my $program ( @{$program_name_paths_ref} ) {

                    unless ( $seen{$program} ) {

                        $seen{$program} = check_binary_in_path(
                            {
                                binary => $program,
                                log    => $log,
                            }
                        );
                    }
                }
            }
        }
    }
}

sub update_sample_info_hash {

##update_sample_info_hash

##Function : Update sample_info with information from pedigree
##Returns  : ""
##Arguments: $sample_info_href, $temp_href, $family_id_ref
##         : $temp_href        => Allowed parameters from pedigre file hash {REF}
##         : $sample_info_href => Info on samples and family hash {REF}
##         : $family_id_ref    => The family ID {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $temp_href;
    my $family_id_ref;

    my $tmpl = {
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        temp_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$temp_href
        },
        family_id_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$family_id_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    foreach my $sample_id ( keys %{ $sample_info_href->{sample} } ) {

        foreach my $key ( keys %{ $sample_info_href->{sample}{$sample_id} } ) {

            if ( exists( $temp_href->{sample}{$sample_id}{$key} ) )
            { #Previous run information, which should be updated using pedigree from current analysis

                $temp_href->{sample}{$sample_id}{$key} =
                  delete( $sample_info_href->{sample}{$sample_id}{$key} )
                  ;    #Required to update info
            }
            else {

                $temp_href->{sample}{$sample_id}{$key} =
                  $sample_info_href->{sample}{$sample_id}{$key};
            }
        }
    }
    %{$sample_info_href} = %$temp_href
      ; #Copy hash with updated keys from what was in sample_info (should be only pedigree %allowed_entries)
}

sub update_to_absolute_path {

##update_to_absolute_path

##Function : Change relative path to absolute path for certain parameter_names
##Returns  : ""
##Arguments: $parameter_href
##         : $parameter_href => Parameter hash {REF}

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $parameter_href;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Set::File qw(set_absolute_path);

    ## Adds dynamic aggregate information from definitions to parameter hash
    # Collect all path that should be made absolute
    set_dynamic_parameter(
        {
            parameter_href => \%parameter,
            aggregates_ref => [q{update_path:absolute_path}],
        }
    );

    foreach my $parameter_name (
        @{ $parameter_href->{dynamic_parameter}{absolute_path} } )
    {

        my $parameter_type =
          $parameter_href->{$parameter_name}{data_type};    #Alias

        if ( $parameter_type eq "ARRAY" ) {                 #Array reference

            foreach my $parameter_value (
                @{ $parameter_href->{$parameter_name}{value} } )
            {

                if ( $parameter_value ne "nocmd_input" ) {

                    my $element_separator_ref =
                      \$parameter_href->{$parameter_name}{element_separator}
                      ;    #Alias - Find what seperates array
                    my @seperate_elements =
                      split( $$element_separator_ref, $parameter_value )
                      ; #Split into seperate elements if written in 1 string on cmd
                    my @absolute_path_elements;

                    foreach my $element (@seperate_elements) {

                        ## Find aboslute path for supplied path or croaks and exists if path does not exists
                        push(
                            @absolute_path_elements,
                            set_absolute_path(
                                {
                                    path           => $element,
                                    parameter_name => $parameter_name,
                                }
                            )
                        );
                    }
                    $parameter_value = join( ",", @absolute_path_elements )
                      ;    #Replace original input with abolute path entries
                }
            }
        }
        elsif ( $parameter_type eq "HASH" ) {    #Hash reference

            my $parameter_value_ref =
              $parameter_href->{$parameter_name}{value};    #Alias

            foreach my $key ( keys %$parameter_value_ref )
            {    #Cannot use each since we are updating key

                if (   ( defined($parameter_value_ref) )
                    && ( $parameter_value_ref->{$key} ne "nocmd_input" ) )
                {

                    ## Find aboslute path for supplied path or croaks and exists if path does not exists
                    my $updated_key = set_absolute_path(
                        {
                            path           => $key,
                            parameter_name => $parameter_name,
                        }
                    );
                    $parameter_value_ref->{$updated_key} =
                      delete( $parameter_value_ref->{$key} );
                }
            }
        }
        elsif ( $parameter_type eq "SCALAR" ) {    #Scalar - not a reference

            my $parameter_value =
              $parameter_href->{$parameter_name}{value};    #Alias

            if (   ( defined($parameter_value) )
                && ( $parameter_value ne "nocmd_input" ) )
            {                                               #Scalar

                $parameter_value = set_absolute_path(
                    {
                        path           => $parameter_value,
                        parameter_name => $parameter_name,
                    }
                );
                $parameter_href->{$parameter_name}{value} = $parameter_value;
            }
        }
    }
}

sub add_to_sample_info {

##add_to_sample_info

##Function : Adds parameter info to sample_info
##Returns  : ""
##Arguments: $active_parameter_href, $sample_info_href, $file_info_href, $family_id_ref
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $sample_info_href      => Info on samples and family hash {REF}
##         : $file_info_href        => File info hash {REF}
##         : $family_id_ref         => The family_id_ref {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $human_genome_reference_ref;
    my $outdata_dir;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        human_genome_reference_ref => {
            default =>
              \$arg_href->{active_parameter_href}{human_genome_reference},
            strict_type => 1,
            store       => \$human_genome_reference_ref
        },
        outdata_dir => {
            default     => $arg_href->{active_parameter_href}{outdata_dir},
            strict_type => 1,
            store       => \$outdata_dir
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( exists( $active_parameter_href->{analysis_type} ) ) {

        $sample_info_href->{analysis_type} =
          $active_parameter_href->{analysis_type};
    }
    if ( exists( $active_parameter_href->{expected_coverage} ) ) {

        $sample_info_href->{expected_coverage} =
          $active_parameter_href->{expected_coverage};
    }
    if ( exists $active_parameter_href->{gatk_path} ) {

        my $gatk_version;
        if ( $active_parameter_href->{gatk_path} =~ /GenomeAnalysisTK-([^,]+)/ )
        {

            $gatk_version = $1;
        }
        else {    #Fall back on actually calling program

            my $jar_path = catfile( $active_parameter_href->{gatk_path},
                "GenomeAnalysisTK.jar" );
            $gatk_version = (`java -jar $jar_path --version 2>&1`);
            chomp $gatk_version;
        }
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => 'gatk',
                version          => $gatk_version,
            }
        );
    }
    if ( exists $active_parameter_href->{picardtools_path} )
    {    #To enable addition of version to sample_info

        my $picardtools_version;
        if ( $active_parameter_href->{picardtools_path} =~
            /picard-tools-([^,]+)/ )
        {

            $picardtools_version = $1;
        }
        else {    #Fall back on actually calling program

            my $jar_path = catfile( $active_parameter_href->{picardtools_path},
                q{picard.jar} );
            $picardtools_version =
              (`java -jar $jar_path CreateSequenceDictionary --version 2>&1`);
            chomp $picardtools_version;
        }

        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => 'picardtools',
                version          => $picardtools_version,
            }
        );
    }
    my @sambamba_programs =
      ( "pbwa_mem", "psambamba_depth", "markduplicates_sambamba_markdup" );
    foreach my $program (@sambamba_programs) {

        if (   ( exists $active_parameter_href->{$program} )
            && ( $active_parameter_href->{$program} == 1 ) )
        {

            if ( !$active_parameter_href->{dry_run_all} ) {

                my $regexp =
                  q?perl -nae 'if($_=~/sambamba\s(\S+)/) {print $1;last;}'?;
                my $sambamba_version = (`sambamba 2>&1 | $regexp`);
                chomp $sambamba_version;
                add_program_outfile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        program_name     => 'sambamba',
                        version          => $sambamba_version,
                    }
                );
                last;    #Only need to check once
            }
        }
    }
    if ( exists $active_parameter_href->{pcnvnator} )
    {                    #To enable addition of version to sample_info

        if (   ( $active_parameter_href->{pcnvnator} == 1 )
            && ( !$active_parameter_href->{dry_run_all} ) )
        {

            my $regexp =
              q?perl -nae 'if($_=~/CNVnator\s+(\S+)/) {print $1;last;}'?;
            my $cnvnator_version = (`cnvnator 2>&1 | $regexp`);
            chomp $cnvnator_version;
            add_program_outfile_to_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    program_name     => 'cnvnator',
                    version          => $cnvnator_version,
                }
            );
        }
    }
    if ( defined($$human_genome_reference_ref) )
    {    #To enable addition of version to sample_info

        $sample_info_href->{human_genome_build}{path} =
          $$human_genome_reference_ref;
        $sample_info_href->{human_genome_build}{source} =
          $file_info_href->{human_genome_reference_source};
        $sample_info_href->{human_genome_build}{version} =
          $file_info_href->{human_genome_reference_version};
    }
    if ( exists( $active_parameter_href->{pedigree_file} ) ) {

        $sample_info_href->{pedigree_file}{path} =
          $active_parameter_href->{pedigree_file}
          ;    #Add pedigree_file to sample_info
        $sample_info_href->{pedigree_file_analysis}{path} =
          catfile( $outdata_dir, $$family_id_ref, "qc_pedigree.yaml" )
          ;    #Add pedigree_file info used in this analysis to sample_info_file
    }
    if ( exists( $active_parameter_href->{log_file} ) ) {

        my $path = dirname( dirname( $active_parameter_href->{log_file} ) );
        $sample_info_href->{log_file_dir} =
          $path;    #Add log_file_dir to SampleInfoFile
        $sample_info_href->{last_log_file_path} =
          $active_parameter_href->{log_file};
    }
}

sub check_config_vs_definition_file {

##check_config_vs_definition_file

##Function : Compare keys from config and definitions file
##Returns  : ""
##Arguments: $reference_href, $comparison_href
##         : $reference_href  => Reference hash {REF}
##         : $comparison_href => Hash to be compared to reference {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $reference_href;
    my $comparison_href;

    my $tmpl = {
        reference_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$reference_href
        },
        comparison_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$comparison_href
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @allowed_unique_keys =
      ( "vcfparser_outfile_count", $reference_href->{family_id} );
    my @unique;

    foreach my $key ( keys %$reference_href ) {

        unless ( exists( $comparison_href->{$key} ) ) {

            push( @unique, $key );
        }
    }
    foreach my $element (@unique) {

        if ( !( any { $_ eq $element } @allowed_unique_keys ) )
        { #Do not print if allowed_unique_keys that have been created dynamically from previous runs

            warn(   "Found illegal key: "
                  . $element
                  . " in config file that is not defined in definitions.yaml\n"
            );
            exit 1;
        }
    }
}

sub check_vep_directories {

##check_vep_directories

##Function : Compare VEP directory and VEP chache versions
##Returns  : ""
##Arguments: $vep_directory_path_ref, $vep_directory_cache_ref
##         : $vep_directory_path_ref  => VEP directory path {REF}
##         : $vep_directory_cache_ref => VEP cache directory path {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vep_directory_path_ref;
    my $vep_directory_cache_ref;

    my $tmpl = {
        vep_directory_path_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$vep_directory_path_ref
        },
        vep_directory_cache_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$vep_directory_cache_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    if ( $$vep_directory_path_ref =~ /ensembl-tools-release-(\d+)/ ) {

        my $vep_directory_path_version = $1;

        unless ( $$vep_directory_cache_ref =~
            /ensembl-tools-release-$vep_directory_path_version/ )
        {

            print $log->fatal(
                "Differing versions between '-vep_directory_path': "
                  . $$vep_directory_path_ref
                  . " and '-vep_directory_cache': "
                  . $$vep_directory_cache_ref,
                "\n"
            );
            exit 1;
        }
    }

}

sub remove_redundant_files {

##remove_redundant_files

##Function : Removes intermediate files from the MIP analysis depending on set MIP parameters
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $lane_href, $infile_lane_prefix_href, $reduce_io_ref, $sample_id, $insample_directory, $FILEHANDLE, family_id_ref, $outaligner_dir_ref, $call_type
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => File info hash {REF}
##         : $lane_href                  => The lane info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $reduce_io_ref              => Reduce IO - modulates processBlocks {REF}
##         : $sample_id                  => Sample id
##         : $insample_directory         => The directory for in sample files to be removed
##         : $FILEHANDLE                 => Filehandle to write to
##         : $family_id_ref              => Family id {REF}
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $call_type                  => Variant call type

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $outaligner_dir_ref;
    my $call_type;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $lane_href;
    my $infile_lane_prefix_href;
    my $reduce_io_ref;
    my $sample_id;
    my $insample_directory;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        lane_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$lane_href
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        reduce_io_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$reduce_io_ref
        },
        sample_id => { strict_type => 1, store => \$sample_id },
        insample_directory =>
          { strict_type => 1, store => \$insample_directory },
        FILEHANDLE    => { store => \$FILEHANDLE, },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $vcfparser_contigs_ref =
      \@{ $file_info_href->{contigs_size_ordered} };    #Set default

    ## Last modules in each processing block that should have the output data deleted
    my $last_module_bamcalibrationblock    = "pgatk_haplotypecaller";
    my $last_module_variantannotationblock = "psnpeff";

    foreach my $program ( @{ $parameter_href->{dynamic_parameter}{program} } ) {

        if ( $active_parameter_href->{$program} > 0 ) {

            if (
                (
                    defined(
                        $parameter_href->{$program}{remove_redundant_file}
                    )
                )
                && ( $parameter_href->{$program}{remove_redundant_file} eq
                    "yes" )
              )
            {

                if ( defined($sample_id) ) {

                    my $indirectory =
                      $parameter_href->{$program}{$sample_id}{indirectory};
                    my $outfile_tag =
                      $file_info_href->{$sample_id}{$program}{file_tag};

                    ## Single files
                    if ( $parameter_href->{$program}
                        {remove_redundant_file_setting} eq "single" )
                    {    #Infiles for prior to potential merge

                      INFILE:
                        foreach my $infile (
                            @{ $infile_lane_prefix_href->{$sample_id} } )
                        {

                          FILE_ENDINGS:
                            foreach my $file_ending (
                                @{ $parameter_href->{$program}{file_endings} } )
                            {

                                my $file_path;

                                if ( defined($outfile_tag) ) {

                                    $file_path = catfile( $indirectory,
                                        $infile . $outfile_tag . $file_ending );
                                }
                                else {

                                    $file_path = catfile( $indirectory,
                                        $infile . $file_ending );
                                }

                                my $most_complete_ref;

                                if ( $file_ending =~ /vcf|bam/ ) {

                                    ## Detect which most_complete_path to use depending on file_ending
                                    $most_complete_ref =
                                      detect_most_complete_file(
                                        {
                                            sample_info_href =>
                                              $sample_info_href,
                                            file_ending_ref => \$file_ending,
                                            sample_id_ref   => \$sample_id,
                                        }
                                      );
                                }
                                ## Checks if the file is recorded as the "most_complete_bam|vcf". If false writes removal of file(s) to supplied filehandle
                                check_most_complete_and_remove_file(
                                    {
                                        FILEHANDLE        => $FILEHANDLE,
                                        most_complete_ref => $most_complete_ref,
                                        file_path_ref     => \$file_path,
                                        file_ending       => $file_ending,
                                    }
                                );
                            }
                        }
                    }
                    ## Merged files
                    if ( $parameter_href->{$program}
                        {remove_redundant_file_setting} eq "merged" )
                    {    #Merge infiles

                        ## Add merged infile name after merging all BAM files per sample_id
                        my $infile =
                          $file_info_href->{$sample_id}{merged_infile};   #Alias

                        if (   ( !$$reduce_io_ref )
                            || ( $program eq $last_module_bamcalibrationblock )
                          )
                        { #Delete intermediate files or last module in processBlock

                          FILE_ENDINGS:
                            foreach my $file_ending (
                                @{ $parameter_href->{$program}{file_endings} } )
                            {

                              CONTIGS:
                                foreach my $contig (@$vcfparser_contigs_ref) {

                                    my $file_path = catfile( $indirectory,
                                            $infile
                                          . $outfile_tag . "_"
                                          . $contig
                                          . $file_ending );

                                    ## Detect which most_complete_path to use depending on file_ending
                                    my $most_complete_ref =
                                      detect_most_complete_file(
                                        {
                                            sample_info_href =>
                                              $sample_info_href,
                                            file_ending_ref => \$file_ending,
                                            sample_id_ref   => \$sample_id,
                                        }
                                      );

                                    ## Checks if the file is recorded as the "most_complete_bam|vcf". If false writes removal of file(s) to supplied filehandle
                                    check_most_complete_and_remove_file(
                                        {
                                            FILEHANDLE => $FILEHANDLE,
                                            most_complete_ref =>
                                              $most_complete_ref,
                                            file_path_ref => \$file_path,
                                            file_ending   => $file_ending,
                                        }
                                    );
                                }
                            }
                        }
                    }
                }
                else
                {    #Otherwise these files would be removed for every sample_id

                    my $indirectory = $parameter_href->{$program}{indirectory};
                    my $outfile_tag =
                      $file_info_href->{$$family_id_ref}{$program}{file_tag};

                    ## Family files
                    if ( $parameter_href->{$program}
                        {remove_redundant_file_setting} eq "family" )
                    {

                      FILE_ENDINGS:
                        foreach my $file_ending (
                            @{ $parameter_href->{$program}{file_endings} } )
                        {

                            my $file_path = catfile( $indirectory,
                                    $$family_id_ref
                                  . $outfile_tag
                                  . $call_type . "*"
                                  . $file_ending );

                            ## Detect which most_complete_path to use depending on file_ending
                            my $most_complete_ref = detect_most_complete_file(
                                {
                                    sample_info_href => $sample_info_href,
                                    file_ending_ref  => \$file_ending,
                                }
                            );

                            ## Checks if the file is recorded as the "most_complete_bam|vcf". If false writes removal of file(s) to supplied filehandle
                            check_most_complete_and_remove_file(
                                {
                                    FILEHANDLE        => $FILEHANDLE,
                                    most_complete_ref => $most_complete_ref,
                                    file_path_ref     => \$file_path,
                                    file_ending       => $file_ending,
                                }
                            );
                        }
                    }
                    elsif ( $parameter_href->{$program}
                        {remove_redundant_file_setting} eq
                        "variant_annotation" )
                    {

                        if ( ( !$$reduce_io_ref )
                            || ( $program eq
                                $last_module_variantannotationblock ) )
                        { #Delete intermediate files or last module in processBlock

                            $outfile_tag =
                              $file_info_href->{$$family_id_ref}{$program}
                              {file_tag};

                            foreach my $file_ending (
                                @{ $parameter_href->{$program}{file_endings} } )
                            {

                                my $file_path = catfile( $indirectory,
                                        $$family_id_ref
                                      . $outfile_tag
                                      . $call_type . "*"
                                      . $file_ending );

                                ## Detect which most_complete_path to use depending on file_ending
                                my $most_complete_ref =
                                  detect_most_complete_file(
                                    {
                                        sample_info_href => $sample_info_href,
                                        file_ending_ref  => \$file_ending,
                                    }
                                  );

                                ## Checks if the file is recorded as the "most_complete_bam|vcf". If false writes removal of file(s) to supplied filehandle
                                check_most_complete_and_remove_file(
                                    {
                                        FILEHANDLE        => $FILEHANDLE,
                                        most_complete_ref => $most_complete_ref,
                                        file_path_ref     => \$file_path,
                                        file_ending       => $file_ending,
                                    }
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}

sub detect_most_complete_file {

##detect_most_complete_file

##Function : Detect which most_complete_path to use depending on file_ending
##Returns  : ""
##Arguments: $active_parameter_href, $file_ending_ref, $sample_id_ref
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $file_ending_ref       => File ending (.file_ending){REF}
##         : $sample_id_ref         => Sample ID {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $file_ending_ref;
    my $sample_id_ref;

    my $tmpl = {
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_ending_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$file_ending_ref
        },
        sample_id_ref => { store => \$sample_id_ref },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set mostcompletePaths
    my $most_complete_bam_ref;
    my $most_complete_vcf_ref = \$sample_info_href->{vcf_file}{ready_vcf}{path};

    if ( defined($$sample_id_ref) ) {

        $most_complete_bam_ref =
          \$sample_info_href->{sample}{$$sample_id_ref}{most_complete_bam}
          {path};
    }

    ## Decide which mostcompletePaths to use
    if ( $$file_ending_ref eq q{.bam} ) {

        return $most_complete_bam_ref;
    }
    if ( $$file_ending_ref eq ".vcf" ) {

        return $most_complete_vcf_ref;
    }
}

sub remove_array_element {

##remove_array_element

##Function : Removes contigs from supplied contigs_ref
##Returns  : ""
##Arguments: $contigs_ref, $remove_contigs_ref
##         : $contigs_ref        => The select file contigs {REF}
##         : $remove_contigs_ref => Remove this contig

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contigs_ref;
    my $remove_contigs_ref;

    my $tmpl = {
        contigs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$contigs_ref
        },
        remove_contigs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$remove_contigs_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    for ( my $index = 0 ; $index < scalar(@$contigs_ref) ; $index++ ) {

        foreach my $remove_contig (@$remove_contigs_ref) {

            if (   ( $contigs_ref->[$index] eq $remove_contig )
                || ( $contigs_ref->[$index] eq "chr" . $remove_contig ) )
            {

                splice( @$contigs_ref, $index, 1 );  #Remove $element from array
            }
        }
    }
}

sub detect_founders {

##detect_founders

##Function : Detect number of founders (i.e. parents ) based on pedigree file
##Returns  : ""|1
##Arguments: $active_parameter_href,
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $sample_info_href      => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @founders;

    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        my $father_info =
          $sample_info_href->{sample}{$sample_id}{father};    #Alias
        my $mother_info =
          $sample_info_href->{sample}{$sample_id}{mother};    #Alias

        if ( ( defined($father_info) ) && ( $father_info ne 0 ) ) {    #Child

            if (
                any { $_ eq $father_info }
                @{ $active_parameter_href->{sample_ids} }
              )
            {    #If element is part of array

                push( @founders, $father_info );
            }
        }
        if ( ( defined($mother_info) ) && ( $mother_info ne 0 ) ) {    #Child

            if (
                any { $_ eq $mother_info }
                @{ $active_parameter_href->{sample_ids} }
              )
            {    #If element is part of array

                push( @founders, $mother_info );
            }
        }
    }
    return scalar(@founders);
}

sub detect_trio {

##detect_trio

##Function : Detect family constellation based on pedigree file
##Returns  : ""|1
##Arguments: $active_parameter_href, $sample_info_href
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $sample_info_href      => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my %trio;

    if ( scalar( @{ $active_parameter_href->{sample_ids} } ) eq 1 ) {

        $log->info(
            "Found single sample: " . $active_parameter_href->{sample_ids}[0],
            "\n" );
        return;
    }
    elsif ( scalar( @{ $active_parameter_href->{sample_ids} } ) eq 3 ) {

        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $father_info =
              $sample_info_href->{sample}{$sample_id}{father};    #Alias
            my $mother_info =
              $sample_info_href->{sample}{$sample_id}{mother};    #Alias

            if ( ( $father_info ne 0 ) && ( $mother_info ne 0 ) ) {    #Child

                $trio{child} = $sample_id;

                if (
                    any { $_ eq $father_info }
                    @{ $active_parameter_href->{sample_ids} }
                  )
                {    #If element is part of array

                    $trio{father} = $father_info;
                }
                if (
                    any { $_ eq $mother_info }
                    @{ $active_parameter_href->{sample_ids} }
                  )
                {    #If element is part of array

                    $trio{mother} = $mother_info;
                }
            }
        }
        if ( scalar( keys %trio ) == 3 ) {

            $log->info(
                "Found trio: Child = "
                  . $trio{child}
                  . ", Father = "
                  . $trio{father}
                  . ", Mother = "
                  . $trio{mother},
                "\n"
            );
            return 1;
        }
    }
}

sub check_string {

##check_string

##Function : Detect "regexp" in string
##Returns  : ""|1
##Arguments: $string
##         : $string => String to be searched
##         : $regexp => regexp to use on string

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $string;
    my $regexp;

    my $tmpl = {
        string =>
          { required => 1, defined => 1, strict_type => 1, store => \$string },
        regexp =>
          { required => 1, defined => 1, strict_type => 1, store => \$regexp },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( $string =~ /$regexp/ ) {

        return 1;
    }
}

sub check_prioritize_variant_callers {

## Function : Check that all active variant callers have a prioritization order and that the prioritization elements match a supported variant caller.
## Returns  :
## Arguments: $parameter_href        => Parameter hash {REF}
##          : $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $variant_callers_ref   => Variant callers to check {REF}
##          : $parameter_names_ref   => Parameter name list {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $variant_callers_ref;
    my $parameter_names_ref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
            ,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
            ,
        },
        parameter_names_ref => {
            required => 1,
            defined  => 1,
            default  => [],
            store    => \$parameter_names_ref,
        },
        variant_callers_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$variant_callers_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my @priority_calls =
      split( ",", $active_parameter_href->{$$parameter_names_ref} );

    ## No matching variant caller
    my @variant_caller_aliases;

    ## Check that all active variant callers have a priority order
  CALLER:
    foreach my $variant_caller ( @{$variant_callers_ref} ) {

        my $variant_caller_alias =
          $parameter_href->{$variant_caller}{outdir_name};
        push @variant_caller_aliases, $variant_caller_alias;

        ## Only active programs
        if ( $active_parameter_href->{$variant_caller} > 0 ) {

            ## If element is not part of string
            if ( !( any { $_ eq $variant_caller_alias } @priority_calls ) ) {

                $log->fatal( $$parameter_names_ref
                      . q{ does not contain active variant caller: '}
                      . $variant_caller_alias
                      . q{'} );
                exit 1;
            }
        }
        ## Only NOT active programs
        if ( $active_parameter_href->{$variant_caller} == 0 ) {

            ## If element is part of string
            if ( ( any { $_ eq $variant_caller_alias } @priority_calls ) ) {

                $log->fatal( $$parameter_names_ref
                      . q{ contains deactivated variant caller: '}
                      . $variant_caller_alias
                      . q{'} );
                exit 1;
            }
        }
    }

    ## Check that prioritize string contains valid variant call names
    foreach my $prioritize_call (@priority_calls) {

        if ( !( any { $_ eq $prioritize_call } @variant_caller_aliases ) )
        {    #If element is not part of string

            $log->fatal( $$parameter_names_ref . ": '"
                  . $prioritize_call
                  . "' does not match any supported variant caller: '"
                  . join( ",", @variant_caller_aliases )
                  . "'" );
            exit 1;
        }
    }
}

sub view_vcf {

##view_vcf

##Function : Reformat variant calling file and index.
##Returns  : ""
##Arguments: $infile_path, $FILEHANDLE, $outfile_path_prefix, $output_type, $index, $index_type
##         : $infile_path            => Path to infile to compress and index
##         : $FILEHANDLE             => SBATCH script FILEHANDLE to print to
##         : $outfile_path_prefix => Out file path no file_ending {Optional}
##         : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##         : $index                  => Generate index of reformated file
##         : $index_type             => Type of index

    my ($arg_href) = @_;

    ## Default(s)
    my $output_type;
    my $index;
    my $index_type;

    ## Flatten argument(s)
    my $infile_path;
    my $FILEHANDLE;
    my $outfile_path_prefix;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE, },
        outfile_path_prefix =>
          { strict_type => 1, store => \$outfile_path_prefix },
        output_type => {
            default     => q{v},
            allow       => [qw{ b u z v }],
            strict_type => 1,
            store       => \$output_type
        },
        index => {
            default     => 1,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$index
        },
        index_type => {
            default     => q{csi},
            allow       => [ undef, qw{ csi tbi } ],
            strict_type => 1,
            store       => \$index_type
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_view bcftools_index };

    my $outfile_path;
    my %output_type_ending = (
        b => q{.bcf},
        u => q{.bcf},
        z => q{.vcf.gz},
        v => q{.vcf},
    );

    if ( defined $outfile_path_prefix ) {

        $outfile_path =
          $outfile_path_prefix . $output_type_ending{$output_type};
    }

    say {$FILEHANDLE} q{## Reformat variant calling file};
    bcftools_view(
        {
            infile_path  => $infile_path,
            outfile_path => $outfile_path,
            output_type  => $output_type,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

    if ($index) {

        say {$FILEHANDLE} q{## Index};
        bcftools_index(
            {
                infile_path => $outfile_path,
                output_type => $index_type,
                FILEHANDLE  => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} "\n";
    }
    return;
}

sub check_aligner {

##check_aligner

##Function : Check that the correct number of aligners is used in MIP and sets the outaligner_dir flag accordingly.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $broadcasts_ref, $outaligner_dir_ref
##         : $parameter_href        => Parameter hash {REF}
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $broadcasts_ref        => Holds the parameters info for broadcasting later {REF}
##         : $outaligner_dir_ref    => Outaligner_dir used in the analysis {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $outaligner_dir_ref;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $broadcasts_ref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        broadcasts_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$broadcasts_ref
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my %aligner;

    foreach my $aligner ( @{ $parameter_href->{dynamic_parameter}{aligners} } )
    {

        if ( $active_parameter_href->{$aligner} > 0 ) {    #Active aligner

            $aligner{total_active_aligner_count}++;
            push( @{ $aligner{active_aligners} }, $aligner );
            $parameter_href->{active_aligner} =
              $aligner;    #Save the active aligner for downstream use

            if ( $$outaligner_dir_ref eq "not_set_yet" ) {

                $$outaligner_dir_ref = $parameter_href->{$aligner}{outdir_name}
                  ;    #Set outaligner_dir parameter depending on active aligner

                my $info = "Set outaligner_dir to: " . $$outaligner_dir_ref;
                push( @$broadcasts_ref, $info );    #Add info to broadcasts
            }
        }
    }
    if (   ( exists( $aligner{total_active_aligner_count} ) )
        && ( $aligner{total_active_aligner_count} > 1 ) )
    {

        $log->fatal( $USAGE, "\n" );
        $log->fatal(
            "You have activate more than 1 aligner: "
              . join( ", ", @{ $aligner{active_aligners} } )
              . ". MIP currently only supports 1 aligner per analysis.",
            "\n"
        );
        exit 1;
    }
}

sub collect_read_length {

##collect_read_length

##Function : Collect read length from an infile
##Returns  : "readLength"
##Arguments: $directory, $read_file, $file
##         : $directory => Directory of file
##         : $read_file => Command used to read file
##         : $file      => File to parse

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $directory;
    my $read_file_command;
    my $file;

    my $tmpl = {
        directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$directory
        },
        read_file_command => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$read_file_command
        },
        file =>
          { required => 1, defined => 1, strict_type => 1, store => \$file },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $seq_length_regexp =
q?perl -ne 'if ($_!~/@/) {chomp($_);my $seq_length = length($_);print $seq_length;last;}' ?
      ;    #Prints sequence length and exits

    my $pwd = cwd();      #Save current direcory
    chdir($directory);    #Move to sample_id infile directory

    my $ret =
      `$read_file_command $file | $seq_length_regexp;`; #Collect sequence length
    return $ret;
}

sub check_program_mode {

##check_program_mode

##Function : Check correct value for program mode in MIP.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href
##         : $parameter_href        => Parameter hash {REF}
##         : $active_parameter_href => Active parameters for this analysis hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my @allowed_values = ( 0, 1, 2 );

  PROGRAMS:
    foreach my $program ( @{ $parameter_href->{dynamic_parameter}{program} } ) {

        if (
            !(
                any { $_ eq $active_parameter_href->{$program} }
                @allowed_values
            )
          )
        {    #If element is not part of array

            $log->fatal( q{'}
                  . $active_parameter_href->{$program}
                  . q{' Is not an allowed mode for program '--}
                  . $program
                  . q{'. Set to: }
                  . join( "|", @allowed_values ) );
            exit 1;
        }
    }
}

sub update_exome_target_bed {

##update_exome_target_bed

##Function : Update exome_target_bed files with human_genome_reference_source_ref and human_genome_reference_version_ref
##Returns  : ""
##Arguments: $exome_target_bed_file_href, $human_genome_reference_source_ref, human_genome_reference_version_ref
##         : $exome_target_bed_file_href        => ExomeTargetBedTestFile hash {REF}
##         : human_genome_reference_source_ref  => The human genome reference source {REF}
##         : human_genome_reference_version_ref => The human genome reference version {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $exome_target_bed_file_href;
    my $human_genome_reference_source_ref;
    my $human_genome_reference_version_ref;

    my $tmpl = {
        exome_target_bed_file_href =>
          { required => 1, store => \$exome_target_bed_file_href },
        human_genome_reference_source_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$human_genome_reference_source_ref
        },
        human_genome_reference_version_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$human_genome_reference_version_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    foreach my $exome_target_bed_file ( keys %$exome_target_bed_file_href ) {

        my $original_file_name = $exome_target_bed_file;

        if (
            (
                $exome_target_bed_file =~
                s/genome_reference_source/$$human_genome_reference_source_ref/
            )
            && ( $exome_target_bed_file =~
                s/_version/$$human_genome_reference_version_ref/ )
          )
        {    #Replace with actual version

            $exome_target_bed_file_href->{$exome_target_bed_file} =
              delete( $exome_target_bed_file_href->{$original_file_name} )
              ; #The delete operator returns the value being deleted i.e. updating hash key while preserving original info
        }
    }
}

sub check_sample_id_in_parameter_path {

##check_sample_id_in_parameter_path

##Function : Check sample_id provided in hash path parameter is included in the analysis and only represented once
##Returns  : ""
##Tags     : check, sampleids, hash
##Arguments: $active_parameter_href, $sample_ids_ref, $parameter_name
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $sample_ids_ref        => Array to loop in for parameter {REF}
##         : $parameter_names_ref   => Parameter name list {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_ids_ref;
    my $parameter_names_ref;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_ids_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$sample_ids_ref
        },
        parameter_names_ref => {
            required => 1,
            defined  => 1,
            default  => [],
            store    => \$parameter_names_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    foreach my $parameter_name (@$parameter_names_ref)
    {    #Lopp through all hash parameters supplied

        my %seen;    #Hash to test duplicate sample_ids later

        foreach my $key ( keys %{ $active_parameter_href->{$parameter_name} } )
        {

            my @parameter_samples =
              split( ",", $active_parameter_href->{$parameter_name}{$key} );

            foreach my $sample_id (@parameter_samples) {

                $seen{$sample_id}++
                  ;    #Increment instance to check duplicates later

                if ( $seen{$sample_id} > 1 ) {    #Check sample_id are unique

                    $log->fatal(
                        "Sample_id: "
                          . $sample_id
                          . " is not uniqe in '-"
                          . $parameter_name . " '"
                          . $key . "="
                          . join( ",", @parameter_samples ),
                        "\n"
                    );
                    exit 1;
                }
            }
        }
        foreach my $sample_id (@$sample_ids_ref) {

            if ( !( any { $_ eq $sample_id } ( keys %seen ) ) )
            {    #If sample_id is not present in parameter_name hash

                $log->fatal(
                    "Could not detect "
                      . $sample_id
                      . " for '--"
                      . $parameter_name
                      . "'. Provided sample_ids are: "
                      . join( ", ", ( keys %seen ) ),
                    "\n"
                );
                exit 1;
            }
        }
    }
}

sub check_sample_id_in_parameter {

##check_sample_id_in_parameter

##Function : Check sample_id provided in hash parameter is included in the analysis and only represented once
##Returns  : ""
##Tags     : check, sampleids, hash
##Arguments: $active_parameter_href, $sample_ids_ref, $parameter_name
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $sample_ids_ref        => Array to loop in for parameter {REF}
##         : $parameter_names_ref   => Parameter name list {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_ids_ref;
    my $parameter_names_ref;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_ids_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$sample_ids_ref
        },
        parameter_names_ref => {
            required => 1,
            defined  => 1,
            default  => [],
            store    => \$parameter_names_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    foreach my $parameter_name (@$parameter_names_ref)
    {    #Lopp through all hash parameters supplied

        if ( defined( $active_parameter_href->{$parameter_name} ) ) {

            foreach my $sample_id (@$sample_ids_ref) {

                ## Check that a value exists
                if (
                    !defined(
                        $active_parameter_href->{$parameter_name}{$sample_id}
                    )
                  )
                {

                    $log->fatal(
                        "Could not find value for "
                          . $sample_id
                          . " for parameter '--"
                          . $parameter_name . "'",
                        "\n"
                    );
                    exit 1;
                }

                ## If sample_id is not present in parameter_name hash
                if (
                    !(
                        any { $_ eq $sample_id }
                        ( keys %{ $active_parameter_href->{$parameter_name} } )
                    )
                  )
                {

                    $log->fatal(
                        "Could not detect "
                          . $sample_id
                          . " for parameter '--"
                          . $parameter_name
                          . "'. Provided sample_ids for parameter are: "
                          . join(
                            ", ",
                            (
                                keys
                                  %{ $active_parameter_href->{$parameter_name} }
                            )
                          ),
                        "\n"
                    );
                    exit 1;
                }
            }
        }
    }
}

sub replace_iupac {

##replace_iupac

##Function : Replace the IUPAC code in alternative allels with N for input stream and writes to stream.
##Returns  : ""
##Arguments: $FILEHANDLE, $stderr_path, $xargs
##         : $FILEHANDLE  => Sbatch filehandle to write to
##         : $stderr_path => Stderr path to errors write to
##         : $xargs       => Write on xargs format

    my ($arg_href) = @_;

    ## Default(s)
    my $xargs;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $stderr_path;

    my $tmpl = {
        FILEHANDLE => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$FILEHANDLE
        },
        stderr_path => { strict_type => 1, store => \$stderr_path },
        xargs       => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$xargs
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    print {$FILEHANDLE} "| ";           #Pipe
    print {$FILEHANDLE} "perl -nae ";

    if ($xargs) {                       #Add escape char

        print {$FILEHANDLE}
q?\'if($_=~/^#/) {print $_;} else { @F[4] =~ s/W|K|Y|R|S|M/N/g; print join(\"\\\t\", @F), \"\\\n\"; }\' ?
          ; #substitute IUPAC code with N to not break vcf specifications (GRCh38)
    }
    else {

        print {$FILEHANDLE}
q?'if($_=~/^#/) {print $_;} else { @F[4] =~ s/W|K|Y|R|S|M/N/g; print join("\t", @F), "\n"; }' ?
          ; #substitute IUPAC code with N to not break vcf specifications (GRCh38)
    }
    if ($stderr_path) {

        print {$FILEHANDLE} "2>> "
          . $stderr_path
          . " ";    #Redirect output to program specific stderr file
    }
}

sub get_matching_values_key {

##get_matching_values_key

##Function : Return the key if the hash value and query match
##Returns  : "key pointing to matched value"
##Arguments: $active_parameter_href, $query_value_ref, $parameter_name
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $query_value_ref       => The value to query in the hash {REF}
##         : $parameter_name        => MIP parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $query_value_ref;
    my $parameter_name;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        query_value_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$query_value_ref
        },
        parameter_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %reversed = reverse %{ $active_parameter_href->{$parameter_name} }
      ;    #Values are now keys and vice versa

    if ( exists $reversed{$$query_value_ref} ) {

        return $reversed{$$query_value_ref};
    }
}

sub get_user_supplied_info {

##get_user_supplied_info

##Function : Detect if user supplied info on parameters otherwise collected from pedigree
##Returns  : "user_supply_switchHash where 1=user input and 0=no user input"
##Arguments: $parameter_href, active_parameter_href
##         : $parameter_href        => Holds all parameters {REF}
##         : $active_parameter_href => Active parameters for this analysis hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Define what should be checked
    my %user_supply_switch = (
        sample_ids        => 0,
        exome_target_bed  => 0,
        analysis_type     => 0,
        expected_coverage => 0,
    );

    ## Detect user supplied info
    foreach my $parameter ( keys %user_supply_switch ) {

        if ( ref( $parameter_href->{$parameter}{value} ) eq "HASH" ) {

            $user_supply_switch{$parameter} = check_user_supplied_info(
                {
                    active_parameter_href => $active_parameter_href,
                    data_ref => \%{ $parameter_href->{$parameter}{value} },
                    parameter_name => $parameter,
                }
            );
        }
        if ( ref( $parameter_href->{$parameter}{value} ) eq "ARRAY" ) {

            $user_supply_switch{$parameter} = check_user_supplied_info(
                {
                    active_parameter_href => $active_parameter_href,
                    data_ref => \@{ $parameter_href->{$parameter}{value} },
                    parameter_name => $parameter,
                }
            );
        }
        else {

            $user_supply_switch{$parameter} = check_user_supplied_info(
                {
                    active_parameter_href => $active_parameter_href,
                    data_ref       => $parameter_href->{$parameter}{value},
                    parameter_name => $parameter,
                }
            );
        }
    }
    return %user_supply_switch;
}

sub get_pedigree_sample_info {

##get_pedigree_sample_info

##Function : Reformat pedigree keys to plink format and collect sample info to various hashes
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $exom_target_bed_test_file_tracker_href, $pedigree_sample_href, $reformatHashRef, $user_supply_switch_href, $sample_id
##         : $parameter_href                         => Parameter hash {REF}
##         : $active_parameter_href                  => Active parameters for this analysis hash {REF}
##         : $sample_info_href                       => Info on samples and family hash {REF}
##         : $file_info_href                         => The associated reference file endings {REF}
##         : $exom_target_bed_test_file_tracker_href => Collect which sample_ids have used a certain capture_kit
##         : $pedigree_sample_href                   => YAML sample info hash {REF}
##         : $user_supply_switch_href                => The user supplied info switch {REF}
##         : $sample_id                              => Sample ID

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $exom_target_bed_test_file_tracker_href;
    my $pedigree_sample_href;
    my $user_supply_switch_href;
    my $sample_id;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        exom_target_bed_test_file_tracker_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$exom_target_bed_test_file_tracker_href
        },
        pedigree_sample_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$pedigree_sample_href
        },
        user_supply_switch_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$user_supply_switch_href
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Add input to sample_info hash for at sample level
    foreach my $key ( keys %$pedigree_sample_href ) {

        $sample_info_href->{sample}{$sample_id}{$key} =
          $pedigree_sample_href->{$key};

        ## Add sex to dynamic parameters
        if ( $key eq "sex" ) {

            push(
                @{
                    $parameter_href->{dynamic_parameter}
                      { $pedigree_sample_href->{$key} }
                },
                $sample_id
            );

            ## Reformat to plink format
            if ( $pedigree_sample_href->{$key} eq "male" ) {

                $parameter_href->{dynamic_parameter}{$sample_id}{plink_sex} = 1;
            }
            elsif ( $pedigree_sample_href->{$key} eq "female" ) {

                $parameter_href->{dynamic_parameter}{$sample_id}{plink_sex} = 2;
            }
            else {

                $parameter_href->{dynamic_parameter}{$sample_id}{plink_sex} =
                  "other";
            }
        }

        ## Add phenotype to dynamic parameters
        if ( $key eq "phenotype" ) {

            push(
                @{
                    $parameter_href->{dynamic_parameter}
                      { $pedigree_sample_href->{$key} }
                },
                $sample_id
            );

            ## Reformat to plink format
            if ( $pedigree_sample_href->{$key} eq "unaffected" ) {

                $parameter_href->{dynamic_parameter}{$sample_id}
                  {plink_phenotype} = 1;
            }
            elsif ( $pedigree_sample_href->{$key} eq "affected" ) {

                $parameter_href->{dynamic_parameter}{$sample_id}
                  {plink_phenotype} = 2;
            }
            else {

                $parameter_href->{dynamic_parameter}{$sample_id}
                  {plink_phenotype} = 0;
            }
        }
    }

    ## Add analysis_type for each individual
    if ( $sample_info_href->{sample}{$sample_id}{analysis_type} )
    {    #Add analysis_type

        if ( !$user_supply_switch_href->{analysis_type} ) {

            my $analysis_type =
              $sample_info_href->{sample}{$sample_id}{analysis_type};    #Alias
            $active_parameter_href->{analysis_type}{$sample_id} =
              $analysis_type;
        }
    }

    ## Add expected_coverage for each individual
    if ( $sample_info_href->{sample}{$sample_id}{expected_coverage} )
    {    #Add expected_coverage

        if ( !$user_supply_switch_href->{expected_coverage} ) {

            my $expected_coverage =
              $sample_info_href->{sample}{$sample_id}{expected_coverage}; #Alias
            $active_parameter_href->{expected_coverage}{$sample_id} =
              $expected_coverage;
        }
    }

    ## Add capture kit for each individual
    if ( ( $sample_info_href->{sample}{$sample_id}{capture_kit} ) ) {

        my $capture_kit =
          $sample_info_href->{sample}{$sample_id}{capture_kit};           #Alias

        ## Return a capture kit depending on user info
        my $exome_target_bed_file = add_capture_kit(
            {
                file_info_href => $file_info_href,
                supported_capture_kit_href =>
                  $active_parameter_href->{supported_capture_kit},
                capture_kit => $capture_kit,
                user_supplied_parameter_switch =>
                  $user_supply_switch_href->{exome_target_bed},
            }
        );

        if ($exome_target_bed_file) {

            push(
                @{
                    $exom_target_bed_test_file_tracker_href
                      ->{$exome_target_bed_file}
                },
                $sample_id
            );

        }
    }
}

sub check_founder_id {

##check_founder_id

##Function : Check that founder_ids are included in the pedigree info
##Returns  : ""
##Arguments: $pedigree_href, $pedigree_sample_ids_ref
##         : $pedigree_href           => Pedigree info {REF}
##         : $pedigree_sample_ids_ref => Array of pedigree samples {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $pedigree_href;
    my $pedigree_sample_ids_ref;

    my $tmpl = {
        pedigree_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$pedigree_href
        },
        pedigree_sample_ids_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$pedigree_sample_ids_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

  SAMPLE:
    foreach my $pedigree_sample_href ( @{ $pedigree_href->{samples} } ) {

        my @founders =
          ( $pedigree_sample_href->{father}, $pedigree_sample_href->{mother} );

      FOUNDER:
        foreach my $founder (@founders) {

            if ($founder) {

                if ( !( any { $_ eq $founder } @$pedigree_sample_ids_ref ) )
                {    #If element is not part of array

                    $log->fatal( "Could not find founder sample_id: "
                          . $founder
                          . " in pedigree file\n" );
                    exit 1;
                }
            }
        }
    }
}

sub update_mip_reference_path {

##update_mip_reference_path

##Function : Update path for supplied reference(s) associated with parameter that should reside in the mip reference directory to full path.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $parameter_name
##         : $parameter_href        => Parameter hash {REF}
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $parameter_name        => Parameter to update

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $parameter_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        parameter_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $reference_dir_ref = \$active_parameter_href->{reference_dir};

    if ( $parameter_href->{$parameter_name}{data_type} eq "SCALAR"
        && defined( $active_parameter_href->{$parameter_name} ) )
    {

        my ( $volume, $directory, $file_name ) =
          File::Spec->splitpath( $active_parameter_href->{$parameter_name} )
          ;    #Split to restate
        $active_parameter_href->{$parameter_name} = $file_name
          ;  #Restate to allow for changing mip reference directory between runs
        $active_parameter_href->{$parameter_name} =
          catfile( $$reference_dir_ref,
            $active_parameter_href->{$parameter_name} );  #Update original value
    }
    elsif ( $parameter_href->{$parameter_name}{data_type} eq "ARRAY" ) {

        foreach my $file ( @{ $active_parameter_href->{$parameter_name} } ) {

            my ( $volume, $directory, $file_name ) =
              File::Spec->splitpath($file);               #Split to restate
            $file = catfile( $$reference_dir_ref, $file_name )
              ;    #Update original element
        }
    }
    elsif ( $parameter_href->{$parameter_name}{data_type} eq "HASH" ) {

        foreach my $file ( keys %{ $active_parameter_href->{$parameter_name} } )
        {

            my ( $volume, $directory, $file_name ) =
              File::Spec->splitpath($file);    #Split to restate
            $active_parameter_href->{$parameter_name}
              { catfile( $$reference_dir_ref, $file_name ) } =
              delete( $active_parameter_href->{$parameter_name}{$file} )
              ;                                #Update original value
        }
    }
}

sub check_vcfanno_toml {

##check_vcfanno_toml

##Function : Check that the supplied vcfanno toml frequency file match record 'file=' within toml config file
##Returns  : ""
##Arguments: $vcfanno_file_toml, $vcfanno_file_freq
##         : $vcfanno_file_toml => Toml config file
##         : $vcfanno_file_freq => Frequency file recorded inside toml file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vcfanno_file_toml;
    my $vcfanno_file_freq;

    my $tmpl = {
        vcfanno_file_toml => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$vcfanno_file_toml
        },
        vcfanno_file_freq => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$vcfanno_file_freq
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    open( $FILEHANDLE, "<", $vcfanno_file_toml )
      or
      $log->logdie( "Can't open '" . $vcfanno_file_toml . "': " . $! . "\n" );

    while (<$FILEHANDLE>) {

        chomp $_;                          #Remove newline

        if ( $_ =~ /^file="(\S+)"/ ) {

            my $file_path_freq = $1;

            if ( $file_path_freq ne $vcfanno_file_freq ) {

                $log->fatal( "The supplied vcfanno_config_file: "
                      . $vcfanno_file_freq
                      . " does not match record 'file="
                      . $file_path_freq
                      . "' in the sv_vcfanno_config file: "
                      . $vcfanno_file_toml );
                exit 1;
            }
            last;
        }
    }
    close $FILEHANDLE;
}

sub check_snpsift_keys {

##check_snpsift_keys

##Function : Check that the supplied
##Returns  : ""
##Arguments: $snpsift_annotation_files_href, $snpsift_annotation_outinfo_key_href
##         : $snpsift_annotation_files_href       => Snpsift annotation files {REF}
##         : $snpsift_annotation_outinfo_key_href => File and outinfo key to add to vcf {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $snpsift_annotation_files_href;
    my $snpsift_annotation_outinfo_key_href;

    my $tmpl = {
        snpsift_annotation_files_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$snpsift_annotation_files_href
        },
        snpsift_annotation_outinfo_key_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$snpsift_annotation_outinfo_key_href
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    foreach my $file ( keys %$snpsift_annotation_outinfo_key_href ) {

        unless ( exists( $snpsift_annotation_files_href->{$file} ) ) {

            $log->fatal( "The supplied snpsift_annotation_outinfo_key file: "
                  . $file
                  . " does not match any file in '--snpsift_annotation_files'"
            );
            $log->fatal( "Supplied snpsift_annotation_files files:\n"
                  . join( "\n", keys %$snpsift_annotation_files_href ) );
            exit 1;
        }
    }
}

sub check_key_exists_in_hash {

##check_key_exists_in_hash

##Function : Test if key from query hash exists truth hash
##Returns  : ""
##Arguments: $truth_href, $query_href, $parameter_name
##         : $truth_href     => Truth hash
##         : $query_href     => Query hash
##         : $parameter_name => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $truth_href;
    my $query_href;
    my $parameter_name;

    my $tmpl = {
        truth_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$truth_href
        },
        query_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$query_href
        },
        parameter_name =>
          { required => 1, defined => 1, store => \$parameter_name },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    foreach my $key ( keys %$query_href ) {

        if ( !exists( $truth_href->{$key} ) ) {

            $log->fatal( $parameter_name
                  . " key '"
                  . $key
                  . "' - Does not exist as module program parameter in MIP" );
            exit 1;
        }
    }
}

sub check_element_exists_in_hash {

##check_element_exists_in_hash

##Function : Test if element from query array exists truth hash
##Returns  : ""
##Arguments: $truth_href, $queryies, $parameter_name
##         : $truth_href     => Truth hash
##         : $queryies       => Query array
##         : $parameter_name => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $truth_href;
    my $queryies;
    my $parameter_name;

    my $tmpl = {
        truth_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$truth_href
        },
        queryies => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$queryies
        },
        parameter_name =>
          { required => 1, defined => 1, store => \$parameter_name },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    foreach my $element (@$queryies) {

        if ( !exists( $truth_href->{$element} ) ) {

            $log->fatal( $parameter_name
                  . " element '"
                  . $element
                  . "' - Does not exist as module program parameter in MIP" );
            exit 1;
        }
    }
}

##Investigate potential autodie error
if ( $@ and $@->isa("autodie::exception") ) {

    if ( $@->matches("default") ) {

        say "Not an autodie error at all";
    }
    if ( $@->matches("open") ) {

        say "Error from open";
    }
    if ( $@->matches(":io") ) {

        say "Non-open, IO error.\n";
    }
}
elsif ($@) {

    say "A non-autodie exception.";
}
