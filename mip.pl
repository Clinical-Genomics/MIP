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
use Cwd;
use Cwd qw{ abs_path };
use File::Basename qw{ dirname basename fileparse };
use File::Spec::Functions qw{ catdir catfile devnull splitpath };
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
  qw{ check_allowed_temp_directory check_cmd_config_vs_definition_file check_parameter_hash };
use MIP::Check::Path qw{ check_target_bed_file_exist check_parameter_files };
use MIP::Check::Reference
  qw{ check_bwa_prerequisites check_capture_file_prerequisites check_file_endings_to_build check_human_genome_file_endings check_human_genome_prerequisites check_parameter_metafiles check_references_for_vt };
use MIP::File::Format::Pedigree
  qw{ create_fam_file parse_yaml_pedigree_file reload_previous_pedigree_info };
use MIP::File::Format::Yaml qw{ load_yaml write_yaml order_parameter_names };
use MIP::Get::Analysis qw{ get_overall_analysis_type print_program };
use MIP::Get::File qw{ get_select_file_contigs };
use MIP::Log::MIP_log4perl qw{ initiate_logger set_default_log4perl_file };
use MIP::QC::Record qw{ add_gene_panel add_most_complete_vcf };
use MIP::Script::Utils qw{ help };
use MIP::Set::Contigs qw{ set_contigs };
use MIP::Set::Parameter
  qw{ set_config_to_active_parameters set_default_config_dynamic_parameters set_dynamic_parameter set_parameter_reference_dir_path set_parameter_to_broadcast };
use MIP::Update::Contigs qw{ update_contigs_for_run };
use MIP::Update::Parameters
  qw{ update_dynamic_config_parameters update_reference_parameters update_vcfparser_outfile_counter };
use MIP::Update::Path qw{ update_to_absolute_path };
use MIP::Update::Programs
  qw{ update_program_mode_with_dry_run_all update_program_mode update_prioritize_flag };

## Recipes
use MIP::Recipes::Analysis::Gzip_fastq qw{ analysis_gzip_fastq };
use MIP::Recipes::Analysis::Split_fastq_file qw{ analysis_split_fastq_file };
use MIP::Recipes::Analysis::Vt_core qw{ analysis_vt_core };
use MIP::Recipes::Pipeline::Wgs qw{ pipeline_wgs };
use MIP::Recipes::Pipeline::Wts qw{ pipeline_wts };

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

### %parameter holds all defined parameters for MIP
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
check_parameter_hash(
    {
        parameter_href         => \%parameter,
        mandatory_key_href     => \%mandatory_key,
        non_mandatory_key_href => \%non_mandatory_key,
        file_path              => $definitions_file,
    }
);

## Set MIP version
our $VERSION = 'v5.0.12';

## Holds all active parameters
my %active_parameter;

## Directories, files, job_ids and sample_info
my ( %infile, %indir_path, %infile_lane_prefix, %lane,
    %infile_both_strands_prefix, %job_id, %sample_info );

my %file_info = (
    exome_target_bed =>
      [qw{ .infile_list .pad100.infile_list .pad100.interval_list }],

    # BWA human genome reference file endings
    bwa_build_reference => [qw{ .amb .ann .bwt .pac .sa }],

    # Human genome meta files
    human_genome_reference_file_endings => [qw{ .dict .fai }],
);

#### Staging Area
### Get and/or set input parameters

## Print usage help if no arguments
if ( not @ARGV ) {

    help(
        {
            USAGE     => $USAGE,
            exit_code => 0,
        }
    );
}

### User Options
GetOptions(
    q{ifd|infile_dirs:s}   => \%{ $active_parameter{infile_dirs} },
    q{rd|reference_dir:s}  => \$active_parameter{reference_dir},
    q{p|project_id:s}      => \$active_parameter{project_id},
    q{s|sample_ids:s}      => \@{ $active_parameter{sample_ids} },
    q{odd|outdata_dir:s}   => \$active_parameter{outdata_dir},
    q{osd|outscript_dir:s} => \$active_parameter{outscript_dir},
    q{f|family_id:s}       => \$active_parameter{family_id},
    q{sck|supported_capture_kit:s} =>
      \%{ $active_parameter{supported_capture_kit} },
    q{dnr|decompose_normalize_references:s} =>
      \@{ $active_parameter{decompose_normalize_references} },
    q{ped|pedigree_file:s} => \$active_parameter{pedigree_file},
    q{hgr|human_genome_reference:s} =>
      \$active_parameter{human_genome_reference},
    q{al|outaligner_dir:s}    => \$active_parameter{outaligner_dir},
    q{at|analysis_type:s}     => \%{ $active_parameter{analysis_type} },
    q{pl|platform:s}          => \$active_parameter{platform},
    q{ec|expected_coverage:s} => \%{ $active_parameter{expected_coverage} },
    q{sao|sample_origin:s}    => \%{ $active_parameter{sample_origin} },
    q{c|config_file:s}        => \$active_parameter{config_file},
    q{ccp|cluster_constant_path:s} => \$active_parameter{cluster_constant_path},
    q{acp|analysis_constant_path:s} =>
      \$active_parameter{analysis_constant_path},
    q{cfa|config_file_analysis:s} => \$active_parameter{config_file_analysis},
    q{sif|sample_info_file:s}     => \$active_parameter{sample_info_file},
    q{dra|dry_run_all=i}          => \$active_parameter{dry_run_all},
    q{jul|java_use_large_pages=n} => \$active_parameter{java_use_large_pages},
    q{ges|genomic_set:s}          => \$active_parameter{genomic_set},
    q{rio|reduce_io=n}            => \$active_parameter{reduce_io},
    q{riu|replace_iupac=n}        => \$active_parameter{replace_iupac},
    q{ppm|print_program_mode=n}   => \$active_parameter{print_program_mode},
    q{pp|print_program}           => sub {
        ## Force ppm to be read before function call
        GetOptions( q{ppm|print_program_mode=n} =>
              \$active_parameter{print_program_mode} );
        print_program(
            {
                parameter_href     => \%parameter,
                print_program_mode => $active_parameter{print_program_mode},
            }
        );
        exit;
    },
    q{l|log_file:s} => \$active_parameter{log_file},
    q{h|help}       => sub { say STDOUT $USAGE; exit; },
    q{v|version}    => sub {
        say STDOUT $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $VERSION,
          $NEWLINE;
        exit;
    },

    #### Bash
    q{bse|bash_set_errexit:s}   => \$active_parameter{bash_set_errexit},
    q{bsu|bash_set_nounset:s}   => \$active_parameter{bash_set_nounset},
    q{bsp|bash_set_pipefail:s}  => \$active_parameter{bash_set_pipefail},
    q{em|email:s}               => \$active_parameter{email},
    q{emt|email_types:s}        => \@{ $active_parameter{email_types} },
    q{mcn|module_core_number:s} => \%{ $active_parameter{module_core_number} },
    q{mot|module_time:s}        => \%{ $active_parameter{module_time} },
    q{mse|module_source_environment_command:s} =>
      \%{ $active_parameter{module_source_environment_command} },
    q{sen|source_main_environment_commands=s{,}} =>
      \@{ $active_parameter{source_main_environment_commands} },
    q{mcn|max_cores_per_node=n} => \$active_parameter{max_cores_per_node},
    q{nrm|node_ram_memory=n}    => \$active_parameter{node_ram_memory},
    q{tmd|temp_directory:s}     => \$active_parameter{temp_directory},
    q{qos|slurm_quality_of_service=s} =>
      \$active_parameter{slurm_quality_of_service},
    q{psfq|psplit_fastq_file=n} => \$active_parameter{psplit_fastq_file},
    q{sfqrdb|split_fastq_file_read_batch=n} =>
      \$active_parameter{split_fastq_file_read_batch},
    q{pgz|pgzip_fastq=n}         => \$active_parameter{pgzip_fastq},
    q{pfqc|pfastqc=n}            => \$active_parameter{pfastqc},
    q{pcta|pcutadapt=n}          => \$active_parameter{pcutadapt},
    q{pmem|pbwa_mem=n}           => \$active_parameter{pbwa_mem},
    q{memhla|bwa_mem_hla=n}      => \$active_parameter{bwa_mem_hla},
    q{memcrm|bwa_mem_cram=n}     => \$active_parameter{bwa_mem_cram},
    q{memsts|bwa_mem_bamstats=n} => \$active_parameter{bwa_mem_bamstats},
    q{memssm|bwa_sambamba_sort_memory_limit:s} =>
      \$active_parameter{bwa_sambamba_sort_memory_limit},
    q{ptp|picardtools_path:s} => \$active_parameter{picardtools_path},
    q{pptm|ppicardtools_mergesamfiles=n} =>
      \$active_parameter{ppicardtools_mergesamfiles},
    q{pmd|pmarkduplicates=n} => \$active_parameter{pmarkduplicates},
    q{mdpmd|markduplicates_picardtools_markduplicates=n} =>
      \$active_parameter{markduplicates_picardtools_markduplicates},
    q{mdsmd|markduplicates_sambamba_markdup=n} =>
      \$active_parameter{markduplicates_sambamba_markdup},
    q{mdshts|markduplicates_sambamba_markdup_hash_table_size=n} =>
      \$active_parameter{markduplicates_sambamba_markdup_hash_table_size},
    q{mdsols|markduplicates_sambamba_markdup_overflow_list_size=n} =>
      \$active_parameter{markduplicates_sambamba_markdup_overflow_list_size},
    q{mdsibs|markduplicates_sambamba_markdup_io_buffer_size=n} =>
      \$active_parameter{markduplicates_sambamba_markdup_io_buffer_size},
    q{pchs|pchanjo_sexcheck=n} => \$active_parameter{pchanjo_sexcheck},
    q{chslle|chanjo_sexcheck_log_level:s} =>
      \$active_parameter{chanjo_sexcheck_log_level},
    q{psdt|psambamba_depth=n}       => \$active_parameter{psambamba_depth},
    q{sdtmod|sambamba_depth_mode:s} => \$active_parameter{sambamba_depth_mode},
    q{sdtcut|sambamba_depth_cutoffs:s} =>
      \@{ $active_parameter{sambamba_depth_cutoffs} },
    q{sdtbed|sambamba_depth_bed:s} => \$active_parameter{sambamba_depth_bed},
    q{sdtbaq|sambamba_depth_base_quality=n} =>
      \$active_parameter{sambamba_depth_base_quality},
    q{sdtmaq|sambamba_depth_mapping_quality=n} =>
      \$active_parameter{sambamba_depth_mapping_quality},
    q{sdtndu|sambamba_depth_noduplicates=n} =>
      \$active_parameter{sambamba_depth_noduplicates},
    q{sdtfqc|sambamba_depth_quality_control=n} =>
      \$active_parameter{sambamba_depth_quality_control},
    q{pbgc|pbedtools_genomecov=n} => \$active_parameter{pbedtools_genomecov},
    q{bgcmc|bedtools_genomecov_max_coverage=n} =>
      \$active_parameter{bedtools_genomecov_max_coverage},
    q{pptcmm|ppicardtools_collectmultiplemetrics=n} =>
      \$active_parameter{ppicardtools_collectmultiplemetrics},
    q{pptchs|ppicardtools_collecthsmetrics=n} =>
      \$active_parameter{ppicardtools_collecthsmetrics},
    q{extb|exome_target_bed=s}     => \%{ $active_parameter{exome_target_bed} },
    q{prcp|prcovplots=n}           => \$active_parameter{prcovplots},
    q{pcnv|pcnvnator=n}            => \$active_parameter{pcnvnator},
    q{cnvhbs|cnv_bin_size=n}       => \$active_parameter{cnv_bin_size},
    q{pdelc|pdelly_call=n}         => \$active_parameter{pdelly_call},
    q{pdel|pdelly_reformat=n}      => \$active_parameter{pdelly_reformat},
    q{deltyp|delly_types:s}        => \@{ $active_parameter{delly_types} },
    q{delexc|delly_exclude_file:s} => \$active_parameter{delly_exclude_file},
    q{pmna|pmanta=n}               => \$active_parameter{pmanta},
    q{ptid|ptiddit=n}              => \$active_parameter{ptiddit},
    q{tidmsp|tiddit_minimum_number_supporting_pairs=n} =>
      \$active_parameter{tiddit_minimum_number_supporting_pairs},
    q{psvc|psv_combinevariantcallsets=n} =>
      \$active_parameter{psv_combinevariantcallsets},
    q{svcvtd|sv_vt_decompose=n} => \$active_parameter{sv_vt_decompose},
    q{svsvdbmp|sv_svdb_merge_prioritize:s} =>
      \$active_parameter{sv_svdb_merge_prioritize},
    q{svcbtv|sv_bcftools_view_filter=n} =>
      \$active_parameter{sv_bcftools_view_filter},
    q{svcdbq|sv_svdb_query=n} => \$active_parameter{sv_svdb_query},
    q{svcdbqd|sv_svdb_query_db_files:s} =>
      \%{ $active_parameter{sv_svdb_query_db_files} },
    q{svcvan|sv_vcfanno=n}        => \$active_parameter{sv_vcfanno},
    q{svcval|sv_vcfanno_lua:s}    => \$active_parameter{sv_vcfanno_lua},
    q{svcvac|sv_vcfanno_config:s} => \$active_parameter{sv_vcfanno_config},
    q{svcvacf|sv_vcfanno_config_file:s} =>
      \$active_parameter{sv_vcfanno_config_file},
    q{svcvah|sv_vcfannotation_header_lines_file:s} =>
      \$active_parameter{sv_vcfannotation_header_lines_file},
    q{svcgmf|sv_genmod_filter=n} => \$active_parameter{sv_genmod_filter},
    q{svcgfr|sv_genmod_filter_1000g:s} =>
      \$active_parameter{sv_genmod_filter_1000g},
    q{svcgft|sv_genmod_filter_threshold:s} =>
      \$active_parameter{sv_genmod_filter_threshold},
    q{svcbcf|sv_combinevariantcallsets_bcf_file=n} =>
      \$active_parameter{sv_combinevariantcallsets_bcf_file},
    q{psvv|psv_varianteffectpredictor=n} =>
      \$active_parameter{psv_varianteffectpredictor},
    q{svvepf|sv_vep_features:s} => \@{ $active_parameter{sv_vep_features} },
    q{svvepl|sv_vep_plugins:s}  => \@{ $active_parameter{sv_vep_plugins} },
    q{psvvcp|psv_vcfparser=n}   => \$active_parameter{psv_vcfparser},
    q{svvcpvt|sv_vcfparser_vep_transcripts=n} =>
      \$active_parameter{sv_vcfparser_vep_transcripts},
    q{svvcppg|sv_vcfparser_per_gene=n} =>
      \$active_parameter{sv_vcfparser_per_gene},
    q{svvcprff|sv_vcfparser_range_feature_file:s} =>
      \$active_parameter{sv_vcfparser_range_feature_file},
    q{svvcprfa|sv_vcfparser_range_feature_annotation_columns:s} =>
      \@{ $active_parameter{sv_vcfparser_range_feature_annotation_columns} },
    q{svvcpsf|sv_vcfparser_select_file:s} =>
      \$active_parameter{sv_vcfparser_select_file},
    q{svvcpsfm|sv_vcfparser_select_file_matching_column=n} =>
      \$active_parameter{sv_vcfparser_select_file_matching_column},
    q{svvcpsfa|sv_vcfparser_select_feature_annotation_columns:s} =>
      \@{ $active_parameter{sv_vcfparser_select_feature_annotation_columns} },
    q{psvr|psv_rankvariant=n} => \$active_parameter{psv_rankvariant},
    q{svravanr|sv_genmod_annotate_regions:n} =>
      \$active_parameter{sv_genmod_annotate_regions},
    q{svravgft|sv_genmod_models_family_type:s} =>
      \$active_parameter{sv_genmod_models_family_type},
    q{svravrpf|sv_genmod_models_reduced_penetrance_file:s} =>
      \$active_parameter{sv_genmod_models_reduced_penetrance_file},
    q{svravwg|sv_genmod_models_whole_gene=n} =>
      \$active_parameter{sv_genmod_models_whole_gene},
    q{svravrm|sv_rank_model_file:s} => \$active_parameter{sv_rank_model_file},
    q{psvre|psv_reformat=n}         => \$active_parameter{psv_reformat},
    q{svrevbf|sv_rankvariant_binary_file=n} =>
      \$active_parameter{sv_rankvariant_binary_file},
    q{svrergf|sv_reformat_remove_genes_file:s} =>
      \$active_parameter{sv_reformat_remove_genes_file},
    q{pbmp|pbcftools_mpileup=n} => \$active_parameter{pbcftools_mpileup},
    q{pbmpfv|bcftools_mpileup_filter_variant} =>
      \$active_parameter{bcftools_mpileup_filter_variant},
    q{pfrb|pfreebayes=n}        => \$active_parameter{pfreebayes},
    q{gtp|gatk_path:s}          => \$active_parameter{gatk_path},
    q{gll|gatk_logging_level:s} => \$active_parameter{gatk_logging_level},
    q{gbdv|gatk_bundle_download_version:s} =>
      \$active_parameter{gatk_bundle_download_version},
    q{gdco|gatk_downsample_to_coverage=n} =>
      \$active_parameter{gatk_downsample_to_coverage},
    q{gdai|gatk_disable_auto_index_and_file_lock=n} =>
      \$active_parameter{gatk_disable_auto_index_and_file_lock},
    q{pgra|pgatk_realigner=n} => \$active_parameter{pgatk_realigner},
    q{graks|gatk_realigner_indel_known_sites:s} =>
      \@{ $active_parameter{gatk_realigner_indel_known_sites} },
    q{pgbr|pgatk_baserecalibration=n} =>
      \$active_parameter{pgatk_baserecalibration},
    q{gbrcov|gatk_baserecalibration_covariates:s} =>
      \@{ $active_parameter{gatk_baserecalibration_covariates} },
    q{gbrkst|gatk_baserecalibration_known_sites:s} =>
      \@{ $active_parameter{gatk_baserecalibration_known_sites} },
    q{gbrrf|gatk_baserecalibration_read_filters=s} =>
      \@{ $active_parameter{gatk_baserecalibration_read_filters} },
    q{gbrdiq|gatk_baserecalibration_disable_indel_qual=n} =>
      \$active_parameter{gatk_baserecalibration_disable_indel_qual},
    q{gbrsqq|gatk_baserecalibration_static_quantized_quals:s} =>
      \@{ $active_parameter{gatk_baserecalibration_static_quantized_quals} },
    q{pghc|pgatk_haplotypecaller=n} =>
      \$active_parameter{pgatk_haplotypecaller},
    q{ghcann|gatk_haplotypecaller_annotation:s} =>
      \@{ $active_parameter{gatk_haplotypecaller_annotation} },
    q{ghckse|gatk_haplotypecaller_snp_known_set:s} =>
      \$active_parameter{gatk_haplotypecaller_snp_known_set},
    q{ghcscb|gatk_haplotypecaller_no_soft_clipped_bases=n} =>
      \$active_parameter{gatk_haplotypecaller_no_soft_clipped_bases},
    q{ghcpim|gatk_haplotypecaller_pcr_indel_model:s} =>
      \$active_parameter{gatk_haplotypecaller_pcr_indel_model},
    q{pggt|pgatk_genotypegvcfs=n} => \$active_parameter{pgatk_genotypegvcfs},
    q{ggtgrl|gatk_genotypegvcfs_ref_gvcf:s} =>
      \$active_parameter{gatk_genotypegvcfs_ref_gvcf},
    q{ggtals|gatk_genotypegvcfs_all_sites=n} =>
      \$active_parameter{gatk_genotypegvcfs_all_sites},
    q{ggbcf|gatk_concatenate_genotypegvcfs_bcf_file=n} =>
      \$active_parameter{gatk_concatenate_genotypegvcfs_bcf_file},
    q{pgvr|pgatk_variantrecalibration=n} =>
      \$active_parameter{pgatk_variantrecalibration},
    q{gvrann|gatk_variantrecalibration_annotations:s} =>
      \@{ $active_parameter{gatk_variantrecalibration_annotations} },
    q{gvrres|gatk_variantrecalibration_resource_snv=s} =>
      \%{ $active_parameter{gatk_variantrecalibration_resource_snv} },
    q{gvrrei|gatk_variantrecalibration_resource_indel=s} =>
      \%{ $active_parameter{gatk_variantrecalibration_resource_indel} },
    q{gvrstf|gatk_variantrecalibration_snv_tsfilter_level:s} =>
      \$active_parameter{gatk_variantrecalibration_snv_tsfilter_level},
    q{gvritf|gatk_variantrecalibration_indel_tsfilter_level:s} =>
      \$active_parameter{gatk_variantrecalibration_indel_tsfilter_level},
    q{gvrdpa|gatk_variantrecalibration_dp_annotation=n} =>
      \$active_parameter{gatk_variantrecalibration_dp_annotation},
    q{gvrsmg|gatk_variantrecalibration_snv_max_gaussians=n} =>
      \$active_parameter{gatk_variantrecalibration_snv_max_gaussians},
    q{gvrimg|gatk_variantrecalibration_indel_max_gaussians=n} =>
      \$active_parameter{gatk_variantrecalibration_indel_max_gaussians},
    q{gcgpss|gatk_calculategenotypeposteriors_support_set:s} =>
      \$active_parameter{gatk_calculategenotypeposteriors_support_set},
    q{pgcv|pgatk_combinevariantcallsets=n} =>
      \$active_parameter{pgatk_combinevariantcallsets},
    q{gcvgmo|gatk_combinevariants_genotype_merge_option:s} =>
      \$active_parameter{gatk_combinevariants_genotype_merge_option},
    q{gcvpc|gatk_combinevariants_prioritize_caller:s} =>
      \$active_parameter{gatk_combinevariants_prioritize_caller},
    q{gcvbcf|gatk_combinevariantcallsets_bcf_file=n} =>
      \$active_parameter{gatk_combinevariantcallsets_bcf_file},
    q{pgvea|pgatk_variantevalall=n} => \$active_parameter{pgatk_variantevalall},
    q{pgvee|pgatk_variantevalexome=n} =>
      \$active_parameter{pgatk_variantevalexome},
    q{gveedbs|gatk_varianteval_dbsnp:s} =>
      \$active_parameter{gatk_varianteval_dbsnp},
    q{gveedbg|gatk_varianteval_gold:s} =>
      \$active_parameter{gatk_varianteval_gold},
    q{ppvab|pprepareforvariantannotationblock=n} =>
      \$active_parameter{pprepareforvariantannotationblock},
    q{prhc|prhocall=n} => \$active_parameter{prhocall},
    q{rhcf|rhocall_frequency_file:s} =>
      \$active_parameter{rhocall_frequency_file},
    q{pvt|pvt=n}             => \$active_parameter{pvt},
    q{vtddec|vt_decompose=n} => \$active_parameter{vt_decompose},
    q{vtdnor|vt_normalize=n} => \$active_parameter{vt_normalize},
    q{vtunq|vt_uniq=n}       => \$active_parameter{vt_uniq},
    q{vtmaa|vt_missing_alt_allele=n} =>
      \$active_parameter{vt_missing_alt_allele},
    q{pfqf|pfrequency_filter=n} => \$active_parameter{pfrequency_filter},
    q{fqfgmf|frequency_genmod_filter=n} =>
      \$active_parameter{frequency_genmod_filter},
    q{fqfgfr|frequency_genmod_filter_1000g:s} =>
      \$active_parameter{frequency_genmod_filter_1000g},
    q{fqfmaf|frequency_genmod_filter_max_af=n} =>
      \$active_parameter{frequency_genmod_filter_max_af},
    q{fqfgft|frequency_genmod_filter_threshold:s} =>
      \$active_parameter{frequency_genmod_filter_threshold},
    q{pvep|pvarianteffectpredictor=n} =>
      \$active_parameter{pvarianteffectpredictor},
    q{vepp|vep_directory_path:s}  => \$active_parameter{vep_directory_path},
    q{vepc|vep_directory_cache:s} => \$active_parameter{vep_directory_cache},
    q{vepr|vep_reference:n}       => \$active_parameter{vep_reference},
    q{vepf|vep_features:s}        => \@{ $active_parameter{vep_features} },
    q{veppl|vep_plugins:s}        => \@{ $active_parameter{vep_plugins} },
    q{pvcp|pvcfparser=n}          => \$active_parameter{pvcfparser},
    q{vcpvt|vcfparser_vep_transcripts=n} =>
      \$active_parameter{vcfparser_vep_transcripts},
    q{vcprff|vcfparser_range_feature_file:s} =>
      \$active_parameter{vcfparser_range_feature_file},
    q{vcprfa|vcfparser_range_feature_annotation_columns:s} =>
      \@{ $active_parameter{vcfparser_range_feature_annotation_columns} },
    q{vcpsf|vcfparser_select_file:s} =>
      \$active_parameter{vcfparser_select_file},
    q{vcpsfm|vcfparser_select_file_matching_column=n} =>
      \$active_parameter{vcfparser_select_file_matching_column},
    q{vcpsfa|vcfparser_select_feature_annotation_columns:s} =>
      \@{ $active_parameter{vcfparser_select_feature_annotation_columns} },
    q{psne|psnpeff=n}      => \$active_parameter{psnpeff},
    q{snep|snpeff_path:s}  => \$active_parameter{snpeff_path},
    q{sneann|snpeff_ann=n} => \$active_parameter{snpeff_ann},
    q{snegbv|snpeff_genome_build_version:s} =>
      \$active_parameter{snpeff_genome_build_version},
    q{snesaf|snpsift_annotation_files=s} =>
      \%{ $active_parameter{snpsift_annotation_files} },
    q{snesaoi|snpsift_annotation_outinfo_key=s} =>
      \%{ $active_parameter{snpsift_annotation_outinfo_key} },
    q{snesdbnsfp|snpsift_dbnsfp_file:s} =>
      \$active_parameter{snpsift_dbnsfp_file},
    q{snesdbnsfpa|snpsift_dbnsfp_annotations:s} =>
      \@{ $active_parameter{snpsift_dbnsfp_annotations} },
    q{prav|prankvariant=n} => \$active_parameter{prankvariant},
    q{ravgft|genmod_models_family_type:s} =>
      \$active_parameter{genmod_models_family_type},
    q{ravanr|genmod_annotate_regions:n} =>
      \$active_parameter{genmod_annotate_regions},
    q{ravcad|genmod_annotate_cadd_files:s} =>
      \@{ $active_parameter{genmod_annotate_cadd_files} },
    q{ravspi|genmod_annotate_spidex_file:s} =>
      \$active_parameter{genmod_annotate_spidex_file},
    q{ravwg|genmod_models_whole_gene=n} =>
      \$active_parameter{genmod_models_whole_gene},
    q{ravrpf|genmod_models_reduced_penetrance_file:s} =>
      \$active_parameter{genmod_models_reduced_penetrance_file},
    q{ravrm|rank_model_file:s} => \$active_parameter{rank_model_file},
    q{pevab|pendvariantannotationblock=n} =>
      \$active_parameter{pendvariantannotationblock},
    q{evabrgf|endvariantannotationblock_remove_genes_file:s} =>
      \$active_parameter{endvariantannotationblock_remove_genes_file},
    q{ravbf|rankvariant_binary_file=n} =>
      \$active_parameter{rankvariant_binary_file},
    q{pped|ppeddy=n}             => \$active_parameter{ppeddy},
    q{pplink|pplink=n}           => \$active_parameter{pplink},
    q{pvai|pvariant_integrity=n} => \$active_parameter{pvariant_integrity},
    q{pevl|pevaluation=n}        => \$active_parameter{pevaluation},
    q{evlnid|nist_id:s}          => \$active_parameter{nist_id},
    q{evlnhc|nist_high_confidence_call_set:s} =>
      \$active_parameter{nist_high_confidence_call_set},
    q{evlnil|nist_high_confidence_call_set_bed:s} =>
      \$active_parameter{nist_high_confidence_call_set_bed},
    q{pqcc|pqccollect=n} => \$active_parameter{pqccollect},
    q{qccsi|qccollect_sampleinfo_file:s} =>
      \$active_parameter{qccollect_sampleinfo_file},
    q{qccref|qccollect_regexp_file:s} =>
      \$active_parameter{qccollect_regexp_file},
    q{qccske|qccollect_skip_evaluation} =>
      \$active_parameter{qccollect_skip_evaluation},
    q{pmqc|pmultiqc=n}           => \$active_parameter{pmultiqc},
    q{pars|panalysisrunstatus=n} => \$active_parameter{panalysisrunstatus},
    q{psac|psacct=n}             => \$active_parameter{psacct},
    q{sacfrf|sacct_format_fields:s} =>
      \@{ $active_parameter{sacct_format_fields} },
  )
  or help(
    {
        USAGE     => $USAGE,
        exit_code => 1,
    }
  );

##Special case:Enable/activate MIP. Cannot be changed from cmd or config
$active_parameter{mip} = $parameter{mip}{default};

## Change relative path to absolute path for parameter with "update_path: absolute_path" in config
update_to_absolute_path(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter,
    }
);

### Config file
## If config from cmd
if ( exists $active_parameter{config_file}
    && defined $active_parameter{config_file} )
{

    ## Loads a YAML file into an arbitrary hash and returns it.
    my %config_parameter =
      load_yaml( { yaml_file => $active_parameter{config_file}, } );

    ## Remove previous analysis specific info not relevant for current run e.g. log file, which is read from pedigree or cmd
    my @remove_keys = (qw{ log_file });

  KEY:
    foreach my $key (@remove_keys) {

        delete $config_parameter{$key};
    }

## Set config parameters into %active_parameter unless $parameter
## has been supplied on the command line
    set_config_to_active_parameters(
        {
            config_parameter_href => \%config_parameter,
            active_parameter_href => \%active_parameter,
        }
    );

    ## Compare keys from config and cmd (%active_parameter) with definitions file (%parameter)
    check_cmd_config_vs_definition_file(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
        }
    );

    my @config_dynamic_parameters = qw{ analysis_constant_path outaligner_dir };

    ## Replace config parameter with cmd info for config dynamic parameter
    set_default_config_dynamic_parameters(
        {
            parameter_href        => \%parameter,
            active_parameter_href => \%active_parameter,
            parameter_names_ref   => \@config_dynamic_parameters,
        }
    );

    ## Loop through all parameters and update info
  PARAMETER:
    foreach my $parameter_name (@order_parameters) {

        ## Updates the active parameters to particular user/cluster for dynamic config parameters following specifications. Leaves other entries untouched.
        update_dynamic_config_parameters(
            {
                active_parameter_href => \%active_parameter,
                parameter_name        => $parameter_name,
            }
        );
    }
}

## Set the default Log4perl file using supplied dynamic parameters.
$active_parameter{log_file} = set_default_log4perl_file(
    {
        active_parameter_href => \%active_parameter,
        cmd_input             => $active_parameter{log_file},
        script                => $script,
        date                  => $date,
        date_time_stamp       => $date_time_stamp,
    }
);

## Creates log object
my $log = initiate_logger(
    {
        file_path => $active_parameter{log_file},
        log_name  => q{MIP},
    }
);

## Parse pedigree file
## Reads family_id_pedigree file in YAML format. Checks for pedigree data for allowed entries and correct format. Add data to sample_info depending on user info.
# Meta data in YAML format
if ( defined $active_parameter{pedigree_file} ) {

    ## Loads a YAML file into an arbitrary hash and returns it. Load parameters from previous run from sample_info_file
    my %pedigree =
      load_yaml( { yaml_file => $active_parameter{pedigree_file}, } );

    $log->info( q{Loaded: } . $active_parameter{pedigree_file} );

    parse_yaml_pedigree_file(
        {
            parameter_href        => \%parameter,
            active_parameter_href => \%active_parameter,
            sample_info_href      => \%sample_info,
            file_path             => $active_parameter{pedigree_file},
            pedigree_href         => \%pedigree,
        }
    );
}

## Detect if all samples has the same sequencing type and return consensus if reached
$parameter{dynamic_parameter}{consensus_analysis_type} =
  get_overall_analysis_type(
    { analysis_type_href => \%{ $active_parameter{analysis_type} }, } );

### Populate uninitilized active_parameters{parameter_name} with default from parameter
PARAMETER:
foreach my $parameter_name (@order_parameters) {

    ## If hash and set - skip
    next PARAMETER
      if ( ref $active_parameter{$parameter_name} eq qw{HASH}
        && keys %{ $active_parameter{$parameter_name} } );

    ## If array and set - skip
    next PARAMETER
      if ( ref $active_parameter{$parameter_name} eq qw{ARRAY}
        && @{ $active_parameter{$parameter_name} } );

    ## If scalar and set - skip
    next PARAMETER
      if ( defined $active_parameter{$parameter_name}
        and not ref $active_parameter{$parameter_name} );

    ### Special case for parameters that are dependent on other parameters values
    my @custom_default_parameters =
      qw{ exome_target_bed bwa_build_reference analysis_type infile_dirs sample_info_file };

    if ( any { $_ eq $parameter_name } @custom_default_parameters ) {

        set_custom_default_to_active_parameter(
            {
                parameter_href        => \%parameter,
                active_parameter_href => \%active_parameter,
                file_info_href        => \%file_info,
                parameter_name        => $parameter_name,
            }
        );
        next PARAMETER;
    }

    ## Checks and sets user input or default values to active_parameters
    set_default_to_active_parameter(
        {
            parameter_href        => \%parameter,
            active_parameter_href => \%active_parameter,
            associated_programs_ref =>
              \@{ $parameter{$parameter_name}{associated_program} },
            parameter_name => $parameter_name,
        }
    );
}

## Update path for supplied reference(s) associated with parameter that should reside in the mip reference directory to full path
set_parameter_reference_dir_path(
    {
        active_parameter_href => \%active_parameter,
        parameter_name        => q{human_genome_reference},
    }
);

## Detect version and source of the human_genome_reference: Source (hg19 or GRCh).
parse_human_genome_reference(
    {
        file_info_href => \%file_info,
        human_genome_reference_ref =>
          \basename( $active_parameter{human_genome_reference} ),
    }
);

## Update exome_target_bed files with human_genome_reference_source_ref and human_genome_reference_version_ref
update_exome_target_bed(
    {
        exome_target_bed_file_href => $active_parameter{exome_target_bed},
        human_genome_reference_source_ref =>
          \$file_info{human_genome_reference_source},
        human_genome_reference_version_ref =>
          \$file_info{human_genome_reference_version},
    }
);

# Holds all active parameters values for broadcasting
my @broadcasts;

set_parameter_to_broadcast(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter,
        order_parameters_ref  => \@order_parameters,
        broadcasts_ref        => \@broadcasts,
    }
);

## Reference in MIP reference directory
PARAMETER:
foreach my $parameter_name ( keys %parameter ) {

    ## Expect file to be in reference directory
    if ( exists $parameter{$parameter_name}{reference} ) {

        update_reference_parameters(
            {
                active_parameter_href => \%active_parameter,
                associated_programs_ref =>
                  \@{ $parameter{$parameter_name}{associated_program} },
                parameter_name => $parameter_name,
            }
        );
    }
}

### Checks

## Check existance of files and directories
PARAMETER:
foreach my $parameter_name ( keys %parameter ) {

    if ( exists $parameter{$parameter_name}{exists_check} ) {

        check_parameter_files(
            {
                parameter_href        => \%parameter,
                active_parameter_href => \%active_parameter,
                associated_programs_ref =>
                  \@{ $parameter{$parameter_name}{associated_program} },
                parameter_name => $parameter_name,
                parameter_exists_check =>
                  $parameter{$parameter_name}{exists_check},
            }
        );
    }
}

## Updates sample_info hash with previous run pedigree info
reload_previous_pedigree_info(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
        sample_info_file_path => $active_parameter{sample_info_file},
    }
);

## Special case since dict is created with .fastq removed
## Check the existance of associated human genome files
check_human_genome_file_endings(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        parameter_name        => q{human_genome_reference},
    }
);

## Check that supplied target file ends with ".bed" and otherwise croaks
TARGET_FILE:
foreach my $target_bed_file ( keys %{ $active_parameter{exome_target_bed} } ) {

    check_target_bed_file_exist(
        {
            parameter_name => q{exome_target_bed},
            path           => $target_bed_file,
        }
    );
}

## Checks parameter metafile exists
check_parameter_metafiles(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
    }
);

## Update the expected number of outfile after vcfparser
update_vcfparser_outfile_counter(
    { active_parameter_href => \%active_parameter, } );

## Collect select file contigs to loop over downstream
if ( $active_parameter{vcfparser_select_file} ) {

## Collects sequences contigs used in select file
    @{ $file_info{select_file_contigs} } = get_select_file_contigs(
        {
            select_file_path =>
              catfile( $active_parameter{vcfparser_select_file} ),
            log => $log,
        }
    );
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
if ( defined $active_parameter{email} ) {    #Allow no malformed email adress

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
my @parameter_keys_to_check =
  (qw{ module_time module_core_number module_source_environment_command });
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
        parameter_names_ref =>
          [qw{ analysis_type expected_coverage sample_origin }],
        parameter_href => \%parameter,
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

    my %no_dry_run_info = (
        analysisrunstatus => q{not_finished},
        analysis_date     => $date_time_stamp,
        mip_version       => $VERSION,
    );

  KEY_VALUE_PAIR:
    while ( my ( $key, $value ) = each %no_dry_run_info ) {

        $sample_info{$key} = $value;
    }
}

### Build recipes
$log->info(q{[Reference check - Reference prerequisites]});
## Check capture file prerequistes exists
PROGRAM:
foreach
  my $program_name ( @{ $parameter{exome_target_bed}{associated_program} } )
{

    next PROGRAM if ( not $active_parameter{$program_name} );

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

## Check human genome prerequistes exists
PROGRAM:
foreach my $program_name (
    @{ $parameter{human_genome_reference}{associated_program} } )
{

    next PROGRAM if ( $program_name eq q{mip} );

    next PROGRAM if ( not $active_parameter{$program_name} );

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

## Check BWA build prerequisites

if ( $active_parameter{pbwa_mem} ) {

    check_bwa_prerequisites(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            sample_info_href        => \%sample_info,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            program_name            => q{bwa_mem},
            parameter_build_name    => q{bwa_build_reference},
        }
    );
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
if ( $active_parameter{psplit_fastq_file} ) {

    $log->info(q{[Split fastq files in batches]});

  SAMPLE_ID:
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
                program_name          => q{split_fastq_file},
                sequence_read_batch =>
                  $active_parameter{split_fastq_file_read_batch},
            }
        );
    }

    ## End here if this module is turned on
    exit;
}

## GZip of fastq files
if (   $active_parameter{pgzip_fastq}
    && $uncompressed_file_switch eq q{uncompressed} )
{

    $log->info(q{[Gzip for fastq files]});

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
                        program_name            => q{gzip_fastq},
                    }
                );

                # Call once per sample_id
                last INFILES;
            }
        }
    }
}

my $consensus_analysis_type =
  $parameter{dynamic_parameter}{consensus_analysis_type};

### WTS
if ( $consensus_analysis_type eq q{wts} ) {

    $log->info( q{Pipeline analysis type: } . $consensus_analysis_type );

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
}

### WES|WGS
if (   $consensus_analysis_type eq q{wgs}
    || $consensus_analysis_type eq q{wes} )
{

    $log->info( q{Pipeline analysis type: } . $consensus_analysis_type );

    ## Pipeline recipe for wts data
    pipeline_wgs(
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
}

## Write QC for programs used in analysis
# Write SampleInfo to yaml file
if ( $active_parameter{sample_info_file} ) {

    ## Writes a YAML hash to file
    write_yaml(
        {
            yaml_href      => \%sample_info,
            yaml_file_path => $active_parameter{sample_info_file},
        }
    );
    $log->info( q{Wrote: } . $active_parameter{sample_info_file} );
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
    -dnr/--decompose_normalize_references Set the references to be decomposed and normalized (defaults: "gatk_realigner_indel_known_sites", "gatk_baserecalibration_known_sites","gatk_haplotypecaller_snp_known_set", "gatk_variantrecalibration_resource_snv", "gatk_variantrecalibration_resource_indel", "frequency_genmod_filter_1000g", "sv_vcfanno_config_file", "gatk_varianteval_gold", "gatk_varianteval_dbsnp","snpsift_annotation_files")
    -ped/--pedigree_file Meta data on samples (defaults to "")
    -hgr/--human_genome_reference Fasta file for the human genome reference (defaults to "GRCh37_homo_sapiens_-d5-.fasta;1000G decoy version 5")
    -ald/--outaligner_dir Setting which aligner out directory was used for alignment in previous analysis (defaults to "{outdata_dir}{outaligner_dir}")
    -at/--analysis_type Type of analysis to perform (sample_id=analysis_type, defaults to "wgs";Valid entries: "wgs", "wes", "wts")
    -pl/--platform Platform/technology used to produce the reads (defaults to "ILLUMINA")
    -ec/--expected_coverage Expected mean target coverage for analysis (sample_id=expected_coverage, defaults to "")
    -sao/--sample_origin Sample origin for analysis (sample_id=sample_origin, defaults to "")
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
    -mse/--module_source_environment_command Set environment variables specific for each module (Format: module "program_name"="command"
    -sen/--source_main_environment_commands Source main environment command in sbatch scripts (defaults to "")
    -mcn/--max_cores_per_node The maximum number of processor cores per node used in the analysis (defaults to "16")
    -nrm/--node_ram_memory The RAM memory size of the node(s) in GigaBytes (Defaults to 24)
    -tmd/--temp_directory Set the temporary directory for all programs (defaults to "/scratch/SLURM_JOB_ID";supply whole path)
    -em/--email E-mail (defaults to "")
    -emt/--email_types E-mail type (defaults to FAIL (=FAIL);Options: BEGIN (=BEGIN) and/or F (=FAIL) and/or END=(END))
    -qos/--slurm_quality_of_service SLURM quality of service command in sbatch scripts (defaults to "normal")

    ####Programs
    -psfq/--psplit_fastq_file Split fastq files in batches of X reads and exits (defaults to "0" (=no))
      -sfqrdb/--split_fastq_file_read_batch The number of sequence reads to place in each batch (defaults to "25,000,000")
    -pgz/--pgzip_fastq Gzip fastq files (defaults to "0" (=no))
    -pfqc/--pfastqc Sequence quality analysis using FastQC (defaults to "0" (=no))
    -pcta/--pcutadapt trim input reads using cutadapt (defaults to "0" (=no))

    ##BWA
    -pmem/--pbwa_mem Align reads using Bwa Mem (defaults to "0" (=no))
      -memhla/--bwa_mem_hla Apply HLA typing (defaults to "1" (=yes))
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

    ##Bcftools
    -pbmp/--pbcftools_mpileup Variant calling using bcftools mpileup (defaults to "0" (=no))
      -pbmpfv/--bcftools_mpileup_filter_variant (Supply flag to enable)

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
      -gcgpss/--gatk_calculategenotypeposteriors_support_set GATK CalculateGenotypePosteriors support set (defaults to "1000g_sites_GRCh37_phase3_v4_20130502.vcf")
    -pgcv/--pgatk_combinevariantcallsets Combine variant call sets (defaults to "0" (=no))
      -gcvbcf/--gatk_combinevariantcallsets_bcf_file Produce a bcf from the GATK CombineVariantCallSet vcf (defaults to "1" (=yes))
      -gcvgmo/--gatk_combinevariants_genotype_merge_option Type of merge to perform (defaults to "PRIORITIZE")
      -gcvpc/--gatk_combinevariants_prioritize_caller The prioritization order of variant callers.(defaults to ""; comma sep; Options: gatk|bcftools|freebayes)
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
      -fqfgmf/--frequency_genmod_filter Remove common variants from vcf file (defaults to "1" (=yes))
      -fqfgfr/--frequency_genmod_filter_1000g Genmod annotate 1000G reference (defaults to "GRCh37_all_wgs_-phase3_v5b.2013-05-02-.vcf.gz")
      -fqfmaf/--frequency_genmod_filter_max_af Annotate MAX_AF from reference (defaults to "")
      -fqfgft/--frequency_genmod_filter_threshold Threshold for filtering variants (defaults to "0.10")
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
    -pars/--panalysisrunstatus Sets the analysis run status flag to finished in sample_info_file (defaults to "0" (=no))
    -psac/--psacct Generating sbatch script for SLURM info on each submitted job (defaults to "0" (=no))
    -sacfrf/--sacct_format_fields Format and fields of sacct output (defaults to "jobid", "jobname%50", "account", "partition", "alloccpus", "TotalCPU", "elapsed", "start", "end", "state", "exitcode")
END_USAGE
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

            my ( $volume, $directory, $fastq_file ) = splitpath($file);
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

sub set_custom_default_to_active_parameter {

    ## Function : Checks and sets user input or default values to active_parameters.
## Returns  :
## Arguments: $parameter_href        => Holds all parameters {REF}
##          : $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $file_info_href         => File info hash {REF}
##          : $parameter_name        => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $file_info_href;
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
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        parameter_name =>
          { required => 1, defined => 1, store => \$parameter_name, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::Parameter qw{ get_capture_kit };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## If capture kit is not set after cmd, config and reading pedigree
    if ( $parameter_name eq q{exome_target_bed} ) {

        ## Return a default capture kit as user supplied no info
        my $capture_kit = get_capture_kit(
            {
                supported_capture_kit_href =>
                  $parameter_href->{supported_capture_kit},
                capture_kit => q{latest},
            }
        );

        ## Set default
        $active_parameter_href->{exome_target_bed}
          {$capture_kit} = join q{,},
          @{ $active_parameter_href->{sample_ids} };

        $log->warn(
q{Could not detect a supplied capture kit. Will Try to use 'latest' capture kit: }
              . $capture_kit );
        return;
    }
    if ( $parameter_name eq q{bwa_build_reference} ) {

        ## Now we now what human genome reference to build from
        $active_parameter_href->{$parameter_name} =
          $active_parameter_href->{human_genome_reference};

        return;
    }
    ## Build default for analysis_type
    if ( $parameter_name eq q{analysis_type} ) {

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {
            $active_parameter_href->{$parameter_name}{$sample_id} = q{wgs};
        }
        return;
    }
    if ( $parameter_name eq q{infile_dirs} ) {
        ## Build default for infile_dirs

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $path = catfile(
                $active_parameter_href->{cluster_constant_path},
                $active_parameter_href->{family_id},
                $active_parameter_href->{analysis_type}{$sample_id},
                $sample_id,
                q{fastq}
            );

            $active_parameter_href->{$parameter_name}{$path} = $sample_id;
        }
        return;
    }
    if ( $parameter_name eq q{sample_info_file} ) {

        $parameter_href->{sample_info_file}{default} = catfile(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
            $active_parameter_href->{family_id} . q{_qc_sample_info.yaml}
        );

        $parameter_href->{qccollect_sampleinfo_file}{default} =
          $parameter_href->{sample_info_file}{default};
        return;
    }
    return;
}

sub set_default_to_active_parameter {

## Function : Checks and sets user input or default values to active_parameters.
## Returns  :
## Arguments: $parameter_href        => Holds all parameters
##          : $active_parameter_href => Holds all set parameter for analysis
##          : $parameter_name        => Parameter name
##          : $associated_programs   => The parameters program(s) {array, REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $associated_programs_ref;
    my $parameter_name;

    ## Default(s)
    my $family_id;

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
        associated_programs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$associated_programs_ref,
        },
        parameter_name =>
          { required => 1, defined => 1, store => \$parameter_name, },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my %only_wgs = ( gatk_genotypegvcfs_ref_gvcf => 1, );

    ## Alias
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};

    ## Do nothing since parameter is not required unless exome mode is enabled
    return
      if ( exists $only_wgs{$parameter_name}
        && $consensus_analysis_type =~ / wgs /xsm );

    ## Check all programs that use parameter
  ASSOCIATED_PROGRAM:
    foreach my $associated_program ( @{$associated_programs_ref} ) {

        ## Only add active programs parameters
        next ASSOCIATED_PROGRAM
          if ( not defined $active_parameter_href->{$associated_program} );

        next ASSOCIATED_PROGRAM
          if ( not $active_parameter_href->{$associated_program} );

        if ( exists $parameter_href->{$parameter_name}{default} ) {
            ## Default exists

            ## Array reference
            if ( $parameter_href->{$parameter_name}{data_type} eq q{ARRAY} ) {

                push
                  @{ $active_parameter_href->{$parameter_name} },
                  @{ $parameter_href->{$parameter_name}{default} };
            }
            elsif ( $parameter_href->{$parameter_name}{data_type} eq q{HASH} ) {
                ## Hash reference

                $active_parameter_href->{$parameter_name} =
                  $parameter_href->{$parameter_name}{default};
            }
            else {
                ## Scalar

                $active_parameter_href->{$parameter_name} =
                  $parameter_href->{$parameter_name}{default};
            }
            ## Set default - no use in continuing
            return;
        }
        else {
            ## No default

            ## Not mandatory - skip
            return
              if ( exists $parameter_href->{$parameter_name}{mandatory}
                && $parameter_href->{$parameter_name}{mandatory} eq q{no} );

            ## We have a logg object and somewhere to write
            $log->fatal($USAGE);
            $log->fatal( q{Supply '-}
                  . $parameter_name
                  . q{' if you want to run }
                  . $associated_program );
            exit 1;
        }
    }
    return;
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
                if ( any { $_ eq $order_parameter_element } @nowrite ) {
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

    foreach my $sample_id ( @{$sample_ids_ref} ) {

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

        $log->fatal(
q{MIP cannot detect what version of human_genome_reference you have supplied. Please supply the reference on this format: [sourceversion]_[species] e.g. 'GRCh37_homo_sapiens' or 'hg19_homo_sapiens'},
            "\n"
        );
        exit 1;
    }

    ## Removes ".file_ending" in filename.FILENDING(.gz)
    $file_info_href->{human_genome_reference_name_prefix} =
      fileparse( $$human_genome_reference_ref, qr/\.fasta|\.fasta\.gz/ );

    $file_info_href->{human_genome_compressed} =
      check_gzipped( { file_name_ref => $human_genome_reference_ref, } );
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
                  ; #No user supplied cmd info, not defined in config file, ADD vit from pedigree file
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

sub size_sort_select_file_contigs {

## Function : Sorts array depending on reference array. NOTE: Only entries present in reference array will survive in sorted array.
## Returns  : @sorted_contigs
## Arguments: $consensus_analysis_type_ref => Consensus analysis_type {REF}
##          : $file_info_href              => File info hash {REF}
##          : $hash_key_sort_reference     => The hash keys sort reference
##          : $hash_key_to_sort            => The keys to sort

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $consensus_analysis_type_ref;
    my $file_info_href;
    my $hash_key_sort_reference;
    my $hash_key_to_sort;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        consensus_analysis_type_ref => {
            default     => \$$,
            defined     => 1,
            required    => 1,
            store       => \$consensus_analysis_type_ref,
            strict_type => 1,
        },
        hash_key_to_sort => {
            defined     => 1,
            required    => 1,
            store       => \$hash_key_to_sort,
            strict_type => 1,
        },
        hash_key_sort_reference => {
            defined     => 1,
            required    => 1,
            store       => \$hash_key_sort_reference,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Hash qw{ check_element_exist_hash_of_array };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my @sorted_contigs;

    ## Sort the contigs depending on reference array
    if ( $file_info_href->{$hash_key_to_sort} ) {

        foreach my $element ( @{ $file_info_href->{$hash_key_sort_reference} } )
        {

            if (
                not check_element_exist_hash_of_array(
                    {
                        element  => $element,
                        hash_ref => $file_info_href,
                        key      => $hash_key_to_sort,
                    }
                )
              )
            {

                push @sorted_contigs, $element;
            }
        }
    }

    ## Test if all contigs collected from select file was sorted by reference contig array
    if ( @sorted_contigs
        && scalar @{ $file_info_href->{$hash_key_to_sort} } !=
        scalar @sorted_contigs )
    {

        foreach my $element ( @{ $file_info_href->{$hash_key_to_sort} } ) {

            ## If element is not part of array
            if ( not any { $_ eq $element } @sorted_contigs ) {

                ## Special case when analysing wes since Mitochondrial contigs have no baits in exome capture kits
                unless ( $$consensus_analysis_type_ref eq q{wes}
                    && $element =~ /MT$|M$/ )
                {

                    $log->fatal( q{Could not detect '##contig'= }
                          . $element
                          . q{ from meta data header in '-vcfparser_select_file' in reference contigs collected from '-human_genome_reference'}
                    );
                    exit 1;
                }
            }
        }
    }
    return @sorted_contigs;
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

    my $found_male        = 0;
    my $found_female      = 0;
    my $found_other       = 0;
    my $found_other_count = 0;

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
            $found_other_count++;
        }
    }
    return $found_male, $found_female, $found_other, $found_other_count;
}

sub remove_pedigree_elements {

## Function : Removes ALL keys at third level except keys in allowed_entries hash.
## Returns  :
## Arguments: $hash_ref => Hash {REF}

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

    my @allowed_entries =
      qw{ family default_gene_panels sample sample_id sample_name capture_kit sex mother father phenotype sequence_type expected_coverage sample_origin };

  FAMILY_INFO:
    for my $key ( keys %{$hash_ref} ) {

        ## If element is not part of array
        if ( not any { $_ eq $key } @allowed_entries ) {

            delete( $hash_ref->{$key} );
        }
    }

  SAMPLE:
    foreach my $sample_id ( keys %{ $hash_ref->{sample} } ) {

      SAMPLE_INFO:
        for my $pedigree_element ( keys %{ $hash_ref->{sample}{$sample_id} } ) {

            ## If element is not part of array
            if ( not any { $_ eq $pedigree_element } @allowed_entries ) {

                delete( $hash_ref->{sample}{$sample_id}{$pedigree_element} );
            }
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
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };

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

    use MIP::QC::Record qw(add_program_outfile_to_sample_info);

    if ( exists( $active_parameter_href->{analysis_type} ) ) {

        $sample_info_href->{analysis_type} =
          $active_parameter_href->{analysis_type};
    }
    if ( exists( $active_parameter_href->{expected_coverage} ) ) {

        $sample_info_href->{expected_coverage} =
          $active_parameter_href->{expected_coverage};
    }
    if ( exists $active_parameter_href->{sample_origin} ) {

        $sample_info_href->{sample_origin} =
          $active_parameter_href->{sample_origin};
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

        if (   ( defined $active_parameter_href->{$program} )
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

            if ( not defined $$outaligner_dir_ref ) {

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

    if ( defined $$human_genome_reference_source_ref ) {

      EXOME_FILE:
        foreach
          my $exome_target_bed_file ( keys %{$exome_target_bed_file_href} )
        {

            my $original_file_name = $exome_target_bed_file;

            ## Replace with actual version
            if ( $exome_target_bed_file =~
                s/genome_reference_source/$$human_genome_reference_source_ref/
                && $exome_target_bed_file =~
                s/_version/$$human_genome_reference_version_ref/ )
            {

                ## The delete operator returns the value being deleted i.e. updating hash key while preserving original info
                $exome_target_bed_file_href->{$exome_target_bed_file} =
                  delete $exome_target_bed_file_href->{$original_file_name};
            }
        }
    }
    return;
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

## Function : Check sample_id provided in hash parameter is included in the analysis and only represented once
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_href        => Holds all parameters {REF}
##          : $parameter_names_ref   => Parameter name list {REF}
##          : $sample_ids_ref        => Array to loop in for parameter {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_names_ref;
    my $parameter_href;
    my $sample_ids_ref;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        parameter_names_ref => {
            required => 1,
            defined  => 1,
            default  => [],
            store    => \$parameter_names_ref
        },
        sample_ids_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$sample_ids_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Loop through all hash parameters supplied
  PARAMETER:
    foreach my $parameter_name ( @{$parameter_names_ref} ) {

        if ( defined $active_parameter_href->{$parameter_name} ) {

            foreach my $sample_id ( @{$sample_ids_ref} ) {

                ## Check that a value exists
                if (
                    !defined(
                        $active_parameter_href->{$parameter_name}{$sample_id}
                    )
                  )
                {

                    next PARAMETER
                      if ( $parameter_href->{$parameter_name}{mandatory} eq
                        q{no} );

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

##Function : Test if key from query hash exists truth hash
##Returns  :
##Arguments: $parameter_name => Parameter name
##         : $truth_href     => Truth hash
##         : $query_href     => Query hash

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_name;
    my $truth_href;
    my $query_href;

    my $tmpl = {
        parameter_name =>
          { required => 1, defined => 1, store => \$parameter_name },
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

  QUERY_KEY:
    foreach my $key ( keys %{$query_href} ) {

        if ( not exists( $truth_href->{$key} ) ) {

            $log->fatal( $parameter_name
                  . q{ key '}
                  . $key
                  . q{' - Does not exist as module program parameter in MIP} );
            exit 1;
        }
    }
    return;
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

    foreach my $element ( @{$queryies} ) {

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
