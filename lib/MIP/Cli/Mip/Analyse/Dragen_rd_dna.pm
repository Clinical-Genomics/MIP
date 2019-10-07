package MIP::Cli::Mip::Analyse::Dragen_rd_dna;

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
use MooseX::Types::Moose qw{ Str Int HashRef Num Bool ArrayRef };
use Moose::Util::TypeConstraints;

## MIPs lib
use MIP::Main::Analyse qw{ mip_analyse };

our $VERSION = 1.05;

extends(qw{ MIP::Cli::Mip::Analyse });

command_short_description(q{Dragen rare disease DNA analysis});

command_long_description(q{Dragen rare disease DNA analysis on wgs sequence data});

command_usage(q{mip <analyse> <dragen_rd_dna> <case_id> --config <config_file> });

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

    ## Mip analyse rare_disease parameters
    ## CLI commands inheritance
    my @definition_files = (
        catfile( $Bin, qw{ definitions mip_parameters.yaml } ),
        catfile( $Bin, qw{ definitions analyse_parameters.yaml } ),
        catfile( $Bin, qw{ definitions dragen_rd_dna_parameters.yaml } ),
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

    ## Get dependency tree and store in parameter hash
    my %dependency_tree = load_yaml(
        {
            yaml_file =>
              catfile( $Bin, qw{ definitions dragen_rd_dna_initiation_map.yaml } ),
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

        # Human genome meta files
        human_genome_reference_file_endings => [qw{ .dict .fai }],
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
        q{human_genome_reference} => (
            cmd_aliases   => [qw{ hgr }],
            cmd_tags      => [q{Default: grch37_homo_sapiens_-d5-.fasta}],
            documentation => q{Human genome reference},
            is            => q{rw},
            isa           => Str,
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
        q{recipe_core_number} => (
            cmd_aliases   => [qw{ rcn }],
            cmd_tags      => [q{recipe_name=X(cores)}],
            documentation => q{Set the number of cores for each recipe},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{recipe_memory} => (
            cmd_aliases   => [qw{ rm }],
            cmd_tags      => [q{recipe_name=X(G)}],
            documentation => q{Set the memory for each recipe},
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
        q{dragen_dna_align_vc} => (
            cmd_aliases   => [qw{ drgdav }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Dragen align and variant call recipe},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{dragen_analysis_dir} => (
            cmd_aliases   => [qw{ drgad }],
            cmd_flag      => q{dragen_analysis_dir},
            documentation => q{Dragen analysis dir},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{dragen_fastq_list_file_path} => (
            cmd_aliases   => [qw{ drgflf }],
            cmd_flag      => q{dragen_fastq_list_file_path},
            cmd_tags      => [q{Format: csv}],
            documentation => q{Dragen fastq list sample id file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{dragen_hash_ref_dir_path} => (
            cmd_aliases   => [qw{ drghrdp }],
            cmd_flag      => q{dragen_hash_ref_dir_path},
            documentation => q{Dragen hash reference dir path},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{dragen_user_at_hostname} => (
            cmd_aliases   => [qw{ drguah }],
            cmd_flag      => q{dragen_user_at_hostname},
            documentation => q{Dragen user at hostname},
            is            => q{rw},
            isa           => Str,
        )
    );
    option(
        q{dragen_joint_calling} => (
            cmd_aliases   => [qw{ drgdjc }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Dragen joint call recipe},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
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
            cmd_aliases   => [qw{ fqacvac }],
            documentation => q{Frequency vcfanno toml config},
            is            => q{rw},
            isa           => Str,
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
        q{vep_directory_path} => (
            cmd_aliases   => [qw{ vepp }],
            documentation => q{Path to VEP script directory},
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
        q{samtools_subsample_mt} => (
            cmd_aliases   => [qw{ ssmt }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Subsample the mitochondria reads},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    return;
}

1;
