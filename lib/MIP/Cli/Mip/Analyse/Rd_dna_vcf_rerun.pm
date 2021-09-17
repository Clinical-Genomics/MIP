package MIP::Cli::Mip::Analyse::Rd_dna_vcf_rerun;

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

command_short_description(q{Rare disease DNA vcf rerun analysis});

command_long_description(q{Rare disease DNA vcf rerun analysis on wes, wgs or mixed sequence data});

command_usage(q{mip <analyse> <rd_dna_vcf_rerun> <case_id> --config <config_file> });

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

    ## %parameter holds all defined parameters for MIP analyse rd_dna_vcf_rerun
    ## CLI commands inheritance
    my $level     = q{rd_dna_vcf_rerun};
    my %parameter = get_parameter_from_definition_files( { level => $level, } );

    my @rd_dna_vcf_rerun_definition_file_paths =
      get_parameter_definition_file_paths( { level => $level, } );

    ### To write parameters and their values to log in logical order
    ## Adds the order of first level keys from definition files to array
    my @order_parameters = get_order_of_parameters(
        { define_parameters_files_ref => \@rd_dna_vcf_rerun_definition_file_paths, } );

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
            cmd_flag => q{dec_norm_ref},
            cmd_tags => [
q{gatk_baserecalibration_known_sites, gatk_haplotypecaller_snp_known_set, gatk_variantrecalibration_resource_snv, gatk_variantrecalibration_resource_indel, frequency_genmod_filter_1000g, gatk_varianteval_gold, gatk_varianteval_dbsnp}
            ],
            documentation => q{Set the references to be decomposed and normalized},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{gatk_path} => (
            documentation => q{Path to GATK directory},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{human_genome_reference} => (
            cmd_tags      => [q{Default: grch37_homo_sapiens_-d5-.fasta}],
            documentation => q{Human genome reference},
            is            => q{rw},
            isa           => Str,
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
        q{picardtools_path} => (
            cmd_aliases => [qw{ ptp }],
            is          => q{rw},
            isa         => Str,
        )
    );

    option(
        q{sv_vcf_rerun_reformat} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Reformat rerun SV variants},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{sv_vcf_rerun_file} => (
            cmd_flag      => q{sv_vcf_rerun_file},
            cmd_tags      => [q{Format: vcf | bcf}],
            documentation => q{Sv variant calling file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sv_annotate} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Annotate and filter structural variant calls},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{sv_vcfanno_config} => (
            documentation => q{Structural variants vcfanno toml config},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sv_fqa_annotations} => (
            documentation => q{Frequency annotations to use},
            is            => q{rw},
            isa           => ArrayRef,
        )
    );

    option(
        q{sv_fqa_filter} => (
            documentation => q{Frequency annotations to use when filtering},
            is            => q{rw},
            isa           => ArrayRef,
        )
    );

    option(
        q{sv_frequency_filter} => (
            documentation => q{Remove common structural variants from vcf},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{sv_svdb_query} => (
            documentation => q{Annotate structural variants using svdb query},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{sv_svdb_query_db_files} => (
            cmd_tags      => [q{file.vcf=vcf_info_key}],
            documentation => q{Database file(s) for annotation},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{sv_varianteffectpredictor} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Annotate SV variants using VEP},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{sv_vep_features} => (
            cmd_tags => [
q{Default: hgvs, symbol, numbers, sift, polyphen, humdiv, domains, protein, ccds, uniprot, biotype, regulatory, tsl, canonical, per_gene, appris}
            ],
            documentation => q{VEP features},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{sv_vcfparser} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Parse structural variants using vcfParser.pl},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{sv_vcfparser_add_all_mt_var} => (
            cmd_flag      => q{sv_vcfparser_all_mt},
            documentation => q{Add all MT variants in select vcf},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{sv_vcfparser_per_gene} => (
            documentation => q{Keep only most severe consequence per gene},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{sv_vcfparser_range_feature_annotation_columns} => (
            cmd_flag      => q{sv_vcfparser_fac},
            documentation => q{Range annotations feature columns},
            is            => q{rw},
            isa           => ArrayRef [Int],
        )
    );

    option(
        q{sv_vcfparser_range_feature_file} => (
            cmd_flag      => q{sv_vcfparser_rff},
            cmd_tags      => [q{Format: tsv}],
            documentation => q{Range annotations file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sv_vcfparser_select_feature_annotation_columns} => (
            cmd_flag      => q{sv_vcfparser_slt_fac},
            documentation => q{Feature columns to use in annotation},
            is            => q{rw},
            isa           => ArrayRef [Int],
        )
    );

    option(
        q{sv_vcfparser_select_file} => (
            cmd_flag      => q{sv_vcfparser_slt_fl},
            cmd_tags      => [q{Format: tsv; HGNC Symbol required in file}],
            documentation => q{Select file with list of genes to analyse separately},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sv_vcfparser_select_file_matching_column} => (
            cmd_flag      => q{sv_vcfparser_slt_fmc},
            documentation => q{Position of HGNC Symbol column in select file},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{sv_vcfparser_vep_transcripts} => (
            cmd_flag      => q{sv_vcfparser_vtr},
            documentation => q{Parse VEP transcript specific entries},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{sv_rankvariant} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Ranking of annotated SV variants},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{sv_genmod_annotate_regions} => (
            cmd_flag      => q{sv_genmod_ann_reg},
            documentation =>
              q{Use predefined gene annotation supplied with genmod for defining genes},
            is  => q{rw},
            isa => Bool,
        )
    );

    option(
        q{sv_genmod_models_case_type} => (
            cmd_flag      => q{sv_genmod_mod_fam_typ},
            cmd_tags      => [q{Default: mip}],
            documentation => q{Use one of the known setups},
            is            => q{rw},
            isa           => enum( [qw{ped alt cmms mip}] ),
        )
    );

    option(
        q{sv_genmod_models_reduced_penetrance_file} => (
            cmd_flag      => q{sv_genmod_mod_red_pen_f},
            documentation => q{File containing genes with reduced penetrance},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sv_genmod_models_whole_gene} => (
            cmd_flag      => q{sv_genmod_mod_whl_gene},
            documentation => q{Allow compound pairs in intronic regions},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{sv_rank_model_file} => (
            documentation => q{Rank model config file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sv_reformat} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Concatenating files},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{sv_reformat_remove_genes_file} => (
            cmd_flag      => q{sv_reformat_rem_gen_f},
            documentation => q{Remove variants with hgnc_ids from file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{vcf_rerun_reformat} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Reformat rerun SV variants},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{vcf_rerun_file} => (
            cmd_flag      => q{vcf_rerun_file},
            cmd_tags      => [q{Format: vcf | bcf }],
            documentation => q{Variant calling file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{prepareforvariantannotationblock} => (
            cmd_flag      => q{prep_for_var_ann_bl},
            cmd_tags      => [q{Analysis recipe switch}],
            documentation =>
              q{Prepare for variant annotation block by copying and splitting files per contig},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{rhocall_ar} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Rhocall performs annotation of variants in autozygosity regions},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{rhocall_frequency_file} => (
            cmd_tags      => [q{Default: grch37_anon_swegen_snp_-2016-10-19-.tab.gz; tsv}],
            documentation => q{Frequency file for bcftools roh calculation},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{bcftools_norm} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Decompose and normalize},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{bcftools_missing_alt_allele} => (
            documentation => q{Remove missing alternative alleles '*'},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{bcftools_normalize} => (
            documentation => q{Normalize variants},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{variant_annotation} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Annotate vcf},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{vcfanno_config} => (
            documentation => q{SNV/Indel vcfanno toml config},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{frequency_filter} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Filter variants on frequency},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{fqf_annotations} => (
            documentation => q{Frequency annotations to use when filtering },
            is            => q{rw},
            isa           => ArrayRef,
        )
    );

    option(
        q{fqf_bcftools_filter_threshold} => (
            cmd_flag      => q{freq_bcftools_fil_trh},
            cmd_tags      => [q{Default: 0.10}],
            documentation => q{Threshold for filtering variants},
            is            => q{rw},
            isa           => Num,
        )
    );

    option(
        q{cadd_ar} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Annotate variants with CADD},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{cadd_column_names} => (
            documentation => q{Column names in cadd tsv},
            is            => q{rw},
            isa           => ArrayRef,
        )
    );

    option(
        q{cadd_vcf_header_file} => (
            documentation => q{},
            is            => q{rw},
            isa           => Str,
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

    option(
        q{vcfparser_ar} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Parse structural variants using vcfParser.pl},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{vcfparser_add_all_mt_var} => (
            cmd_flag      => q{vcfparser_all_mt},
            documentation => q{Add all MT variants in select vcf},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{vcfparser_pli_score_file} => (
            cmd_tags      => [q{Format: tsv}],
            documentation => q{Gene pLI score file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{vcfparser_range_feature_annotation_columns} => (
            cmd_flag      => q{vcfparser_fac},
            documentation => q{Range annotations feature columns},
            is            => q{rw},
            isa           => ArrayRef [Int],
        )
    );

    option(
        q{vcfparser_range_feature_file} => (
            cmd_flag      => q{vcfparser_rff},
            cmd_tags      => [q{Format: tsv}],
            documentation => q{Range annotations file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{vcfparser_select_file} => (
            cmd_flag      => q{vcfparser_slt_fl},
            cmd_tags      => [q{Format: tsv; HGNC Symbol required in file}],
            documentation => q{Select file with list of genes to analyse separately},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{vcfparser_select_feature_annotation_columns} => (
            cmd_flag      => q{vcfparser_slt_fac},
            documentation => q{Feature columns to use in annotation},
            is            => q{rw},
            isa           => ArrayRef [Int],
        )
    );

    option(
        q{vcfparser_select_file_matching_column} => (
            cmd_flag      => q{vcfparser_slt_fmc},
            documentation => q{Position of HGNC Symbol column in select file},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{vcfparser_vep_transcripts} => (
            cmd_flag      => q{vcfparser_vtr},
            documentation => q{Parse VEP transcript specific entries},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{rankvariant} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Ranking of annotated variants},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{genmod_annotate_regions} => (
            cmd_flag      => q{genmod_ann_reg},
            documentation =>
              q{Use predefined gene annotation supplied with genmod for defining genes},
            is  => q{rw},
            isa => Bool,
        )
    );

    option(
        q{genmod_models_case_type} => (
            cmd_flag      => q{genmod_mod_fam_typ},
            cmd_tags      => [q{Default: mip}],
            documentation => q{Use one of the known setups},
            is            => q{rw},
            isa           => enum( [qw{ped alt cmms mip}] ),
        )
    );

    option(
        q{genmod_models_reduced_penetrance_file} => (
            cmd_flag      => q{genmod_mod_red_pen_f},
            documentation => q{File containing genes with reduced penetrance},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{genmod_models_whole_gene} => (
            cmd_flag      => q{genmod_mod_whl_gene},
            documentation => q{Allow compound pairs in intronic regions},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{rank_model_file} => (
            documentation => q{Rank model config file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{endvariantannotationblock} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{End variant annotation block by concatenating files},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{endvariantannotationblock_remove_genes_file} => (
            cmd_flag      => q{endvarannbl_rem_gen_f},
            documentation => q{Remove variants with hgnc_ids from file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{samtools_subsample_mt} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Subsample the mitochondria reads},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    return;
}

1;
