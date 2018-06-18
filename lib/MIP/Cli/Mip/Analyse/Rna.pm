package MIP::Cli::Mip::Analyse::Rna;

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

our $VERSION = 0.03;

extends(qw{ MIP::Cli::Mip::Analyse });

command_short_description(q{Rna analysis});

command_long_description(q{Rna analysis on wts sequence data});

command_usage(q{mip <analyse> <rna> <family_id> --config <config_file> });

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
    use MIP::Get::Analysis qw{ get_dependency_tree_order print_program };

    ## Mip analyse rna parameters
    ## CLI commands inheritance
    my @definition_files = (
        catfile( $Bin, qw{ definitions mip_parameters.yaml } ),
        catfile( $Bin, qw{ definitions analyse_parameters.yaml } ),
        catfile( $Bin, qw{ definitions rna_parameters.yaml } ),
    );

    ## Non mandatory parameter definition keys to check
    my $non_mandatory_parameter_keys_path =
      catfile( $Bin, qw{ definitions non_mandatory_parameter_keys.yaml } );

    ## Mandatory parameter definition keys to check
    my $mandatory_parameter_keys_path =
      catfile( $Bin, qw{ definitions mandatory_parameter_keys.yaml } );

    ### %parameter holds all defined parameters for MIP
    ### mip analyse rna
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
            ),
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

    ## Order programs - Parsed from initiation file
    my %dependency_tree = load_yaml(
        {
            yaml_file => catfile( $Bin, qw{ definitions rna_initiation.yaml } ),
        }
    );

    my @order_programs;
    get_dependency_tree_order(
        {
            dependency_tree_href => \%dependency_tree,
            programs_ref         => \@order_programs
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
            order_programs_ref    => \@order_programs,
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
        q{java_use_large_pages} => (
            cmd_aliases   => [qw{ jul }],
            documentation => q{Use large page memory},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{module_core_number} => (
            cmd_aliases   => [qw{ mcn }],
            cmd_tags      => [q{program_name=X(cores)}],
            documentation => q{Set the number of cores for each module},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{module_time} => (
            cmd_aliases   => [qw{ mot }],
            cmd_tags      => [q{program_name=time(hours)}],
            documentation => q{Set the time allocation for each module},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{outaligner_dir} => (
            cmd_aliases => [qw{ ald }],
            documentation =>
q{Sets which aligner out directory was used for alignment in previous analysis},
            is  => q{rw},
            isa => enum( [qw{ star }] ),
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
        q{ppicardtools_mergesamfiles} => (
            cmd_aliases   => [qw{ ppms }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Merge bam files using Picardtools},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{psalmon_quant} => (
            cmd_aliases   => [qw{ psqt }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Quantify transcripts using salmon},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );
    option(
        q{salmon_rna_lib_configuration} => (
            cmd_aliases   => [qw{ psqt_bob }],
            cmd_tags      => [q{Default: ISF}],
            documentation => q{Library orientation and strandedness},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{pstar_aln} => (
            cmd_aliases   => [qw{ pstn }],
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
            cmd_aliases => [qw{ stn_asom }],
            cmd_tags    => [q{Default: 10}],
            documentation =>
              q{Minimum overhang (i.e. block size) for spliced alignments},
            is  => q{rw},
            isa => Int,
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
        q{pstar_fusion} => (
            cmd_aliases   => [qw{ pstf }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Detect fusion transcripts with star fusion},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
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
        q{pgatk_splitncigarrreads} => (
            cmd_aliases   => [qw{ pgs }],
            cmd_flag      => q{pgatk_splitncigarreads},
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Split reads that contain Ns in their cigar},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
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
        q{pgatk_asereadcounter} => (
            cmd_aliases   => [qw{ pgae }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Allel specific expression},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{pgatk_variantfiltration} => (
            cmd_aliases   => [qw{ pgvf }],
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

    return;
}

1;
