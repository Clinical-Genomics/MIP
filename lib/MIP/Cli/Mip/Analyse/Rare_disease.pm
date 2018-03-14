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
    my $define_parameters_path =
      catfile( $Bin, qw{ definitions rare_disease_parameters.yaml } );

    ## Non mandatory parameter definition keys to check
    my $non_mandatory_parameter_keys_path = catfile( $Bin,
        qw{ definitions rare_disease_non_mandatory_parameter_keys.yaml } );

    ## Mandatory parameter definition keys to check
    my $mandatory_parameter_keys_path = catfile( $Bin,
        qw{ definitions rare_disease_mandatory_parameter_keys.yaml } );

    ### %parameter holds all defined parameters for MIP
    ### analyse rare_disease
    my %parameter = parse_definition_file(
        {
            define_parameters_path => $define_parameters_path,
            non_mandatory_parameter_keys_path =>
              $non_mandatory_parameter_keys_path,
            mandatory_parameter_keys_path => $mandatory_parameter_keys_path,
        }
    );

    ## Print programs and exit
    if ( $active_parameter{print_programs} ) {

        print_program(
            {
                define_parameters_file => $define_parameters_path,
                parameter_href         => \%parameter,
                print_program_mode     => $active_parameter{print_program_mode},
            }
        );
        exit;
    }

    ### To add/write parameters in the correct order
    ## Adds the order of first level keys from yaml file to array
    my @order_parameters = order_parameter_names(
        {
            file_path => $define_parameters_path,
        }
    );

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

    # write_args(\%parameter);
    #write_args( \%active_parameter );

    #exit;

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
            documentation =>
              q{Set the references to be decomposed and normalized},
            is  => q{rw},
            isa => q{ArrayRef},
        )
    );

    option(
        q{genomic_set} => (
            cmd_aliases   => [qw{ ges }],
            cmd_tags      => [q{sorted BED}],
            documentation => q{Selection of relevant regions post alignment},
            is            => q{ro},
            isa           => q{Str},
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
            isa           => q{Bool},
        )
    );

    option(
        q{bwa_mem_cram} => (
            cmd_aliases   => [qw{ memcrm }],
            documentation => q{Use CRAM-format for additional output file},
            is            => q{rw},
            isa           => q{Bool},
        )
    );

    option(
        q{bwa_mem_bamstats} => (
            cmd_aliases   => [qw{ memsts }],
            documentation => q{Collect statistics from BAM files},
            is            => q{rw},
            isa           => q{Bool},
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
            isa => q{Str},
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
            isa           => q{Bool},
        )
    );

    option(
        q{markduplicates_sambamba_markdup} => (
            cmd_aliases   => [qw{ mdsmd }],
            cmd_flag      => q{sambamba_markdup},
            documentation => q{Markduplicates using Sambamba markduplicates},
            is            => q{rw},
            isa           => q{Bool},
        )
    );

    option(
        q{markduplicates_sambamba_markdup} => (
            cmd_aliases   => [qw{ mdsmd }],
            cmd_flag      => q{sambamba_markdup},
            documentation => q{Markduplicates using Sambamba markduplicates},
            is            => q{rw},
            isa           => q{Bool},
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
            isa => q{Int},
        )
    );

    option(
        q{markduplicates_sambamba_markdup_overflow_list_size} => (
            cmd_aliases   => [qw{ mdsols }],
            cmd_flag      => q{sba_mdup_ols},
            cmd_tags      => [q{Default: 200000}],
            documentation => q{Sambamba size of the overflow list},
            is            => q{rw},
            isa           => q{Int},
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
            isa => q{Int},
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
            isa           => q{ArrayRef[Int]},
        )
    );

    option(
        q{sambamba_depth_bed} => (
            cmd_aliases   => [qw{ sdtbed }],
            documentation => q{Reference bed file},
            is            => q{rw},
            isa           => q{Str},
        )
    );

    option(
        q{sambamba_depth_base_quality} => (
            cmd_aliases   => [qw{ sdtbaq }],
            cmd_flag      => q{sba_depth_bq},
            cmd_tags      => [q{Default: 10}],
            documentation => q{Do not count bases with lower base quality},
            is            => q{rw},
            isa           => q{Int},
        )
    );

    option(
        q{sambamba_depth_mapping_quality} => (
            cmd_aliases   => [qw{ sdtmaq }],
            cmd_flag      => q{sba_depth_mq},
            cmd_tags      => [q{Default: 10}],
            documentation => q{Do not count reads with lower mapping quality},
            is            => q{rw},
            isa           => q{Int},
        )
    );

    option(
        q{sambamba_depth_noduplicates} => (
            cmd_aliases => [qw{ sdtndu }],
            cmd_flag    => q{sba_depth_nod},
            documentation =>
              q{Do not include duplicates in coverage calculation},
            is  => q{rw},
            isa => q{Bool},
        )
    );

    option(
        q{sambamba_depth_quality_control} => (
            cmd_aliases => [qw{ sdtfqc }],
            cmd_flag    => q{sba_depth_qc},
            documentation =>
              q{Do not include reads with failed quality control},
            is  => q{rw},
            isa => q{Bool},
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
            isa           => q{Int},
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
            isa           => q{Int},
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
            isa => q{Str},
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
            isa           => q{Int},
        )
    );

    option(
        q{tiddit_bin_size} => (
            cmd_aliases   => [qw{ tidbin }],
            cmd_tags      => [q{Default: 500}],
            documentation => q{Size of coverage bins in calculation},
            is            => q{rw},
            isa           => q{Int},
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
            isa           => q{Bool},
        )
    );

    option(
        q{sv_svdb_merge_prioritize} => (
            cmd_aliases => [qw{ svsvdbmp }],
            documentation =>
              q{Prioritization order of structural variant callers},
            is  => q{rw},
            isa => q{Str},
        )
    );

    option(
        q{sv_bcftools_view_filter} => (
            cmd_aliases => [qw{ svcbtv }],
            documentation =>
              q{Include structural variants with PASS in FILTER column},
            is  => q{rw},
            isa => q{Bool},
        )
    );

    option(
        q{sv_svdb_query} => (
            cmd_aliases   => [qw{ svcdbq }],
            documentation => q{Annotate structural variants using svdb query},
            is            => q{rw},
            isa           => q{Bool},
        )
    );

    option(
        q{sv_svdb_query_db_files} => (
            cmd_aliases   => [qw{ svcdbqd }],
            cmd_tags      => [q{file.vcf=vcf_info_key}],
            documentation => q{Database file(s) for annotation},
            is            => q{rw},
            isa           => q{HashRef},
        )
    );

    option(
        q{sv_vcfanno} => (
            cmd_aliases   => [qw{ svcvan }],
            documentation => q{Annotate structural variants},
            is            => q{rw},
            isa           => q{Bool},
        )
    );

    option(
        q{sv_vcfanno_lua} => (
            cmd_aliases   => [qw{ svcval }],
            documentation => q{VcfAnno lua postscripting file},
            is            => q{rw},
            isa           => q{Str},
        )
    );

    option(
        q{sv_vcfanno_config} => (
            cmd_aliases   => [qw{ svcvac }],
            documentation => q{VcfAnno toml config},
            is            => q{rw},
            isa           => q{Str},
        )
    );

    option(
        q{sv_vcfanno_config_file} => (
            cmd_aliases => [qw{ svcvacf }],
            cmd_tags =>
              [q{Default: GRCh37_all_sv_-phase3_v2.2013-05-02-.vcf.gz}],
            documentation => q{Annotation file within vcfAnno config toml file},
            is            => q{rw},
            isa           => q{Str},
        )
    );

    option(
        q{sv_vcfannotation_header_lines_file} => (
            cmd_aliases => [qw{ svcvah }],
            cmd_flag    => q{sv_vcfanno_hlf},
            documentation =>
              q{Adjust for postscript by adding required header lines to vcf},
            is  => q{rw},
            isa => q{Str},
        )
    );

    option(
        q{sv_genmod_filter} => (
            cmd_aliases   => [qw{ svcgmf }],
            documentation => q{Remove common structural variants from vcf},
            is            => q{rw},
            isa           => q{Bool},
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
            isa => q{Str},
        )
    );

    option(
        q{sv_genmod_filter_threshold} => (
            cmd_aliases   => [qw{ svcgft }],
            cmd_tags      => [q{Default: 0.10}],
            documentation => q{Threshold for filtering structural variants},
            is            => q{rw},
            isa           => q{Num},
        )
    );

    option(
        q{sv_combinevariantcallsets_bcf_file} => (
            cmd_aliases => [qw{ svcbcf }],
            documentation =>
              q{Produce a bcf from the CombineStructuralVariantCallSet vcf},
            is  => q{rw},
            isa => q{Bool},
        )
    );

    return;
}

sub write_args {

    my ($arg_href) = @_;

    # do something
    use Data::Dumper;

    #    say STDERR $arg_href->{pbwa_mem};
    #    foreach my $sample ( @{ $arg_href->{sample_ids} } ) {
    #        say STDERR $sample;
    #    }
    print Dumper($arg_href);
    return;
}

1;
