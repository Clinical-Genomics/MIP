package MIP::Cli::Mip::Analyse::Cancer;

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

command_short_description(q{Cancer analysis});

command_long_description(q{Cancer analysis on panel, wes or wgs sequence data});

command_usage(q{mip <analyse> <cancer> <family_id> --config <config_file>});

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

    ## Mip analyse cancer parameters
    my @definition_files = (
        catfile( $Bin, qw{ definitions mip_parameters.yaml } ),
        catfile( $Bin, qw{ definitions analyse_parameters.yaml } ),
        catfile( $Bin, qw{ definitions cancer_parameters.yaml } ),
    );

    ## Non mandatory parameter definition keys to check
    my $non_mandatory_parameter_keys_path =
      catfile( $Bin, qw{ definitions non_mandatory_parameter_keys.yaml } );

    ## Mandatory parameter definition keys to check
    my $mandatory_parameter_keys_path =
      catfile( $Bin, qw{ definitions mandatory_parameter_keys.yaml } );

    ### %parameter holds all defined parameters for MIP
    ### mip analyse cancer
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
    ## Adds the order of first level keys from yaml file to array
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
            cmd_tags    => [q{gatk_baserecalibration_known_sites}],
            documentation =>
              q{Set the references to be decomposed and normalized},
            is  => q{rw},
            isa => ArrayRef [Str],
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
            cmd_aliases => [qw{ gdai }],
            cmd_flag    => q{gatk_dis_auto_ind_fl},
            documentation =>
              q{Disable auto index creation and locking when reading rods},
            is  => q{rw},
            isa => Bool,
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
        q{gatk_logging_level} => (
            cmd_aliases   => [qw{ gll }],
            cmd_tags      => [q{Default: INFO}],
            documentation => q{Set the GATK log level},
            is            => q{rw},
            isa           => enum( [qw{ DEBUG INFO ERROR FATAL }] ),
        )
    );

    option(
        q{gatk_path} => (
            cmd_aliases   => [qw{ gtp }],
            documentation => q{Path to GATK},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{human_genome_reference} => (
            cmd_aliases   => [qw{ hgr }],
            cmd_tags      => [q{Default: GRCh37_homo_sapiens_-d5-.fasta}],
            documentation => q{Human genome reference},
            is            => q{rw},
            isa           => Str,
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
            isa => enum( [qw{ bwa }] ),
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
            cmd_aliases   => [qw{ sao }],
            cmd_tags      => [q{sample_id=sample_origin}],
            documentation => q{Sample origin for analysis},
            is            => q{rw},
            isa           => HashRef,
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
        q{psambamba_depth} => (
            cmd_aliases   => [qw{ psdt }],
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
        q{pvardict} => (
            cmd_aliases   => [qw{ pvrd }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Variant calling using Vardict},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{vrd_af_threshold} => (
            cmd_aliases   => [qw{ vdraf }],
            cmd_tags      => [q{Default: 0.01}],
            documentation => q{AF threshold for variant calling},
            is            => q{rw},
            isa           => Num,
        )
    );

    option(
        q{vrd_chrom_start} => (
            cmd_aliases   => [qw{ vrdcs }],
            cmd_tags      => [q{Default: 1}],
            documentation => q{Column for chromosome in the output},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{vrd_input_bed_file} => (
            cmd_aliases   => [qw{ vrdbed }],
            documentation => q{Infile path for region info bed file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{vrd_max_mm} => (
            cmd_aliases   => [qw{ vrdmm }],
            cmd_tags      => [q{Default: 4.5}],
            documentation => q{Maximum mean mismatches allowed},
            is            => q{rw},
            isa           => Num,
        )
    );

    option(
        q{vrd_max_pval} => (
            cmd_aliases   => [qw{ vrdmp }],
            cmd_tags      => [q{Default: 0.9}],
            documentation => q{Maximum p-value, set to 0 to keep all variants},
            is            => q{rw},
            isa           => Num,
        )
    );

    option(
        q{vrd_region_end} => (
            cmd_aliases   => [qw{ vrdre }],
            cmd_tags      => [q{Default: 3}],
            documentation => q{Column for region end position in the output},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{vrd_region_start} => (
            cmd_aliases   => [qw{ vrdrs }],
            cmd_tags      => [q{Default: 2}],
            documentation => q{Column for region start position in the output},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{vrd_segment_annotn} => (
            cmd_aliases   => [qw{ vrdsa }],
            cmd_tags      => [q{Default: 4}],
            documentation => q{Column for segment annotation in the output},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{vrd_somatic_only} => (
            cmd_aliases   => [qw{ vrdso }],
            cmd_tags      => [q{Default: no}],
            documentation => q{Output only candidate somatic},
            is            => q{rw},
            isa           => Str,
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

    return;
}

1;
