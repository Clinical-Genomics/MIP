package MIP::Recipes::Pipeline::Rare_disease;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Spec::Functions qw{ catdir catfile };

## CPANM
use Readonly;

##MIPs lib/
use MIP::Delete::List qw{ delete_male_contig };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ pipeline_rare_disease };
}

## Constants
Readonly my $SPACE => q{ };

sub pipeline_rare_disease {

## Function : Pipeline recipe for wes and or wgs data analysis.
## Returns  :

## Arguments: $parameter_href          => Parameter hash {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $indir_path_href         => Indirectory hash {REF}
##          : $infile_href             => Infile hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $lane_href               => The lane info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner dir used in the analysis
##          : $log                     => Log object to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $indir_path_href;
    my $infile_href;
    my $infile_lane_prefix_href;
    my $lane_href;
    my $job_id_href;
    my $log;

    ## Default(s)
    my $outaligner_dir;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        indir_path_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$indir_path_href,
            strict_type => 1,
        },
        infile_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_href,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        lane_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$lane_href,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Recipes
    use MIP::Recipes::Analysis::Analysisrunstatus
      qw{ analysis_analysisrunstatus };
    use MIP::Recipes::Analysis::Bamcalibrationblock
      qw{ analysis_bamcalibrationblock };
    use MIP::Recipes::Analysis::Bcftools_mpileup
      qw { analysis_bcftools_mpileup };
    use MIP::Recipes::Analysis::Bedtools_genomecov
      qw{ analysis_bedtools_genomecov };
    use MIP::Recipes::Analysis::Bwa_mem qw{ analysis_bwa_mem };
    use MIP::Recipes::Build::Human_genome_prerequisites
      qw{ build_human_genome_prerequisites };
    use MIP::Recipes::Analysis::Chanjo_sex_check
      qw{ analysis_chanjo_sex_check };
    use MIP::Recipes::Analysis::Cnvnator qw{ analysis_cnvnator };
    use MIP::Recipes::Analysis::Delly_call qw{ analysis_delly_call };
    use MIP::Recipes::Analysis::Delly_reformat qw{ analysis_delly_reformat };
    use MIP::Recipes::Analysis::Endvariantannotationblock
      use MIP::Recipes::Analysis::Expansionhunter
      qw{ analysis_expansionhunter };
    use MIP::Recipes::Analysis::Fastqc qw{ analysis_fastqc };
    use MIP::Recipes::Analysis::Freebayes qw { analysis_freebayes_calling };
    use MIP::Recipes::Analysis::Frequency_filter
      qw{ analysis_frequency_filter };
    use MIP::Recipes::Analysis::Gatk_baserecalibration
      qw{ analysis_gatk_baserecalibration };
    use MIP::Recipes::Analysis::Gatk_combinevariantcallsets
      qw{ analysis_gatk_combinevariantcallsets };
    use MIP::Recipes::Analysis::Gatk_concatenate_genotypegvcfs
      qw{ analysis_gatk_concatenate_genotypegvcfs };
    use MIP::Recipes::Analysis::Gatk_genotypegvcfs
      qw{ analysis_gatk_genotypegvcfs };
    use MIP::Recipes::Analysis::Gatk_haplotypecaller
      qw{ analysis_gatk_haplotypecaller };
    use MIP::Recipes::Analysis::Gatk_realigner qw{ analysis_gatk_realigner };
    use MIP::Recipes::Analysis::Gatk_variantevalall
      qw{ analysis_gatk_variantevalall };
    use MIP::Recipes::Analysis::Gatk_variantevalexome
      qw{ analysis_gatk_variantevalexome };
    use MIP::Recipes::Analysis::Gatk_variantrecalibration
      qw{ analysis_gatk_variantrecalibration_wgs analysis_gatk_variantrecalibration_wes };
    use MIP::Recipes::Analysis::Manta qw{ analysis_manta };
    use MIP::Recipes::Analysis::Markduplicates qw{ analysis_markduplicates };
    use MIP::Recipes::Analysis::Mip_vcfparser
      qw{ analysis_mip_vcfparser analysis_sv_vcfparser };
    use MIP::Recipes::Analysis::Multiqc qw{ analysis_multiqc };
    use MIP::Recipes::Analysis::Peddy qw{ analysis_peddy };
    use MIP::Recipes::Analysis::Picardtools_collecthsmetrics
      qw{ analysis_picardtools_collecthsmetrics };
    use MIP::Recipes::Analysis::Picardtools_collectmultiplemetrics
      qw{ analysis_picardtools_collectmultiplemetrics };
    use MIP::Recipes::Analysis::Picardtools_genotypeconcordance
      qw{ analysis_picardtools_genotypeconcordance };
    use MIP::Recipes::Analysis::Picardtools_mergesamfiles
      qw{ analysis_picardtools_mergesamfiles };
    use MIP::Recipes::Analysis::Plink qw{ analysis_plink };
    use MIP::Recipes::Analysis::Prepareforvariantannotationblock
      qw{ analysis_prepareforvariantannotationblock };
    use MIP::Recipes::Analysis::Qccollect qw{ analysis_qccollect };
    use MIP::Recipes::Analysis::Rankvariant
      qw{ analysis_rankvariant analysis_rankvariant_unaffected analysis_sv_rankvariant analysis_sv_rankvariant_unaffected };
    use MIP::Recipes::Analysis::Rcoverageplots qw{ analysis_rcoverageplots };
    use MIP::Recipes::Analysis::Rhocall qw{ analysis_rhocall_annotate };
    use MIP::Recipes::Analysis::Rtg_vcfeval qw{ analysis_rtg_vcfeval  };
    use MIP::Recipes::Analysis::Sacct qw{ analysis_sacct };
    use MIP::Recipes::Analysis::Sambamba_depth qw{ analysis_sambamba_depth };
    use MIP::Recipes::Analysis::Samtools_subsample_MT
      qw{ analysis_samtools_subsample_MT };
    use MIP::Recipes::Analysis::Sv_reformat qw{ analysis_sv_reformat };
    use MIP::Recipes::Analysis::Snpeff qw{ analysis_snpeff };
    use MIP::Recipes::Analysis::Sv_combinevariantcallsets
      qw{ analysis_sv_combinevariantcallsets };
    use MIP::Recipes::Analysis::Tiddit qw{ analysis_tiddit };
    use MIP::Recipes::Analysis::Variantannotationblock
      qw{ analysis_variantannotationblock };
    use MIP::Recipes::Analysis::Variant_integrity
      qw{ analysis_variant_integrity };
    use MIP::Recipes::Analysis::Vcf2cytosure qw{ analysis_vcf2cytosure };
    use MIP::Recipes::Analysis::Vep qw{ analysis_vep analysis_vep_sv };
    use MIP::Recipes::Analysis::Vt qw{ analysis_vt };
    use MIP::Recipes::Build::Rare_disease qw{build_rare_disease_meta_files};

    ### Build recipes
    $log->info(q{[Reference check - Reference prerequisites]});

    build_rare_disease_meta_files(
        {
            parameter_href          => $parameter_href,
            active_parameter_href   => $active_parameter_href,
            sample_info_href        => $sample_info_href,
            file_info_href          => $file_info_href,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            log                     => $log,
        }
    );

    ### Analysis recipes
    if ( $active_parameter_href->{pfastqc} ) {

        $log->info(q{[Fastqc]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $fastqc_program_name = q{fastqc};

            my $outsample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $fastqc_program_name );
            analysis_fastqc(
                {
                    parameter_href        => $parameter_href,
                    active_parameter_href => $active_parameter_href,
                    sample_info_href      => $sample_info_href,
                    infiles_ref           => \@{ $infile_href->{$sample_id} },

                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    insample_directory      => $indir_path_href->{$sample_id},
                    outsample_directory     => $outsample_directory,
                    sample_id               => $sample_id,
                    program_name            => $fastqc_program_name,
                }
            );
        }
    }

    if ( $active_parameter_href->{pbwa_mem} ) {

        $log->info(q{[BWA Mem]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $outsample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );

            analysis_bwa_mem(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infiles_ref             => \@{ $infile_href->{$sample_id} },
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    insample_directory      => $indir_path_href->{$sample_id},
                    outsample_directory     => $outsample_directory,
                    sample_id               => $sample_id,
                    program_name            => q{bwa_mem},
                }
            );
        }
    }
## Run consecutive models
    if ( $active_parameter_href->{reduce_io} ) {

        ## Enable as program
        $active_parameter_href->{pbamcalibrationblock} = 1;

        $log->info(q{[Bamcalibrationblock]});

        analysis_bamcalibrationblock(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                lane_href               => $lane_href,
                job_id_href             => $job_id_href,
                outaligner_dir => $active_parameter_href->{outaligner_dir},
                program_name   => q{bamcalibrationblock},
                log            => $log,
            }
        );

    }
    else {

        ## Always run merge even for single samples to rename them correctly for standardised downstream processing.
        if ( $active_parameter_href->{ppicardtools_mergesamfiles} ) {

            $log->info(q{[Picardtools mergesamfiles]});

          SAMPLE_ID:
            foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } )
            {

                my $insample_directory =
                  catdir( $active_parameter_href->{outdata_dir},
                    $sample_id, $active_parameter_href->{outaligner_dir} );
                my $outsample_directory =
                  catdir( $active_parameter_href->{outdata_dir},
                    $sample_id, $active_parameter_href->{outaligner_dir} );

                analysis_picardtools_mergesamfiles(
                    {
                        parameter_href          => $parameter_href,
                        active_parameter_href   => $active_parameter_href,
                        sample_info_href        => $sample_info_href,
                        file_info_href          => $file_info_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        lane_href               => $lane_href,
                        job_id_href             => $job_id_href,
                        insample_directory      => $insample_directory,
                        outsample_directory     => $outsample_directory,
                        sample_id               => $sample_id,
                        program_name            => q{picardtools_mergesamfiles},
                    }
                );
            }
        }

        if ( $active_parameter_href->{pmarkduplicates} ) {

            $log->info(q{[Markduplicates]});

          SAMPLE_ID:
            foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } )
            {

                ## Assign directories
                my $insample_directory =
                  catdir( $active_parameter_href->{outdata_dir},
                    $sample_id, $active_parameter_href->{outaligner_dir} );
                my $outsample_directory =
                  catdir( $active_parameter_href->{outdata_dir},
                    $sample_id, $active_parameter_href->{outaligner_dir} );

                analysis_markduplicates(
                    {
                        parameter_href          => $parameter_href,
                        active_parameter_href   => $active_parameter_href,
                        sample_info_href        => $sample_info_href,
                        file_info_href          => $file_info_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        job_id_href             => $job_id_href,
                        insample_directory      => $insample_directory,
                        outsample_directory     => $outsample_directory,
                        sample_id               => $sample_id,
                        program_name            => q{markduplicates},
                    }
                );
            }
        }

        if ( $active_parameter_href->{pgatk_realigner} ) {

            $log->info(q{[GATK realignertargetcreator/indelrealigner]});

          SAMPLE_ID:
            foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } )
            {

                ## Assign directories
                my $insample_directory =
                  catdir( $active_parameter_href->{outdata_dir},
                    $sample_id, $active_parameter_href->{outaligner_dir} );
                my $outsample_directory =
                  catdir( $active_parameter_href->{outdata_dir},
                    $sample_id, $active_parameter_href->{outaligner_dir} );

                analysis_gatk_realigner(
                    {
                        parameter_href          => $parameter_href,
                        active_parameter_href   => $active_parameter_href,
                        sample_info_href        => $sample_info_href,
                        file_info_href          => $file_info_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        job_id_href             => $job_id_href,
                        sample_id               => $sample_id,
                        insample_directory      => $insample_directory,
                        outsample_directory     => $outsample_directory,
                        program_name            => q{gatk_realigner},
                    }
                );
            }
        }

        if ( $active_parameter_href->{pgatk_baserecalibration} ) {

            $log->info(q{[GATK baserecalibrator/printreads]});

          SAMPLE_ID:
            foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } )
            {

                ## Assign directories
                my $insample_directory =
                  catdir( $active_parameter_href->{outdata_dir},
                    $sample_id, $active_parameter_href->{outaligner_dir} );
                my $outsample_directory =
                  catdir( $active_parameter_href->{outdata_dir},
                    $sample_id, $active_parameter_href->{outaligner_dir} );

                analysis_gatk_baserecalibration(
                    {
                        parameter_href          => $parameter_href,
                        active_parameter_href   => $active_parameter_href,
                        sample_info_href        => $sample_info_href,
                        file_info_href          => $file_info_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        job_id_href             => $job_id_href,
                        sample_id               => $sample_id,
                        insample_directory      => $insample_directory,
                        outsample_directory     => $outsample_directory,
                        program_name            => q{gatk_baserecalibration},
                    }
                );
            }
        }
    }
    if ( $active_parameter_href->{pchanjo_sexcheck} ) {

        $log->info(q{[Chanjo sexcheck]});

      SAMPLE_IDS:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );
            my $outsample_directory = catdir(
                $active_parameter_href->{outdata_dir},    $sample_id,
                $active_parameter_href->{outaligner_dir}, q{coveragereport}
            );

            analysis_chanjo_sex_check(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    sample_id               => $sample_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    program_name            => q{chanjo_sexcheck},
                }
            );
        }
    }
    if ( $active_parameter_href->{psamtools_subsample_mt} ) {

        $log->info(q{[Subsample MT]});

      SAMPLE_IDS:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );
            my $outsample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );

            analysis_samtools_subsample_MT(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    sample_id               => $sample_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    program_name            => q{samtools_subsample_mt},
                }
            );
        }
    }
    if ( $active_parameter_href->{psambamba_depth} ) {

        $log->info(q{[Sambamba depth]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );
            my $outsample_directory = catdir(
                $active_parameter_href->{outdata_dir},    $sample_id,
                $active_parameter_href->{outaligner_dir}, q{coveragereport}
            );

            analysis_sambamba_depth(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    sample_id               => $sample_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    program_name            => q{sambamba_depth},
                }
            );
        }
    }

    if ( $active_parameter_href->{pbedtools_genomecov} ) {

        $log->info(q{[Bedtools genomecov]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Assign directories
            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );
            my $outsample_directory = catdir(
                $active_parameter_href->{outdata_dir},    $sample_id,
                $active_parameter_href->{outaligner_dir}, q{coveragereport}
            );

            analysis_bedtools_genomecov(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    sample_id               => $sample_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    program_name            => q{bedtools_genomecov},
                }
            );
        }
    }
    if ( $active_parameter_href->{ppicardtools_collectmultiplemetrics} ) {

        $log->info(q{[Picardtools collectmultiplemetrics]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Assign directories
            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );
            my $outsample_directory = catdir(
                $active_parameter_href->{outdata_dir},    $sample_id,
                $active_parameter_href->{outaligner_dir}, q{coveragereport}
            );

            analysis_picardtools_collectmultiplemetrics(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    sample_id               => $sample_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    program_name => q{picardtools_collectmultiplemetrics},
                }
            );
        }
    }
    if ( $active_parameter_href->{ppicardtools_collecthsmetrics} ) {

        $log->info(q{[Picardtools collecthsmetrics]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Assign directories
            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );
            my $outsample_directory = catdir(
                $active_parameter_href->{outdata_dir},    $sample_id,
                $active_parameter_href->{outaligner_dir}, q{coveragereport}
            );

            analysis_picardtools_collecthsmetrics(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    sample_id               => $sample_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    program_name            => q{picardtools_collecthsmetrics},
                }
            );
        }
    }
## Run Rcovplot scripts
    if ( $active_parameter_href->{prcovplots} ) {

        if ( $active_parameter_href->{pbedtools_genomecov} > 0 ) {

            $log->info(q{[Rcovplots]});

            my $program_name = lc q{rcovplots};

          SAMPLE_ID:
            foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } )
            {

                ## Assign directories
                my $insample_directory =
                  catdir( $active_parameter_href->{outdata_dir},
                    $sample_id, $active_parameter_href->{outaligner_dir} );
                my $outsample_directory = catdir(
                    $active_parameter_href->{outdata_dir},
                    $sample_id, $active_parameter_href->{outaligner_dir},
                    q{coveragereport}
                );
                analysis_rcoverageplots(
                    {
                        parameter_href          => $parameter_href,
                        active_parameter_href   => $active_parameter_href,
                        sample_info_href        => $sample_info_href,
                        file_info_href          => $file_info_href,
                        lane_href               => $lane_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        job_id_href             => $job_id_href,
                        sample_id               => $sample_id,
                        insample_directory      => $insample_directory,
                        outsample_directory     => $outsample_directory,
                        program_name            => q{rcovplots},
                    }
                );
            }
        }
    }
    if ( $active_parameter_href->{pcnvnator} ) {

        $log->info(q{[CNVnator]});

        my $cnvnator_program_name = q{cnvnator};

        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Assign directories
            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );

            my $outsample_directory = catdir(
                $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir},
                $cnvnator_program_name
            );

            analysis_cnvnator(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    sample_id               => $sample_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    program_name            => $cnvnator_program_name,
                }
            );
        }
    }
    if ( $active_parameter_href->{pdelly_call} ) {

        $log->info(q{[Delly_call]});

        my $delly_program_name = q{delly_call};

        my $program_outdirectory_name =
          $parameter_href->{ q{p} . $delly_program_name }{outdir_name};

        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );
            my $outsample_directory = catdir(
                $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir},
                $program_outdirectory_name
            );

            analysis_delly_call(
                {
                    active_parameter_href   => $active_parameter_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    insample_directory      => $insample_directory,
                    job_id_href             => $job_id_href,
                    outsample_directory     => $outsample_directory,
                    parameter_href          => $parameter_href,
                    program_name            => $delly_program_name,
                    sample_id               => $sample_id,
                    sample_info_href        => $sample_info_href,
                }
            );
        }
    }
    if ( $active_parameter_href->{pdelly_reformat} ) {

        $log->info(q{[Delly_reformat]});

        my $delly_refrm_program_name = q{delly_reformat};

        my $program_outdirectory_name =
          $parameter_href->{ q{p} . $delly_refrm_program_name }{outdir_name};

        my $outfamily_directory = catfile(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
            $active_parameter_href->{outaligner_dir},
            $program_outdirectory_name
        );

        analysis_delly_reformat(
            {
                active_parameter_href   => $active_parameter_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                outfamily_directory     => $outfamily_directory,
                parameter_href          => $parameter_href,
                program_name            => $delly_refrm_program_name,
                sample_info_href        => $sample_info_href,
            }
        );
    }
    if ( $active_parameter_href->{pmanta} ) {

        $log->info(q{[Manta]});

        my $manta_program_name = q{manta};

        my $outfamily_directory = catfile(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
            $active_parameter_href->{outaligner_dir},
            $manta_program_name,
        );

        analysis_manta(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                program_name            => $manta_program_name,
                outfamily_directory     => $outfamily_directory,
            }
        );
    }
    if ( $active_parameter_href->{ptiddit} ) {

        $log->info(q{[Tiddit]});

        my $tiddit_program_name = q{tiddit};

        my $outfamily_directory = catfile(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
            $active_parameter_href->{outaligner_dir},
            $tiddit_program_name,
        );

        analysis_tiddit(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                program_name            => $tiddit_program_name,
                outfamily_directory     => $outfamily_directory,
            }
        );
    }
    if ( $active_parameter_href->{pexpansionhunter} ) {

        $log->info(q{[ExpansionHunter]});

        my $expansionhunter_program_name = q{expansionhunter};

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Assign directories
            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );

            my $outsample_directory = catdir(
                $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir},
                $expansionhunter_program_name
            );

            analysis_expansionhunter(
                {
                    active_parameter_href   => $active_parameter_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    insample_directory      => $insample_directory,
                    job_id_href             => $job_id_href,
                    outsample_directory     => $outsample_directory,
                    parameter_href          => $parameter_href,
                    program_name            => $cnvnator_program_name,
                    sample_id               => $sample_id,
                    sample_info_href        => $sample_info_href,
                }
            );
        }
    }
    if ( $active_parameter_href->{psv_combinevariantcallsets} ) {

        $log->info(q{[SV combinevariantcallsets]});

        analysis_sv_combinevariantcallsets(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                program_name            => q{sv_combinevariantcallsets},
            }
        );
    }
    if ( $active_parameter_href->{psv_varianteffectpredictor} ) {

        $log->info(q{[SV varianteffectpredictor]});

        analysis_vep_sv(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                contigs_ref             => \@{ $file_info_href->{contigs} },
                program_name            => q{sv_varianteffectpredictor},
            }
        );
    }
    if ( $active_parameter_href->{psv_vcfparser} ) {

        $log->info(q{[SV vcfparser]});

        analysis_sv_vcfparser(
            {
                active_parameter_href   => $active_parameter_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                parameter_href          => $parameter_href,
                program_name            => q{sv_vcfparser},
                sample_info_href        => $sample_info_href,
            }
        );
    }
    if ( $active_parameter_href->{psv_rankvariant} ) {

        $log->info(q{[SV rankvariant]});

        if ( defined $parameter_href->{dynamic_parameter}{unaffected}
            && @{ $parameter_href->{dynamic_parameter}{unaffected} } eq
            @{ $active_parameter_href->{sample_ids} } )
        {

            $log->warn(
q{Only unaffected sample(s) in pedigree - skipping genmod 'models', 'score' and 'compound'}
            );

            analysis_sv_rankvariant_unaffected(
                {
                    active_parameter_href   => $active_parameter_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    parameter_href          => $parameter_href,
                    program_name            => q{sv_rankvariant},
                    sample_info_href        => $sample_info_href,
                }
            );
        }
        else {
            analysis_sv_rankvariant(
                {
                    active_parameter_href   => $active_parameter_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    parameter_href          => $parameter_href,
                    program_name            => q{sv_rankvariant},
                    sample_info_href        => $sample_info_href,
                }
            );
        }
    }
    if ( $active_parameter_href->{psv_reformat} ) {

        $log->info(q{[SV reformat]});

        analysis_sv_reformat(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                program_name            => q{sv_reformat},
            }
        );
    }
    if ( $active_parameter_href->{pvcf2cytosure} ) {

        $log->info(q{[Vcf2cytosure]});

        my $v2cs_program_name = q{vcf2cytosure};

        my $outfamily_directory = catfile(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
            $active_parameter_href->{outaligner_dir},
            $v2cs_program_name,
        );

        analysis_vcf2cytosure(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                outfamily_directory     => $outfamily_directory,
                program_name            => $v2cs_program_name,
            }
        );
    }
    if ( $active_parameter_href->{pbcftools_mpileup} ) {

        $log->info(q{[Bcftools mpileup]});

        my $program_outdirectory_name =
          $parameter_href->{pbcftools_mpileup}{outdir_name};

        my $outfamily_directory = catfile(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
            $active_parameter_href->{outaligner_dir},
            $program_outdirectory_name
        );

        analysis_bcftools_mpileup(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                program_name            => q{bcftools_mpileup},
                outfamily_directory     => $outfamily_directory,
            }
        );
    }
    if ( $active_parameter_href->{pfreebayes} ) {

        $log->info(q{[Freebayes]});

        my $program_outdirectory_name =
          $parameter_href->{pfreebayes}{outdir_name};

        my $outfamily_directory = catfile(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
            $active_parameter_href->{outaligner_dir},
            $program_outdirectory_name
        );

        analysis_freebayes_calling(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                program_name            => q{freebayes},
                outfamily_directory     => $outfamily_directory,
            }
        );
    }
    if ( $active_parameter_href->{pgatk_haplotypecaller} ) {

        $log->info(q{[GATK haplotypecaller]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );
            my $outsample_directory = catdir(
                $active_parameter_href->{outdata_dir},    $sample_id,
                $active_parameter_href->{outaligner_dir}, q{gatk}
            );

            analysis_gatk_haplotypecaller(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    sample_id               => $sample_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    program_name            => q{gatk_haplotypecaller},
                }
            );
        }
    }
    if ( $active_parameter_href->{pgatk_genotypegvcfs} ) {

        $log->info(q{[GATK genotypegvcfs]});

        my $family_analysis_directory = catfile(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
            $active_parameter_href->{outaligner_dir}, q{gatk},
        );

        my $outfamily_file_directory = catdir(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
        );

        analysis_gatk_genotypegvcfs(
            {
                parameter_href           => $parameter_href,
                active_parameter_href    => $active_parameter_href,
                sample_info_href         => $sample_info_href,
                file_info_href           => $file_info_href,
                infile_lane_prefix_href  => $infile_lane_prefix_href,
                job_id_href              => $job_id_href,
                program_name             => q{gatk_genotypegvcfs},
                family_id                => $active_parameter_href->{family_id},
                outfamily_directory      => $family_analysis_directory,
                outfamily_file_directory => $outfamily_file_directory,
            }
        );

        $log->info(q{[GATK concatenate genotypegvcfs files]});

        analysis_gatk_concatenate_genotypegvcfs(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                infamily_directory      => $family_analysis_directory,
                outfamily_directory     => $family_analysis_directory,
                program_name            => q{gatk_genotypegvcfs},
            }
        );
    }
    if ( $active_parameter_href->{pgatk_variantrecalibration} ) {

        $log->info(q{[GATK variantrecalibrator/applyrecalibration]});

        my $program_outdirectory_name =
          $parameter_href->{pgatk_variantrecalibration}{outdir_name};

        my $infamily_directory = catfile(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
            $active_parameter_href->{outaligner_dir},
            $program_outdirectory_name
        );
        my $outfamily_directory = $infamily_directory;
        my $consensus_analysis_type =
          $parameter_href->{dynamic_parameter}{consensus_analysis_type};

        if ( $consensus_analysis_type eq q{wes} ) {

            analysis_gatk_variantrecalibration_wes(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    infamily_directory      => $infamily_directory,
                    outfamily_directory     => $outfamily_directory,
                    program_name            => q{gatk_variantrecalibration},
                }
            );
        }
        else {

            ## WGS and WES/WGS
            analysis_gatk_variantrecalibration_wgs(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    infamily_directory      => $infamily_directory,
                    outfamily_directory     => $outfamily_directory,
                    program_name            => q{gatk_variantrecalibration},
                }
            );
        }
    }
    if ( $active_parameter_href->{pgatk_combinevariantcallsets} ) {

        $log->info(q{[GATK combinevariantcallsets]});

        my $outfamily_directory = catfile(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
            $active_parameter_href->{outaligner_dir}
        );
        analysis_gatk_combinevariantcallsets(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                outfamily_directory     => $outfamily_directory,
                program_name            => q{gatk_combinevariantcallsets},
            }
        );
    }
    if ( $active_parameter_href->{ppeddy} ) {

        $log->info(q{[Peddy]});

        my $peddy_program_name = q{peddy};

        my $infamily_directory = catdir(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
            $active_parameter_href->{outaligner_dir}
        );

        my $outfamily_directory = catfile(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
            $active_parameter_href->{outaligner_dir},
            q{casecheck},
            $peddy_program_name
        );

        analysis_peddy(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                program_name            => $peddy_program_name,
                infamily_directory      => $infamily_directory,
                outfamily_directory     => $outfamily_directory,
            }
        );
    }
    if ( $active_parameter_href->{pplink} ) {

        $log->info(q{[Plink]});

        my $plink_program_name = q{plink};

        my $infamily_directory = catdir(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
            $active_parameter_href->{outaligner_dir}
        );

        my $outfamily_directory = catfile(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
            $active_parameter_href->{outaligner_dir},
            q{casecheck},
            $plink_program_name
        );

        analysis_plink(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                program_name            => $plink_program_name,
                infamily_directory      => $infamily_directory,
                outfamily_directory     => $outfamily_directory,
            }
        );
    }
    if ( $active_parameter_href->{pvariant_integrity} ) {

        $log->info(q{[Variant_integrity]});

        my $variant_inty_program_name = q{variant_integrity};

        my $infamily_directory = catdir(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
            $active_parameter_href->{outaligner_dir}
        );
        my $outfamily_directory = catfile(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
            $active_parameter_href->{outaligner_dir},
            q{casecheck},
            $variant_inty_program_name
        );

        analysis_variant_integrity(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                program_name            => $variant_inty_program_name,
                infamily_directory      => $infamily_directory,
                outfamily_directory     => $outfamily_directory,
            }
        );
    }
    if ( $active_parameter_href->{prtg_vcfeval} ) {

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            if ( $sample_id =~ /$active_parameter_href->{nist_id}/sxm ) {

                $log->info(q{[Rtg evaluation]});

                my $rtg_program_name = q{rtg_vcfeval};

                ## Assign directories
                my $infamily_directory = catdir(
                    $active_parameter_href->{outdata_dir},
                    $active_parameter_href->{family_id},
                    $active_parameter_href->{outaligner_dir}
                );
                my $outfamily_directory = catfile(
                    $active_parameter_href->{outdata_dir},
                    $active_parameter_href->{family_id},
                    $active_parameter_href->{outaligner_dir},
                    $rtg_program_name
                );
                analysis_rtg_vcfeval(
                    {
                        parameter_href          => $parameter_href,
                        active_parameter_href   => $active_parameter_href,
                        sample_info_href        => $sample_info_href,
                        file_info_href          => $file_info_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        job_id_href             => $job_id_href,
                        sample_id               => $sample_id,
                        call_type               => q{BOTH},
                        infamily_directory      => $infamily_directory,
                        outfamily_directory     => $outfamily_directory,
                        program_name            => $rtg_program_name,
                    }
                );
            }
        }
    }
    if ( $active_parameter_href->{pevaluation} ) {

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            if ( $sample_id =~ /$active_parameter_href->{nist_id}/sxm ) {

                $log->info(q{[Evaluation]});

                my $evolution_program_name = q{evaluation};

                ## Assign directories
                my $infamily_directory = catdir(
                    $active_parameter_href->{outdata_dir},
                    $active_parameter_href->{family_id},
                    $active_parameter_href->{outaligner_dir}
                );
                my $outfamily_directory = catfile(
                    $active_parameter_href->{outdata_dir},
                    $active_parameter_href->{family_id},
                    $active_parameter_href->{outaligner_dir},
                    $evolution_program_name
                );
                analysis_picardtools_genotypeconcordance(
                    {
                        parameter_href          => $parameter_href,
                        active_parameter_href   => $active_parameter_href,
                        sample_info_href        => $sample_info_href,
                        file_info_href          => $file_info_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        job_id_href             => $job_id_href,
                        sample_id               => $sample_id,
                        call_type               => q{BOTH},
                        infamily_directory      => $infamily_directory,
                        outfamily_directory     => $outfamily_directory,
                        program_name            => $evolution_program_name,
                    }
                );
            }
        }
    }
    if ( $active_parameter_href->{pgatk_variantevalall} ) {

        $log->info(q{[GATK variantevalall]});

        my $variant_eval_all_program_name = q{gatk_variantevalall};
      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Assign directories
            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );

            my $outsample_directory = catdir(
                $active_parameter_href->{outdata_dir},
                $sample_id,
                $active_parameter_href->{outaligner_dir},
                $variant_eval_all_program_name
            );

            analysis_gatk_variantevalall(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    sample_id               => $sample_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    program_name            => $variant_eval_all_program_name,
                }
            );
        }
    }

### If no males or other remove contig Y from all downstream analysis
    my @file_info_contig_keys = (qw{ contigs_size_ordered contigs });

  KEY:
    foreach my $key (@file_info_contig_keys) {

        ## Removes contig_names from contigs array if no male or 'other' found
        @{ $file_info_href->{$key} } = delete_male_contig(
            {
                contigs_ref => \@{ $file_info_href->{$key} },
                found_male  => $active_parameter_href->{found_male},
            }
        );
    }

    if ( $active_parameter_href->{reduce_io} ) {

        $active_parameter_href->{pvariantannotationblock} =
          1;    # Enable as program

        $log->info(q{[Variantannotationblock]});

        analysis_variantannotationblock(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                outaligner_dir => $active_parameter_href->{outaligner_dir},
                call_type      => q{BOTH},
                program_name   => q{variantannotationblock},
            }
        );
    }
    else {

        if ( $active_parameter_href->{pprepareforvariantannotationblock} ) {

            $log->info(q{[Prepareforvariantannotationblock]});

            analysis_prepareforvariantannotationblock(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    call_type               => q{BOTH},
                    program_name => q{prepareforvariantannotationblock},
                }
            );
        }

        if ( $active_parameter_href->{prhocall} ) {

            $log->info(q{[Rhocall]});

            my $infamily_directory = catdir(
                $active_parameter_href->{outdata_dir},
                $active_parameter_href->{family_id},
                $active_parameter_href->{outaligner_dir}
            );
            my $outfamily_directory = $infamily_directory;

            analysis_rhocall_annotate(
                {
                    active_parameter_href   => $active_parameter_href,
                    call_type               => q{BOTH},
                    file_info_href          => $file_info_href,
                    infamily_directory      => $infamily_directory,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    outfamily_directory     => $outfamily_directory,
                    parameter_href          => $parameter_href,
                    program_name            => q{rhocall},
                    sample_info_href        => $sample_info_href,
                }
            );
        }

        if ( $active_parameter_href->{pvt} ) {

            $log->info(q{[Vt]});

            my $infamily_directory = catdir(
                $active_parameter_href->{outdata_dir},
                $active_parameter_href->{family_id},
                $active_parameter_href->{outaligner_dir}
            );
            my $outfamily_directory = $infamily_directory;

            analysis_vt(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    infamily_directory      => $infamily_directory,
                    job_id_href             => $job_id_href,
                    call_type               => q{BOTH},
                    outfamily_directory     => $outfamily_directory,
                    program_name            => q{vt},
                }
            );
        }
        if ( $active_parameter_href->{pfrequency_filter} ) {

            $log->info(q{[Frequency_filter]});

            my $infamily_directory = catdir(
                $active_parameter_href->{outdata_dir},
                $active_parameter_href->{family_id},
                $active_parameter_href->{outaligner_dir}
            );
            my $outfamily_directory = $infamily_directory;

            analysis_frequency_filter(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    infamily_directory      => $infamily_directory,
                    job_id_href             => $job_id_href,
                    call_type               => q{BOTH},
                    outfamily_directory     => $outfamily_directory,
                    program_name            => q{frequency_filter},
                }
            );
        }

        if ( $active_parameter_href->{pvarianteffectpredictor} ) {

            $log->info(q{[Varianteffectpredictor]});

            analysis_vep(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    call_type               => q{BOTH},
                    program_name            => q{varianteffectpredictor},
                }
            );
        }
        if ( $active_parameter_href->{pvcfparser} ) {

            $log->info(q{[Vcfparser]});

            my $infamily_directory = catdir(
                $active_parameter_href->{outdata_dir},
                $active_parameter_href->{family_id},
                $active_parameter_href->{outaligner_dir}
            );
            my $outfamily_directory = $infamily_directory;

            analysis_mip_vcfparser(
                {
                    active_parameter_href   => $active_parameter_href,
                    call_type               => q{BOTH},
                    file_info_href          => $file_info_href,
                    infamily_directory      => $infamily_directory,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    outfamily_directory     => $outfamily_directory,
                    parameter_href          => $parameter_href,
                    program_name            => q{vcfparser},
                    sample_info_href        => $sample_info_href,
                }
            );
        }

        if ( $active_parameter_href->{psnpeff} ) {

            $log->info(q{[Snpeff]});

            my $infamily_directory = catdir(
                $active_parameter_href->{outdata_dir},
                $active_parameter_href->{family_id},
                $active_parameter_href->{outaligner_dir}
            );
            my $outfamily_directory = $infamily_directory;

            analysis_snpeff(
                {
                    active_parameter_href   => $active_parameter_href,
                    call_type               => q{BOTH},
                    file_info_href          => $file_info_href,
                    infamily_directory      => $infamily_directory,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    outfamily_directory     => $outfamily_directory,
                    parameter_href          => $parameter_href,
                    program_name            => q{snpeff},
                    sample_info_href        => $sample_info_href,
                }
            );
        }
        if ( $active_parameter_href->{prankvariant} ) {

            $log->info(q{[Rankvariant]});

            if ( defined $parameter_href->{dynamic_parameter}{unaffected}
                && @{ $parameter_href->{dynamic_parameter}{unaffected} } eq
                @{ $active_parameter_href->{sample_ids} } )
            {

                $log->warn(
q{Only unaffected sample in pedigree - skipping genmod 'models', 'score' and 'compound'}
                );

                analysis_rankvariant_unaffected(
                    {
                        parameter_href          => $parameter_href,
                        active_parameter_href   => $active_parameter_href,
                        sample_info_href        => $sample_info_href,
                        file_info_href          => $file_info_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        job_id_href             => $job_id_href,
                        call_type               => q{BOTH},
                        program_name            => q{rankvariant},
                    }
                );
            }
            else {

                analysis_rankvariant(
                    {
                        parameter_href          => $parameter_href,
                        active_parameter_href   => $active_parameter_href,
                        sample_info_href        => $sample_info_href,
                        file_info_href          => $file_info_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        job_id_href             => $job_id_href,
                        call_type               => q{BOTH},
                        program_name            => q{rankvariant},
                    }
                );
            }
        }
        if ( $active_parameter_href->{pendvariantannotationblock} ) {

            $log->info(q{[Endvariantannotationblock]});

            analysis_endvariantannotationblock(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    call_type               => q{BOTH},
                    program_name            => q{endvariantannotationblock},
                }
            );
        }
    }

    if ( $active_parameter_href->{pgatk_variantevalexome} ) {

        $log->info(q{[GATK variantevalexome]});

        my $variant_eval_exome_program_name = q{gatk_variantevalexome};

        ## Assign directories
        my $infamily_directory = catdir(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id},
            $active_parameter_href->{outaligner_dir}
        );

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Assign directories
            my $outsample_directory = catdir(
                $active_parameter_href->{outdata_dir},
                $sample_id,
                $active_parameter_href->{outaligner_dir},
                $variant_eval_exome_program_name
            );

            analysis_gatk_variantevalexome(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    sample_id               => $sample_id,
                    infamily_directory      => $infamily_directory,
                    outsample_directory     => $outsample_directory,
                    program_name            => $variant_eval_exome_program_name,
                }
            );
        }
    }
    if ( $active_parameter_href->{pqccollect} ) {

        $log->info(q{[Qccollect]});

        my $outfamily_directory = catdir(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id}
        );

        analysis_qccollect(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                program_name            => q{qccollect},
                outfamily_directory     => $outfamily_directory,
            }
        );
    }

    if ( $active_parameter_href->{pmultiqc} ) {

        $log->info(q{[Multiqc]});

        analysis_multiqc(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                program_name            => q{multiqc},
            }
        );
    }

    if ( $active_parameter_href->{panalysisrunstatus} ) {

        $log->info(q{[Analysis run status]});

        analysis_analysisrunstatus(
            {
                active_parameter_href   => $active_parameter_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                parameter_href          => $parameter_href,
                program_name            => q{analysisrunstatus},
                sample_info_href        => $sample_info_href,
            }
        );
    }
    if ( $active_parameter_href->{psacct} ) {

        $log->info(q{[Sacct]});

        analysis_sacct(
            {
                active_parameter_href   => $active_parameter_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                parameter_href          => $parameter_href,
                program_name            => q{sacct},
                sample_info_href        => $sample_info_href,
            }
        );
    }
    return;
}

1;
