package MIP::Recipes::Analysis::Gatk_variantrecalibration;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw { uniq };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $ASTERISK $AMPERSAND $COLON $DASH $DOT $LOG_NAME $NEWLINE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ analysis_gatk_variantrecalibration_wes analysis_gatk_variantrecalibration_wgs };

}

sub analysis_gatk_variantrecalibration_wes {

## Function : GATK VariantRecalibrator/ApplyRecalibration analysis recipe for wes and panel data
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        profile_base_command => {
            default     => q{sbatch},
            store       => \$profile_base_command,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Contigs qw{ delete_contig_elements };
    use MIP::File_info qw{ get_io_files };
    use MIP::Pedigree qw{ create_fam_file };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bcftools qw{ bcftools_norm };
    use MIP::Program::Gatk
      qw{ gatk_applyvqsr gatk_calculategenotypeposteriors gatk_selectvariants gatk_variantrecalibrator };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constants
    Readonly my $JAVA_MEMORY_ALLOCATION => 10;
    Readonly my $MAX_GAUSSIAN_LEVEL     => 4;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my %io = get_io_files(
        {
            id             => $case_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
        }
    );

    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path        = $io{in}{file_path};

    my $enable_indel_max_gaussians_filter =
      $active_parameter_href->{gatk_variantrecalibration_indel_max_gaussians};
    my $enable_snv_max_gaussians_filter =
      $active_parameter_href->{gatk_variantrecalibration_snv_max_gaussians};
    my $referencefile_path  = $active_parameter_href->{human_genome_reference};
    my $resource_indel_href = $active_parameter_href->{gatk_variantrecalibration_resource_indel};
    my $resource_snv_href   = $active_parameter_href->{gatk_variantrecalibration_resource_snv};
    my %recipe              = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $recipe{job_id_chain},
                id                     => $case_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => [$infile_name_prefix],
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );
    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my $outfile_path        = $io{out}{file_path};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe{core_number},
            directory_id          => $case_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe{memory},
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            temp_directory        => $temp_directory,
        }
    );

    ### Split to enable submission to &sample_info_qc later
    #my ( $volume, $directory, $stderr_file ) =
    #  splitpath( $recipe_info_path . $DOT . q{stderr.txt} );

    ## Create .fam file to be used in variant calling analyses
    my $fam_file_path = catfile( $outdir_path_prefix, $case_id . $DOT . q{fam} );
    create_fam_file(
        {
            case_id          => $case_id,
            execution_mode   => q{system},
            fam_file_path    => $fam_file_path,
            filehandle       => $filehandle,
            sample_ids_ref   => $active_parameter_href->{sample_ids},
            sample_info_href => $sample_info_href,
        }
    );

    ## Check if "--pedigree" and "--pedigreeValidationType" should be included in analysis
    my %commands = gatk_pedigree_flag(
        {
            fam_file_path => $fam_file_path,
        }
    );

    ### GATK VariantRecalibrator
    ## Set mode to be used in variant recalibration
    # Exome will be processed using mode BOTH since there are to few INDELS
    # to use in the recalibration model even though using 30 exome BAMS in
    # Haplotypecaller step.
    my $mode = q{BOTH};

    my $select_infile_path;
    my $norm_infile_path;

    say {$filehandle} q{## GATK VariantRecalibrator};

    ## Get parameters
    my $max_gaussian_level;
    my @ts_tranches = @{ $active_parameter_href->{gatk_variantrecalibration_ts_tranches} };
    my @annotations = @{ $active_parameter_href->{gatk_variantrecalibration_annotations} };

    ### Special case: Not to be used with hybrid capture
    ## Removes an element from array and return new array while leaving orginal contigs_ref untouched
    @annotations = delete_contig_elements(
        {
            contigs_ref        => \@annotations,
            remove_contigs_ref => [qw{ DP }],
        }
    );
    my @snv_resources = _build_gatk_resource_command( { resources_href => $resource_snv_href, } );
    my @indel_resources =
      _build_gatk_resource_command( { resources_href => $resource_indel_href, } );

    # Create distinct set i.e. no duplicates.
    my @resources = uniq( @indel_resources, @snv_resources );

    ## Use hard filtering
    if (   $enable_snv_max_gaussians_filter
        || $enable_indel_max_gaussians_filter )
    {

        $max_gaussian_level = $MAX_GAUSSIAN_LEVEL;
    }

    my $recal_file_path = $outfile_path_prefix . $DOT . q{intervals};
    gatk_variantrecalibrator(
        {
            annotations_ref      => \@annotations,
            filehandle           => $filehandle,
            infile_path          => $infile_path,
            java_use_large_pages => $active_parameter_href->{java_use_large_pages},
            max_gaussian_level   => $max_gaussian_level,
            memory_allocation    => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
            mode                 => $mode,
            outfile_path         => $recal_file_path,
            referencefile_path   => $referencefile_path,
            resources_ref        => \@resources,
            rscript_file_path    => $recal_file_path . $DOT . q{plots.R},
            temp_directory       => $temp_directory,
            tranches_file_path   => $recal_file_path . $DOT . q{tranches},
            ts_tranches_ref      => \@ts_tranches,
            verbosity            => $active_parameter_href->{gatk_logging_level},
        }
    );
    say {$filehandle} $NEWLINE;

    ## GATK ApplyVQSR
    say {$filehandle} q{## GATK ApplyVQSR};

    ## Get parameters
    my $ts_filter_level;
    ## Exome and panel analysis use combined reference for more power

    ## Infile genotypegvcfs combined vcf which used reference gVCFs to create combined vcf file
    $ts_filter_level = $active_parameter_href->{gatk_variantrecalibration_snv_tsfilter_level};

    my $apply_vqsr_outfile_path = $outfile_path_prefix . $UNDERSCORE . q{apply} . $outfile_suffix;
    gatk_applyvqsr(
        {
            filehandle           => $filehandle,
            infile_path          => $infile_path,
            java_use_large_pages => $active_parameter_href->{java_use_large_pages},
            memory_allocation    => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
            mode                 => $mode,
            outfile_path         => $apply_vqsr_outfile_path,
            recal_file_path      => $recal_file_path,
            referencefile_path   => $referencefile_path,
            temp_directory       => $temp_directory,
            tranches_file_path   => $recal_file_path . $DOT . q{tranches},
            ts_filter_level      => $ts_filter_level,
            verbosity            => $active_parameter_href->{gatk_logging_level},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Set infiles for next step(s)
    $select_infile_path = $apply_vqsr_outfile_path;
    $norm_infile_path   = $apply_vqsr_outfile_path;

    if ( not $active_parameter_href->{gatk_variantrecalibration_keep_unnormalised} ) {

        ## Bcftools norm, left-align and normalize indels, split multiallelics
        my $norm_outfile_path =
          $outfile_path_prefix . $UNDERSCORE . q{normalized} . $outfile_suffix;
        bcftools_norm(
            {
                filehandle      => $filehandle,
                infile_path     => $norm_infile_path,
                multiallelic    => $DASH,
                outfile_path    => $norm_outfile_path,
                output_type     => q{v},
                reference_path  => $referencefile_path,
                stderrfile_path => $outfile_path_prefix . $UNDERSCORE . q{normalized.stderr},
            }
        );
        ## Set outfile path for next step
        $select_infile_path = $norm_outfile_path;

        say {$filehandle} $NEWLINE;
    }

    ### GATK SelectVariants

    ## Removes all genotype information for exome ref and recalulates meta-data info for remaining samples in new file.
    # Exome analysis

    say {$filehandle} q{## GATK SelectVariants};

    gatk_selectvariants(
        {
            filehandle           => $filehandle,
            exclude_non_variants => 1,
            infile_path          => $select_infile_path,
            java_use_large_pages => $active_parameter_href->{java_use_large_pages},
            memory_allocation    => q{Xmx2g},
            outfile_path         => $outfile_path,
            referencefile_path   => $referencefile_path,
            sample_names_ref     => \@{ $active_parameter_href->{sample_ids} },
            temp_directory       => $temp_directory,
            verbosity            => $active_parameter_href->{gatk_logging_level},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Genotype refinement
    if (    $parameter_href->{cache}{trio}
        and $active_parameter_href->{gatk_calculategenotypeposteriors} )
    {

        say {$filehandle} q{## GATK CalculateGenotypePosteriors};

        my $calculategt_outfile_path =
          $outfile_path_prefix . $UNDERSCORE . q{refined} . $outfile_suffix;
        gatk_calculategenotypeposteriors(
            {
                filehandle                 => $filehandle,
                infile_path                => $outfile_path,
                java_use_large_pages       => $active_parameter_href->{java_use_large_pages},
                memory_allocation          => q{Xmx6g},
                num_ref_samples_if_no_call =>
                  $active_parameter_href->{gatk_num_reference_samples_if_no_call},
                outfile_path                 => $calculategt_outfile_path,
                pedigree                     => $commands{pedigree},
                referencefile_path           => $referencefile_path,
                supporting_callset_file_path =>
                  $active_parameter_href->{gatk_calculate_genotype_call_set},
                temp_directory => $temp_directory,
                verbosity      => $active_parameter_href->{gatk_logging_level},
            }
        );
        say {$filehandle} $NEWLINE;

        ## Change name of file to accomodate downstream
        gnu_mv(
            {
                filehandle   => $filehandle,
                infile_path  => $calculategt_outfile_path,
                outfile_path => $outfile_path,
            }
        );
        say {$filehandle} $NEWLINE;
    }

    if ( not $active_parameter_href->{gatk_variantrecalibration_keep_unnormalised} ) {

        ## BcfTools norm, Left-align and normalize indels, split multiallelics
        my $selected_norm_outfile_path =
          $outfile_path_prefix . $UNDERSCORE . q{selected_normalized} . $outfile_suffix;
        bcftools_norm(
            {
                filehandle     => $filehandle,
                infile_path    => $outfile_path,
                multiallelic   => $DASH,
                outfile_path   => $selected_norm_outfile_path,
                output_type    => q{v},
                reference_path => $referencefile_path,
            }
        );
        say {$filehandle} $NEWLINE;

        ## Change name of file to accomodate downstream
        gnu_mv(
            {
                filehandle   => $filehandle,
                infile_path  => $selected_norm_outfile_path,
                outfile_path => $outfile_path,
            }
        );
        say {$filehandle} $NEWLINE;
    }
    close $filehandle;

    if ( $recipe{mode} == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_case},
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
                log                               => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_ids_ref     => \@{ $active_parameter_href->{sample_ids} },
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub analysis_gatk_variantrecalibration_wgs {

## Function : GATK VariantRecalibrator/ApplyRecalibration analysis recipe for wgs data
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $profile_base_command;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        profile_base_command => {
            default     => q{sbatch},
            store       => \$profile_base_command,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Contigs qw{ delete_contig_elements };
    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Pedigree qw{ create_fam_file gatk_pedigree_flag };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bcftools qw{ bcftools_norm };
    use MIP::Program::Gatk
      qw{ gatk_applyvqsr gatk_calculategenotypeposteriors gatk_selectvariants gatk_variantrecalibrator };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info
      qw{ set_recipe_outfile_in_sample_info set_processing_metafile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constants
    Readonly my $JAVA_MEMORY_ALLOCATION               => 24;
    Readonly my $MAX_GAUSSIAN_LEVEL_INDEL             => 4;
    Readonly my $MAX_GAUSSIAN_LEVEL_SNV               => 6;
    Readonly my $MAX_GAUSSIAN_LEVEL_SNV_SINGLE_SAMPLE => 4;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my %io = get_io_files(
        {
            id             => $case_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
        }
    );

    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path        = $io{in}{file_path};

    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my $enable_indel_max_gaussians_filter =
      $active_parameter_href->{gatk_variantrecalibration_indel_max_gaussians};
    my $enable_snv_max_gaussians_filter =
      $active_parameter_href->{gatk_variantrecalibration_snv_max_gaussians};
    my $referencefile_path  = $active_parameter_href->{human_genome_reference};
    my $resource_indel_href = $active_parameter_href->{gatk_variantrecalibration_resource_indel};
    my $resource_snv_href   = $active_parameter_href->{gatk_variantrecalibration_resource_snv};
    my %recipe              = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $recipe{job_id_chain},
                id                     => $case_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => [$infile_name_prefix],
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );
    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my $outfile_path        = $io{out}{file_path};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe{core_number},
            directory_id          => $case_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe{memory},
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            temp_directory        => $temp_directory,
        }
    );

    ## Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) =
      splitpath( $recipe_info_path . $DOT . q{stderr.txt} );

    ## Create .fam file to be used in variant calling analyses
    my $fam_file_path = catfile( $outdir_path_prefix, $case_id . $DOT . q{fam} );
    create_fam_file(
        {
            case_id          => $case_id,
            execution_mode   => q{system},
            fam_file_path    => $fam_file_path,
            filehandle       => $filehandle,
            sample_ids_ref   => $active_parameter_href->{sample_ids},
            sample_info_href => $sample_info_href,
        }
    );

    ## Check if "--pedigree" should be included in analysis
    my %commands = gatk_pedigree_flag(
        {
            fam_file_path => $fam_file_path,
        }
    );

    ### GATK VariantRecalibrator
    ## Set mode to be used in variant recalibration
    # SNP and INDEL will be recalibrated successively in the same file
    # because when you specify eg SNP mode, the indels are emitted
    # without modification, and vice-versa.
    my @modes = qw{ SNP INDEL };

  MODE:
    foreach my $mode (@modes) {

        say {$filehandle} q{## GATK VariantRecalibrator};

        ## Get parameters
        my $max_gaussian_level;
        my $varrecal_infile_path;
        my @ts_tranches = @{ $active_parameter_href->{gatk_variantrecalibration_ts_tranches} };

        if ( $mode eq q{SNP} ) {

            $varrecal_infile_path = $infile_path;

            ## Use hard filtering
            if ($enable_snv_max_gaussians_filter) {

                ## Use fewer Gaussians for single sample cases
                if ( scalar @{ $active_parameter_href->{sample_ids} } == 1 ) {

                    $max_gaussian_level = $MAX_GAUSSIAN_LEVEL_SNV_SINGLE_SAMPLE;
                }
                else {

                    $max_gaussian_level = $MAX_GAUSSIAN_LEVEL_SNV;
                }
            }
        }

        ## Use created recalibrated snp vcf as input
        if ( $mode eq q{INDEL} ) {

            $varrecal_infile_path = $outfile_path_prefix . $DOT . q{SNV} . $outfile_suffix;

            ## Use hard filtering
            if ($enable_indel_max_gaussians_filter) {

                $max_gaussian_level = $MAX_GAUSSIAN_LEVEL_INDEL;
            }
        }

        my @annotations = @{ $active_parameter_href->{gatk_variantrecalibration_annotations} };

        ## Special case: Not to be used with hybrid capture. NOTE: Disabled when analysing wes + wgs in the same run
        if ( $consensus_analysis_type ne q{wgs} ) {

            ## Removes an element from array and return new array while leaving orginal contigs_ref untouched
            @annotations = delete_contig_elements(
                {
                    contigs_ref        => \@annotations,
                    remove_contigs_ref => [qw{ DP }],
                }
            );
        }

        my @resources;
        if ( $mode eq q{SNP} ) {

            @resources = _build_gatk_resource_command( { resources_href => $resource_snv_href, } );
        }
        if ( $mode eq q{INDEL} ) {

            @resources =
              _build_gatk_resource_command( { resources_href => $resource_indel_href, } );

            ## MQ annotation not part of GATK BP for INDEL
            @annotations = delete_contig_elements(
                {
                    contigs_ref        => \@annotations,
                    remove_contigs_ref => [qw{ MQ }],
                }
            );
        }

        my $recal_file_path = $outfile_path_prefix . $DOT . q{intervals};
        gatk_variantrecalibrator(
            {
                annotations_ref       => \@annotations,
                filehandle            => $filehandle,
                infile_path           => $varrecal_infile_path,
                java_use_large_pages  => $active_parameter_href->{java_use_large_pages},
                max_gaussian_level    => $max_gaussian_level,
                memory_allocation     => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
                mode                  => $mode,
                outfile_path          => $recal_file_path,
                referencefile_path    => $referencefile_path,
                resources_ref         => \@resources,
                rscript_file_path     => $recal_file_path . $DOT . q{plots.R},
                temp_directory        => $temp_directory,
                tranches_file_path    => $recal_file_path . $DOT . q{tranches},
                ts_tranches_ref       => \@ts_tranches,
                trust_all_polymorphic =>
                  $active_parameter_href->{gatk_variantrecalibration_trust_all_polymorphic},
                verbosity => $active_parameter_href->{gatk_logging_level},
            }
        );
        say {$filehandle} $NEWLINE;

        ## GATK ApplyVQSR
        say {$filehandle} q{## GATK ApplyVQSR};

        ## Get parameters
        my $applyvqsr_infile_path;
        my $applyvqsr_outfile_path;
        my $ts_filter_level;

        if ( $mode eq q{SNP} ) {

            $applyvqsr_infile_path  = $varrecal_infile_path;
            $applyvqsr_outfile_path = $outfile_path_prefix . $DOT . q{SNV} . $outfile_suffix;
            $ts_filter_level =
              $active_parameter_href->{gatk_variantrecalibration_snv_tsfilter_level};
        }

        ## Use created recalibrated snp vcf as input
        if ( $mode eq q{INDEL} ) {

            $applyvqsr_infile_path  = $outfile_path_prefix . $DOT . q{SNV} . $outfile_suffix;
            $applyvqsr_outfile_path = $outfile_path;
            $ts_filter_level =
              $active_parameter_href->{gatk_variantrecalibration_indel_tsfilter_level};
        }

        gatk_applyvqsr(
            {
                filehandle           => $filehandle,
                infile_path          => $applyvqsr_infile_path,
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                memory_allocation    => q{Xmx10g},
                mode                 => $mode,
                outfile_path         => $applyvqsr_outfile_path,
                recal_file_path      => $recal_file_path,
                referencefile_path   => $referencefile_path,
                temp_directory       => $temp_directory,
                tranches_file_path   => $recal_file_path . $DOT . q{tranches},
                ts_filter_level      => $ts_filter_level,
                verbosity            => $active_parameter_href->{gatk_logging_level},
            }
        );
        say {$filehandle} $NEWLINE;
    }

    ## GenotypeRefinement
    if (    $parameter_href->{cache}{trio}
        and $active_parameter_href->{gatk_calculategenotypeposteriors} )
    {

        say {$filehandle} q{## GATK CalculateGenotypePosteriors};

        my $calculategt_outfile_path =
          $outfile_path_prefix . $UNDERSCORE . q{refined} . $outfile_suffix;
        gatk_calculategenotypeposteriors(
            {
                filehandle                 => $filehandle,
                infile_path                => $outfile_path,
                java_use_large_pages       => $active_parameter_href->{java_use_large_pages},
                memory_allocation          => q{Xmx6g},
                num_ref_samples_if_no_call =>
                  $active_parameter_href->{gatk_num_reference_samples_if_no_call},
                outfile_path                 => $calculategt_outfile_path,
                pedigree                     => $commands{pedigree},
                referencefile_path           => $referencefile_path,
                supporting_callset_file_path =>
                  $active_parameter_href->{gatk_calculate_genotype_call_set},
                temp_directory => $temp_directory,
                verbosity      => $active_parameter_href->{gatk_logging_level},
            }
        );
        say {$filehandle} $NEWLINE;

        ## Change name of file to accomodate downstream
        gnu_mv(
            {
                filehandle   => $filehandle,
                infile_path  => $calculategt_outfile_path,
                outfile_path => $outfile_path,
            }
        );
        say {$filehandle} $NEWLINE;
    }

    if ( not $active_parameter_href->{gatk_variantrecalibration_keep_unnormalised} ) {

        ## BcfTools norm, Left-align and normalize indels, split multiallelics
        my $bcftools_outfile_path =
          $outfile_path_prefix . $UNDERSCORE . q{normalized} . $outfile_suffix;
        bcftools_norm(
            {
                filehandle     => $filehandle,
                infile_path    => $outfile_path,
                multiallelic   => $DASH,
                output_type    => q{v},
                outfile_path   => $bcftools_outfile_path,
                reference_path => $referencefile_path,
            }
        );
        say {$filehandle} $NEWLINE;

        ## Change name of file to accomodate downstream
        gnu_mv(
            {
                filehandle   => $filehandle,
                infile_path  => $bcftools_outfile_path,
                outfile_path => $outfile_path,
            }
        );
        say {$filehandle} $NEWLINE;
    }
    close $filehandle;

    if ( $recipe{mode} == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_case},
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
                log                               => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_ids_ref     => \@{ $active_parameter_href->{sample_ids} },
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub _build_gatk_resource_command {

## Function : Build resources in the correct format for GATK
## Returns  :
## Arguments: $resources_href => Resources to build comand for {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $resources_href;

    my $tmpl = {
        resources_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$resources_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @built_resources;

  RESOURCE:
    while ( my ( $file, $string ) = each %{$resources_href} ) {

        ## Build resource string
        push @built_resources, $string . $SPACE . $file;
    }
    return @built_resources;
}

1;
