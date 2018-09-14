package MIP::Recipes::Analysis::Gatk_variantrecalibration;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw { uniq };
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ analysis_gatk_variantrecalibration_wgs analysis_gatk_variantrecalibration_wes };

}

## Constants
Readonly my $ASTERIX    => q{*};
Readonly my $AMPERSAND  => q{&};
Readonly my $COLON      => q{:};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_gatk_variantrecalibration_wgs {

## Function : GATK VariantRecalibrator/ApplyRecalibration analysis recipe for wgs data
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => Variant call type
##          : $family_id               => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $sample_info_href;

    ## Default(s)
    my $call_type;
    my $family_id;
    my $outaligner_dir;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        call_type =>
          { default => q{BOTH}, store => \$call_type, strict_type => 1, },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
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
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
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

    use MIP::Delete::List qw{ delete_contig_elements };
    use MIP::File::Format::Pedigree qw{ create_fam_file gatk_pedigree_flag };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Bcftools qw{ bcftools_norm };
    use MIP::Program::Variantcalling::Gatk
      qw{ gatk_variantrecalibrator gatk_applyvqsr gatk_selectvariants gatk_calculategenotypeposteriors };
    use MIP::QC::Record
      qw{ add_program_outfile_to_sample_info add_processing_metafile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };

    ## Constants
    Readonly my $MAX_GAUSSIAN_LEVEL_INDEL             => 4;
    Readonly my $MAX_GAUSSIAN_LEVEL_SNV               => 6;
    Readonly my $MAX_GAUSSIAN_LEVEL_SNV_SINGLE_SAMPLE => 4;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $program_mode = $active_parameter_href->{$program_name};

    ## Alias
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $enable_indel_max_gaussians_filter =
      $active_parameter_href->{gatk_variantrecalibration_indel_max_gaussians};
    my $enable_snv_max_gaussians_filter =
      $active_parameter_href->{gatk_variantrecalibration_snv_max_gaussians};
    my $gatk_jar =
      catfile( $active_parameter_href->{gatk_path}, q{GenomeAnalysisTK.jar} );
    my $job_id_chain       = $parameter_href->{$program_name}{chain};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $resource_indel_href =
      $active_parameter_href->{gatk_variantrecalibration_resource_indel};
    my $resource_snv_href =
      $active_parameter_href->{gatk_variantrecalibration_resource_snv};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Assign directories
    my $program_outdirectory_name =
      $parameter_href->{$program_name}{outdir_name};

    ## For ".fam" file
    my $infamily_directory = catfile( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir, $program_outdirectory_name );
    my $outfamily_directory = $infamily_directory;
    my $outfamily_file_directory =
      catdir( $active_parameter_href->{outdata_dir}, $family_id );

    ## Used downstream
    $parameter_href->{$program_name}{indirectory} = $outfamily_directory;

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$family_id}{gatk_genotypegvcfs}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$program_name}{file_tag};

    ## Files
    my $infile_prefix  = $family_id . $infile_tag . $call_type;
    my $outfile_prefix = $family_id . $outfile_tag . $call_type;

    ## Paths
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ### Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            program_name   => q{gatk_genotypegvcfs},
            suffix_key     => q{outfile_suffix},
        }
    );

    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            file_suffix    => $parameter_href->{$program_name}{outfile_suffix},
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            call_type             => $call_type,
            core_number           => $core_number,
            directory_id          => $family_id,
            FILEHANDLE            => $FILEHANDLE,
            job_id_href           => $job_id_href,
            log                   => $log,
            process_time          => $time,
            program_directory =>
              catfile( $outaligner_dir, $program_outdirectory_name ),
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ## Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) =
      splitpath( $program_info_path . $DOT . q{stderr.txt} );

    ## Create .fam file to be used in variant calling analyses
    my $fam_file_path =
      catfile( $outfamily_file_directory, $family_id . $DOT . q{fam} );
    create_fam_file(
        {
            active_parameter_href => $active_parameter_href,
            fam_file_path         => $fam_file_path,
            FILEHANDLE            => $FILEHANDLE,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Check if "--pedigree" should be included in analysis
    my %commands = gatk_pedigree_flag(
        {
            fam_file_path => $fam_file_path,
            program_name  => $program_name,
        }
    );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => catfile(
                $infamily_directory, $infile_prefix . $infile_suffix . $ASTERIX
            ),
            outfile_path => $temp_directory
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ### GATK VariantRecalibrator
    ## Set mode to be used in variant recalibration
# SNP and INDEL will be recalibrated successively in the same file because when you specify eg SNP mode, the indels are emitted without modification, and vice-versa.
    my @modes = qw{ SNP INDEL };

  MODE:
    foreach my $mode (@modes) {

        say {$FILEHANDLE} q{## GATK VariantRecalibrator};

        ##Get parameters
        my $infile_path;
        my $max_gaussian_level;

        if ( $mode eq q{SNP} ) {

            $infile_path = $file_path_prefix . $infile_suffix;

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

            $infile_path =
              $outfile_path_prefix . $DOT . q{SNV} . $infile_suffix;

            ## Use hard filtering
            if ($enable_indel_max_gaussians_filter) {

                $max_gaussian_level = $MAX_GAUSSIAN_LEVEL_INDEL;
            }
        }

        my @annotations =
          @{ $active_parameter_href->{gatk_variantrecalibration_annotations} };
        if ( $consensus_analysis_type ne q{wgs} ) {

            ### Special case: Not to be used with hybrid capture. NOTE: Disabled when analysing wes + wgs in the same run
            ## Removes an element from array and return new array while leaving orginal elements_ref untouched
            @annotations = delete_contig_elements(
                {
                    elements_ref       => \@annotations,
                    remove_contigs_ref => [qw{ DP }],
                }
            );
        }

        my @resources;
        if ( $mode eq q{SNP} ) {

            @resources = _build_gatk_resource_command(
                { resources_href => $resource_snv_href, } );
        }
        if ( $mode eq q{INDEL} ) {

            @resources = _build_gatk_resource_command(
                { resources_href => $resource_indel_href, } );
        }

        my $recal_file_path = $file_path_prefix . $DOT . q{intervals};
        gatk_variantrecalibrator(
            {
                annotations_ref => \@annotations,
                FILEHANDLE      => $FILEHANDLE,
                infile_path     => $infile_path,
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                max_gaussian_level    => $max_gaussian_level,
                memory_allocation     => q{Xmx24g},
                mode                  => $mode,
                outfile_path          => $recal_file_path,
                referencefile_path    => $referencefile_path,
                resources_ref         => \@resources,
                rscript_file_path     => $recal_file_path . $DOT . q{plots.R},
                temp_directory        => $temp_directory,
                tranches_file_path    => $recal_file_path . $DOT . q{tranches},
                trust_all_polymorphic => $active_parameter_href
                  ->{gatk_variantrecalibration_trust_all_polymorphic},
                verbosity => $active_parameter_href->{gatk_logging_level},
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## GATK ApplyVQSR
        say {$FILEHANDLE} q{## GATK ApplyVQSR};

        ## Get parameters
        my $outfile_path;
        my $ts_filter_level;

        if ( $mode eq q{SNP} ) {

            $infile_path = $file_path_prefix . $infile_suffix;
            $outfile_path =
              $outfile_path_prefix . $DOT . q{SNV} . $outfile_suffix;
            $ts_filter_level = $active_parameter_href
              ->{gatk_variantrecalibration_snv_tsfilter_level};
        }

        ## Use created recalibrated snp vcf as input
        if ( $mode eq q{INDEL} ) {

            $infile_path =
              $outfile_path_prefix . $DOT . q{SNV} . $outfile_suffix;
            $outfile_path    = $outfile_path_prefix . $outfile_suffix;
            $ts_filter_level = $active_parameter_href
              ->{gatk_variantrecalibration_indel_tsfilter_level};
        }

        gatk_applyvqsr(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $infile_path,
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                memory_allocation  => q{Xmx10g},
                mode               => $mode,
                outfile_path       => $outfile_path,
                recal_file_path    => $recal_file_path,
                referencefile_path => $referencefile_path,
                temp_directory     => $temp_directory,
                tranches_file_path => $recal_file_path . $DOT . q{tranches},
                ts_filter_level    => $ts_filter_level,
                verbosity => $active_parameter_href->{gatk_logging_level},
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## GenotypeRefinement
    if ( $parameter_href->{dynamic_parameter}{trio} ) {

        say {$FILEHANDLE} q{## GATK CalculateGenotypePosteriors};

        my $outfile_path =
          $outfile_path_prefix . $UNDERSCORE . q{refined} . $outfile_suffix;
        gatk_calculategenotypeposteriors(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $outfile_path_prefix . $outfile_suffix,
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                memory_allocation            => q{Xmx6g},
                outfile_path                 => $outfile_path,
                pedigree                     => $commands{pedigree},
                referencefile_path           => $referencefile_path,
                supporting_callset_file_path => $active_parameter_href
                  ->{gatk_calculategenotypeposteriors_support_set},
                temp_directory => $temp_directory,
                verbosity      => $active_parameter_href->{gatk_logging_level},
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Change name of file to accomodate downstream
        gnu_mv(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $outfile_path,
                outfile_path => $outfile_path_prefix . $outfile_suffix,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## BcfTools norm, Left-align and normalize indels, split multiallelics
    my $bcftools_outfile_path =
      $outfile_path_prefix . $UNDERSCORE . q{normalized} . $outfile_suffix;
    bcftools_norm(
        {
            FILEHANDLE     => $FILEHANDLE,
            infile_path    => $outfile_path_prefix . $outfile_suffix,
            multiallelic   => q{-},
            output_type    => q{v},
            outfile_path   => $bcftools_outfile_path,
            reference_path => $referencefile_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Change name of file to accomodate downstream
    gnu_mv(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $bcftools_outfile_path,
            outfile_path => $outfile_path_prefix . $outfile_suffix,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    my @outfiles = (
        $outfile_path_prefix . $outfile_suffix . $ASTERIX,
        $file_path_prefix . $DOT . q{intervals} . $DOT . q{tranches.pdf},
    );

  FILE:
    foreach my $outfile (@outfiles) {

        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $outfile,
                outfile_path => $outfamily_directory,
            }
        );
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        ## Collect QC metadata info for later use
        my $program_outfile_path =
          catfile( $outfamily_directory, $outfile_prefix . $outfile_suffix );
        add_program_outfile_to_sample_info(
            {
                path             => $program_outfile_path,
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );

        # Used to find order of samples in qccollect downstream
        add_program_outfile_to_sample_info(
            {
                path             => $program_outfile_path,
                program_name     => q{pedigree_check},
                sample_info_href => $sample_info_href,
            }
        );

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                sbatch_file_name => $file_path,
            }
        );
    }
    return;
}

sub analysis_gatk_variantrecalibration_wes {

## Function : GATK VariantRecalibrator/ApplyRecalibration analysis recipe for wes data
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => Variant call type
##          : $family_id               => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $call_type;
    my $outaligner_dir;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        call_type =>
          { default => q{BOTH}, store => \$call_type, strict_type => 1, },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
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
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
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

    use MIP::Delete::List qw{ delete_contig_elements };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Bcftools qw{ bcftools_norm };
    use MIP::Program::Variantcalling::Gatk
      qw{ gatk_variantrecalibrator gatk_applyvqsr gatk_selectvariants gatk_calculategenotypeposteriors };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };

    ## Constants
    Readonly my $MAX_GAUSSIAN_LEVEL => 4;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set program mode
    my $program_mode = $active_parameter_href->{$program_name};

    ## Alias
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $enable_indel_max_gaussians_filter =
      $active_parameter_href->{gatk_variantrecalibration_indel_max_gaussians};
    my $enable_snv_max_gaussians_filter =
      $active_parameter_href->{gatk_variantrecalibration_snv_max_gaussians};
    my $gatk_jar =
      catfile( $active_parameter_href->{gatk_path}, q{GenomeAnalysisTK.jar} );
    my $job_id_chain       = $parameter_href->{$program_name}{chain};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $resource_indel_href =
      $active_parameter_href->{gatk_variantrecalibration_resource_indel};
    my $resource_snv_href =
      $active_parameter_href->{gatk_variantrecalibration_resource_snv};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Assign directories
    my $program_outdirectory_name =
      $parameter_href->{$program_name}{outdir_name};

    my $infamily_directory = catfile( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir, $program_outdirectory_name );
    my $outfamily_directory = $infamily_directory;
    ## For ".fam" file
    my $outfamily_file_directory =
      catdir( $active_parameter_href->{outdata_dir}, $family_id );

    ## Used downstream
    $parameter_href->{$program_name}{indirectory} = $outfamily_directory;

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$family_id}{gatk_genotypegvcfs}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$program_name}{file_tag};

    ## Files
    my $infile_prefix  = $family_id . $infile_tag . $call_type;
    my $outfile_prefix = $family_id . $outfile_tag . $call_type;

    ## Paths
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ### Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            program_name   => q{gatk_genotypegvcfs},
            suffix_key     => q{outfile_suffix},
        }
    );

    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            file_suffix    => $parameter_href->{$program_name}{outfile_suffix},
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            call_type             => $call_type,
            core_number           => $core_number,
            directory_id          => $family_id,
            FILEHANDLE            => $FILEHANDLE,
            job_id_href           => $job_id_href,
            log                   => $log,
            process_time          => $time,
            program_directory =>
              catfile( $outaligner_dir, $program_outdirectory_name ),
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ## Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) =
      splitpath( $program_info_path . $DOT . q{stderr.txt} );

    ## Create .fam file to be used in variant calling analyses
    my $fam_file_path =
      catfile( $outfamily_file_directory, $family_id . $DOT . q{fam} );
    create_fam_file(
        {
            active_parameter_href => $active_parameter_href,
            fam_file_path         => $fam_file_path,
            FILEHANDLE            => $FILEHANDLE,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Check if "--pedigree" and "--pedigreeValidationType" should be included in analysis
    my %commands = gatk_pedigree_flag(
        {
            fam_file_path => $fam_file_path,
            program_name  => $program_name,
        }
    );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => catfile(
                $infamily_directory, $infile_prefix . $infile_suffix . $ASTERIX
            ),
            outfile_path => $temp_directory
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ### GATK VariantRecalibrator
    ## Set mode to be used in variant recalibration
# Exome will be processed using mode BOTH since there are to few INDELS to use in the recalibration model even though using 30 exome BAMS in Haplotypecaller step.
    my @modes = q{BOTH};

  MODE:
    foreach my $mode (@modes) {

        say {$FILEHANDLE} q{## GATK VariantRecalibrator};

        ##Get parameters
        my $max_gaussian_level;
        my $infile_path;

        ## Exome analysis use combined reference for more power
# Infile HaplotypeCaller combined vcf which used reference gVCFs to create combined vcf (30> samples gCVFs)
        $infile_path = $file_path_prefix . $infile_suffix;

        my @annotations =
          @{ $active_parameter_href->{gatk_variantrecalibration_annotations} };
        ### Special case: Not to be used with hybrid capture. NOTE: Disabled when analysing wes + wgs in the same run
        ## Removes an element from array and return new array while leaving orginal elements_ref untouched
        @annotations = delete_contig_elements(
            {
                elements_ref       => \@annotations,
                remove_contigs_ref => [qw{ DP }],
            }
        );
        my @snv_resources = _build_gatk_resource_command(
            { resources_href => $resource_snv_href, } );
        my @indel_resources = _build_gatk_resource_command(
            { resources_href => $resource_indel_href, } );

        # Create distinct set i.e. no duplicates.
        my @resources = uniq( @indel_resources, @snv_resources );

        ## Use hard filtering
        if (   $enable_snv_max_gaussians_filter
            || $enable_indel_max_gaussians_filter )
        {

            $max_gaussian_level = $MAX_GAUSSIAN_LEVEL;
        }

        my $recal_file_path = $file_path_prefix . $DOT . q{intervals};
        gatk_variantrecalibrator(
            {
                annotations_ref => \@annotations,
                FILEHANDLE      => $FILEHANDLE,
                infile_path     => $infile_path,
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                max_gaussian_level => $max_gaussian_level,
                memory_allocation  => q{Xmx10g},
                mode               => $mode,
                outfile_path       => $recal_file_path,
                referencefile_path => $referencefile_path,
                resources_ref      => \@resources,
                rscript_file_path  => $recal_file_path . $DOT . q{plots.R},
                temp_directory     => $temp_directory,
                tranches_file_path => $recal_file_path . $DOT . q{tranches},
                verbosity => $active_parameter_href->{gatk_logging_level},
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## GATK ApplyVQSR
        say {$FILEHANDLE} q{## GATK ApplyVQSR};

        ## Get parameters
        my $outfile_path;
        my $ts_filter_level;
        ## Exome analysis use combined reference for more power

        ## Infile genotypegvcfs combined vcf which used reference gVCFs to create combined vcf file
        $infile_path = $file_path_prefix . $infile_suffix;
        $outfile_path =
          $outfile_path_prefix . $UNDERSCORE . q{filtered} . $outfile_suffix;
        $ts_filter_level = $active_parameter_href
          ->{gatk_variantrecalibration_snv_tsfilter_level};

        gatk_applyvqsr(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $infile_path,
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                memory_allocation  => q{Xmx10g},
                mode               => $mode,
                outfile_path       => $outfile_path,
                recal_file_path    => $recal_file_path,
                referencefile_path => $referencefile_path,
                temp_directory     => $temp_directory,
                tranches_file_path => $recal_file_path . $DOT . q{tranches},
                ts_filter_level    => $ts_filter_level,
                verbosity => $active_parameter_href->{gatk_logging_level},
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## BcfTools norm, Left-align and normalize indels, split multiallelics
    bcftools_norm(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $outfile_path_prefix
              . $UNDERSCORE
              . q{filtered}
              . $outfile_suffix,
            multiallelic => q{-},
            outfile_path => $outfile_path_prefix
              . $UNDERSCORE
              . q{filtered_normalized}
              . $outfile_suffix,
            output_type     => q{v},
            reference_path  => $referencefile_path,
            stderrfile_path => $outfile_path_prefix
              . $UNDERSCORE
              . q{filtered_normalized.stderr},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ### GATK SelectVariants

    ## Removes all genotype information for exome ref and recalulates meta-data info for remaining samples in new file.
    # Exome analysis

    say {$FILEHANDLE} q{## GATK SelectVariants};

    gatk_selectvariants(
        {
            FILEHANDLE           => $FILEHANDLE,
            exclude_non_variants => 1,
            infile_path          => $outfile_path_prefix
              . $UNDERSCORE
              . q{filtered_normalized}
              . $outfile_suffix,
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            memory_allocation  => q{Xmx2g},
            outfile_path       => $outfile_path_prefix . $outfile_suffix,
            referencefile_path => $referencefile_path,
            sample_names_ref   => \@{ $active_parameter_href->{sample_ids} },
            temp_directory     => $temp_directory,
            verbosity          => $active_parameter_href->{gatk_logging_level},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## GenotypeRefinement
    if ( $parameter_href->{dynamic_parameter}{trio} ) {

        say {$FILEHANDLE} q{## GATK CalculateGenotypePosteriors};

        my $outfile_path =
          $outfile_path_prefix . $UNDERSCORE . q{refined} . $outfile_suffix;
        gatk_calculategenotypeposteriors(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $outfile_path_prefix . $outfile_suffix,
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                memory_allocation            => q{Xmx6g},
                outfile_path                 => $outfile_path,
                pedigree                     => $commands{pedigree},
                referencefile_path           => $referencefile_path,
                supporting_callset_file_path => $active_parameter_href
                  ->{gatk_calculategenotypeposteriors_support_set},
                temp_directory => $temp_directory,
                verbosity      => $active_parameter_href->{gatk_logging_level},
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Change name of file to accomodate downstream
        gnu_mv(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $outfile_path,
                outfile_path => $outfile_path_prefix . $outfile_suffix,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## BcfTools norm, Left-align and normalize indels, split multiallelics
    bcftools_norm(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $outfile_path_prefix . $outfile_suffix,
            multiallelic => q{-},
            outfile_path => $outfile_path_prefix
              . $UNDERSCORE
              . q{normalized}
              . $outfile_suffix,
            output_type    => q{v},
            reference_path => $referencefile_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Change name of file to accomodate downstream
    gnu_mv(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $outfile_path_prefix
              . $UNDERSCORE
              . q{normalized}
              . $outfile_suffix,
            outfile_path => $outfile_path_prefix . $outfile_suffix,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    my @outfiles = (
        $outfile_path_prefix . $outfile_suffix . $ASTERIX,
        $file_path_prefix . $DOT . q{intervals} . $DOT . q{tranches.pdf},
    );

  FILE:
    foreach my $outfile (@outfiles) {

        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $outfile,
                outfile_path => $outfamily_directory,
            }
        );
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        ## Collect QC metadata info for later use
        my $program_outfile_path =
          catfile( $outfamily_directory, $outfile_prefix . $outfile_suffix );
        add_program_outfile_to_sample_info(
            {
                path             => $program_outfile_path,
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );

        # Used to find order of samples in qccollect downstream
        add_program_outfile_to_sample_info(
            {
                path             => $program_outfile_path,
                program_name     => q{pedigree_check},
                sample_info_href => $sample_info_href,
            }
        );

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                sbatch_file_name => $file_path,
            }
        );
    }
    return;
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
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$resources_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @built_resources;

    while ( my ( $file, $string ) = each %{$resources_href} ) {

        push @built_resources, $string . $COLON . $file;
    }
    return @built_resources;
}

1;

