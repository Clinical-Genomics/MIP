package MIP::Recipes::Analysis::Plink;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile splitpath};
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_plink };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $DASH       => q{-};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $PIPE       => q{|};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_plink {

## Function : Tests sample for correct relatives (only performed for samples with relatives defined in pedigree file) performed on sequence data.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => Variant call type
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infamily_directory      => In family directory
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outfamily_directory     => Out family directory
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $infamily_directory;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $outfamily_directory;
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
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        call_type => {
            default     => q{BOTH},
            strict_type => 1,
            store       => \$call_type,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infamily_directory,
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
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
        },
        outfamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfamily_directory,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Pedigree qw{ create_fam_file };
    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_family_dead_end };
    use MIP::Program::Variantcalling::Bcftools
      qw(bcftools_view bcftools_annotate);
    use MIP::Program::Variantcalling::Plink
      qw{ plink_calculate_inbreeding plink_check_sex_chroms plink_create_mibs plink_fix_fam_ped_map_freq plink_sex_check plink_variant_pruning };
    use MIP::QC::Record
      qw{ add_program_outfile_to_sample_info add_program_metafile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Program::Variantcalling::Vt qw(vt_uniq);

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Unpack parameters
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $human_genome_reference_version =
      $file_info_href->{human_genome_reference_version};
    my $human_genome_reference_source =
      $file_info_href->{human_genome_reference_source};
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    my $program_directory =
      catfile( $outaligner_dir, q{casecheck}, $program_name );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            call_type                       => $call_type,
            core_number                     => $core_number,
            directory_id                    => $family_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            process_time                    => $time,
            program_directory               => $program_directory,
            program_name                    => $program_name,
            source_environment_commands_ref => [$source_environment_cmd],
        }
    );

    # Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $program_info_file ) =
      splitpath($program_info_path);

    # To enable submission to &sample_info_qc later
    my $stderr_file = $program_info_file . q{.stderr.txt};

    # To enable submission to &sample_info_qc later
    my $stdout_file = $program_info_file . q{.stdout.txt};

    ## Paths
    my $outfamily_file_directory =
      catfile( $active_parameter_href->{outdata_dir}, $family_id );

    ## Files
    my $infile_tag =
      $file_info_href->{$family_id}{pgatk_combinevariantcallsets}{file_tag};
    my $infile_prefix = $family_id . $infile_tag . $call_type;

    ## Paths
    my $file_path_prefix = catfile( $temp_directory, $infile_prefix );

    ## Get infile_suffix from baserecalibration jobid chain
    my $infile_suffix = get_file_suffix(
        {
            jobid_chain =>
              $parameter_href->{pgatk_combinevariantcallsets}{chain},
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

    my $family_file =
      catfile( $outfamily_file_directory, $family_id . $DOT . q{fam} );

    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            fam_file_path         => $family_file,
            FILEHANDLE            => $FILEHANDLE,
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => catfile(
                $infamily_directory,
                $infile_prefix . $infile_suffix . $ASTERISK
            ),
            outfile_path => $temp_directory
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Prepare input
    say {$FILEHANDLE} q{## Remove indels using bcftools};
    bcftools_view(
        {
            exclude_types_ref => [q{indels}],
            FILEHANDLE        => $FILEHANDLE,
            infile_path       => $file_path_prefix . $infile_suffix,
            outfile_path      => $file_path_prefix
              . $UNDERSCORE
              . q{no_indels}
              . $infile_suffix,
            output_type => q{v},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Create uniq IDs and remove duplicate variants};
    bcftools_annotate(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $file_path_prefix
              . $UNDERSCORE
              . q{no_indels}
              . $infile_suffix,
            output_type    => q{v},
            remove_ids_ref => [q{ID}],
            set_id         => q?+'%CHROM:%POS:%REF:%ALT'?,
        }
    );

    print {$FILEHANDLE} $PIPE . $SPACE;

    ## Drops duplicate variants that appear later in the the VCF file
    vt_uniq(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $DASH,
            outfile_path => $file_path_prefix
              . $UNDERSCORE
              . q{no_indels_ann_uniq}
              . $infile_suffix,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    my $outfile_prefix =
      catfile( $temp_directory, $family_id . $UNDERSCORE . q{data} );

    ### Plink variant pruning and creation of unique Ids

    Readonly my $INDEP_WINDOW_SIZE   => 50;
    Readonly my $INDEP_STEP_SIZE     => 5;
    Readonly my $INDEP_VIF_THRESHOLD => 2;

    say {$FILEHANDLE} q{## Create pruning set and uniq IDs};
    plink_variant_pruning(
        {
            FILEHANDLE          => $FILEHANDLE,
            const_fid           => $family_id,
            outfile_prefix      => $outfile_prefix,
            make_bed            => 1,
            indep               => 1,
            indep_step_size     => $INDEP_STEP_SIZE,
            indep_vif_threshold => $INDEP_VIF_THRESHOLD,
            indep_window_size   => $INDEP_WINDOW_SIZE,
            set_missing_var_ids => q?@:#[?
              . $human_genome_reference_version
              . q?]\$1,\$2?,
            vcffile_path => $file_path_prefix
              . $UNDERSCORE
              . q{no_indels_ann_uniq}
              . $infile_suffix,
            vcf_half_call  => q{haploid},
            vcf_require_gt => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE}
      q{## Update Plink fam. Create ped and map file and frequency report};

    my $binary_fileset_prefix =
      catfile( $temp_directory, $family_id . $UNDERSCORE . q{data} );

    plink_fix_fam_ped_map_freq(
        {
            binary_fileset_prefix => $binary_fileset_prefix,
            FILEHANDLE            => $FILEHANDLE,
            fam_file_path         => $family_file,
            freqx                 => 1,
            make_just_fam         => 1,
            outfile_prefix        => $outfile_prefix,
            recode                => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    $outfile_prefix = catfile( $outfamily_directory, $family_id );

    # Only perform if more than 1 sample
    if ( scalar @{ $active_parameter_href->{sample_ids} } > 1 ) {

        say {$FILEHANDLE} q{## Calculate inbreeding coefficients per family};
        plink_calculate_inbreeding(
            {
                binary_fileset_prefix => $binary_fileset_prefix,
                extract_file => $binary_fileset_prefix . $DOT . q{prune.in},
                FILEHANDLE   => $FILEHANDLE,
                inbreeding_coefficients => 1,
                het                     => 1,
                outfile_prefix          => $outfile_prefix,
                small_sample            => 1,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        say {$FILEHANDLE} q{## Create Plink .mibs per family};
        plink_create_mibs(
            {
                cluster        => 1,
                FILEHANDLE     => $FILEHANDLE,
                map_file_path  => $binary_fileset_prefix . $DOT . q{map},
                matrix         => 1,
                outfile_prefix => $outfile_prefix,
                ped_file_path  => $binary_fileset_prefix . $DOT . q{ped},
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ### Plink sex-check
    ## Get parameters
    my $genome_build;

    if ( $human_genome_reference_source eq q{GRCh} ) {

        $genome_build = q{b} . $human_genome_reference_version;
    }
    else {

        $genome_build = q{hg} . $human_genome_reference_version;
    }

    $outfile_prefix =
      catfile( $temp_directory, $family_id . $UNDERSCORE . q{data} );
    Readonly my $CHR_X_NUMBER => 23;
    Readonly my $CHR_Y_NUMBER => 24;

    plink_check_sex_chroms(
        {
            binary_fileset_prefix => $binary_fileset_prefix,
            FILEHANDLE            => $FILEHANDLE,
            no_fail               => 1,
            make_bed              => 1,
            outfile_prefix        => $outfile_prefix . $UNDERSCORE . q{unsplit},
            regions_ref           => [ $CHR_X_NUMBER, $CHR_Y_NUMBER ],
            split_x               => $genome_build,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    Readonly my $FEMALE_MAX_F => 0.2;
    Readonly my $MALE_MIN_F   => 0.75;

    ## Get parameters
    my $sex_check_min_f;
    if ( $consensus_analysis_type eq q{wes} ) {
        $sex_check_min_f = $FEMALE_MAX_F . $SPACE . $MALE_MIN_F;
    }
    my $extract_file;

    if ( scalar @{ $active_parameter_href->{sample_ids} } > 1 ) {

        $extract_file = $binary_fileset_prefix . $DOT . q{prune.in};
    }

    $outfile_prefix = catfile( $outfamily_directory, $family_id );

    plink_sex_check(
        {
            binary_fileset_prefix => $binary_fileset_prefix
              . $UNDERSCORE
              . q{unsplit},
            extract_file       => $extract_file,
            FILEHANDLE         => $FILEHANDLE,
            outfile_prefix     => $outfile_prefix,
            read_freqfile_path => $binary_fileset_prefix . $DOT . q{frqx},
            sex_check_min_f    => $sex_check_min_f,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    if ( $mip_program_mode == 1 ) {

        # Only perform if more than 1 sample
        if ( scalar @{ $active_parameter_href->{sample_ids} } > 1 ) {
            ## Collect QC metadata info for later use
            add_program_outfile_to_sample_info(
                {
                    path => catfile(
                        $outfamily_directory, $family_id . $DOT . q{het}
                    ),
                    program_name     => q{inbreeding_factor},
                    sample_info_href => $sample_info_href,
                }
            );

            ## Collect QC metadata info for later use
            add_program_outfile_to_sample_info(
                {
                    path => catfile(
                        $outfamily_directory, $family_id . $DOT . q{mibs}
                    ),
                    program_name     => q{relation_check},
                    sample_info_href => $sample_info_href,
                }
            );
        }

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                path => catfile(
                    $outfamily_directory, $family_id . $DOT . q{sexcheck}
                ),
                program_name     => q{plink_sexcheck},
                sample_info_href => $sample_info_href,
            }
        );

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                path             => catfile( $directory, $stdout_file ),
                program_name     => q{plink2},
                sample_info_href => $sample_info_href,
            }
        );

        slurm_submit_job_sample_id_dependency_family_dead_end(
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

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});
    return;
}

1;
