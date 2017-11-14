package MIP::Recipes::Analysis::Plink;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile splitpath};
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;

## CPANM
use autodie qw{ :all };
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

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

## Function : ests sample for correct relatives (only performed for samples with relatives defined in pedigree file) performed on sequence data.
## Returns  :
## Arguments: $parameter_href          => Parameter hash {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $program_name            => Program name
##          : $family_id               => Family id
##          : $infamily_directory      => In family directory
##          : $outfamily_directory     => Out family directory
##          : $temp_directory          => Temporary directory
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $call_type               => Variant call type

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;
    my $infamily_directory;
    my $outfamily_directory;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
    my $temp_directory;
    my $call_type;

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
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        temp_directory => {
            required    => $arg_href->{active_parameter_href}{temp_directory},
            defined     => 1,
            strict_type => 1,
            store       => \$temp_directory,
        },
        infamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infamily_directory,
        },
        outfamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfamily_directory,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
        },
        call_type => {
            default    => q{BOTH},
            strict_type => 1,
            store       => \$call_type,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ### SET IMPORT IN ALPHABETIC ORDER
    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::Program::Variantcalling::Bcftools
      qw(bcftools_view bcftools_annotate);
    use MIP::Program::Variantcalling::Plink qw{ plink_calculate_inbreeding plink_check_sex_chroms plink_create_mibs plink_fix_fam_ped_map_freq plink_sex_check plink_variant_pruning };
    use Program::Variantcalling::Vt qw(vt_uniq);
    use MIP::QC::Record
      qw{ add_program_outfile_to_sample_info add_program_metafile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Unpack parameters
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    my $program_directory = catfile(lc($outaligner_dir), q{casecheck}, $program_name );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $family_id,
            program_name          => $program_name,
            program_directory     => $program_directory,
            call_type   => $call_type,
            core_number => $core_number,
            process_time => $time,
        }
    );

    #Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $program_info_file ) = splitpath($program_info_path);

    #To enable submission to &sample_info_qc later
    my $stderr_file = $program_info_file . q{.stderr.txt};

    #To enable submission to &sample_info_qc later
    my $stdout_file = $program_info_file . q{.stdout.txt};

    ## Paths
    my $outfamily_file_directory =
      catfile( $active_parameter_href->{outdata_dir}, $family_id );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$family_id}{pgatk_combinevariantcallsets}
      {file_tag};
    my $infile_prefix = $family_id . $infile_tag . $call_type;
    my $file_path_prefix = catfile( $temp_directory, $infile_prefix );

    ## Get infile_suffix from baserecalibration jobid chain
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
            jobid_chain =>
              $parameter_href->{pgatk_combinevariantcallsets}{chain},
        }
    );

    my $family_file =
      catfile( $outfamily_file_directory, $family_id . q{.fam} );

    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            parameter_href        => $parameter_href,
            active_parameter_href => $active_parameter_href,
            sample_info_href      => $sample_info_href,
            FILEHANDLE            => $FILEHANDLE,
            fam_file_path         => $family_file,
        }
    );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => catfile(
                $infamily_directory, $infile_prefix . $infile_suffix . $ASTERISK
            ),
            outfile_path => $temp_directory
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Prepare input
    say {$FILEHANDLE} q{## Remove indels using bcftools};
    bcftools_view(
        {
            infile_path  => $file_path_prefix . $infile_suffix,
            outfile_path => $file_path_prefix . $UNDERSCORE . q{no_indels} . $infile_suffix,
            output_type  => q{v},
            exclude_types_ref => [q{indels}],
            FILEHANDLE        => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Create uniq IDs and remove duplicate variants};
    bcftools_annotate(
        {
            remove_ids_ref => [q{ID}],
            set_id         => q?+'%CHROM:%POS:%REF:%ALT'?,
            infile_path    => $file_path_prefix . $UNDERSCORE . q{no_indels} . $infile_suffix,
            output_type    => q{v},
            FILEHANDLE     => $FILEHANDLE,
        }
    );

    print {$FILEHANDLE} $PIPE . $SPACE ;

    ## Drops duplicate variants that appear later in the the VCF file
    vt_uniq(
        {
            infile_path  => $DASH,
            outfile_path => $file_path_prefix
              . $UNDERSCORE . q{no_indels_ann_uniq}
              . $infile_suffix,
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ### Plink variant pruning and creation of unique Ids
    say {$FILEHANDLE} q{## Create pruning set and uniq IDs};
    plink_variant_pruning(
        {
            vcffile_path => $file_path_prefix
              . $UNDERSCORE . q{no_indels_ann_uniq}
              . $infile_suffix,
            outfile_prefix =>
              catfile( $temp_directory, $family_id . $UNDERSCORE . q{data} ),
            vcf_require_gt      => 1,
            vcf_half_call       => q{haploid},
            set_missing_var_ids => q?@:#[?
              . $file_info_href->{human_genome_reference_version}
              . q?]\$1,\$2?,
            const_fid           => $family_id,
            make_bed            => 1,
            indep               => 1,
            indep_window_size   => 50,
            indep_step_size     => 5,
            indep_vif_threshold => 2,
            FILEHANDLE          => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE}
      q{## Update Plink fam. Create ped and map file and frequency report};

    my $binary_fileset_prefix =  catfile( $temp_directory, $family_id . $UNDERSCORE . q{data} );
    my $outfile_prefix = catfile( $temp_directory, $family_id . $UNDERSCORE . q{data} );

    plink_fix_fam_ped_map_freq(
        {
            binary_fileset_prefix => $binary_fileset_prefix ,
            fam_file_path => $family_file,
            make_just_fam => 1,
            recode        => 1,
            freqx         => 1,
            outfile_prefix => $outfile_prefix,
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    $outfile_prefix = catfile( $outfamily_directory, $family_id );

    # Only perform if more than 1 sample
    if ( scalar( @{ $active_parameter_href->{sample_ids} } ) > 1 ) {

        say {$FILEHANDLE} "## Calculate inbreeding coefficients per family";
        plink_calculate_inbreeding(
            {
                binary_fileset_prefix => $binary_fileset_prefix,
                outfile_prefix => $outfile_prefix,
                het                     => 1,
                small_sample            => 1,
                inbreeding_coefficients => 1,
                extract_file            => $binary_fileset_prefix . $DOT . q{prune.in},
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        say {$FILEHANDLE} q{## Create Plink .mibs per family};
        plink_create_mibs(
            {
                ped_file_path => $binary_fileset_prefix . $DOT . q{ped},
                map_file_path => $binary_fileset_prefix . $DOT . q{map},
                cluster => 1,
                matrix  => 1,
                outfile_prefix => $outfile_prefix,
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ### Plink sex-check
    ## Get parameters
    my $genome_build;
    if ( $file_info_href->{human_genome_reference_source} eq q{GRCh} ) {

        $genome_build = q{b} . $file_info_href->{human_genome_reference_version};
    }
    else {

        $genome_build =
          q{hg} . $file_info_href->{human_genome_reference_version};
    }



    plink_check_sex_chroms(
        {
            regions_ref => [ 23, 24 ],
            split_x     => $genome_build,
            no_fail     => 1,
            make_bed    => 1,
            binary_fileset_prefix => $binary_fileset_prefix,
            outfile_prefix => $outfile_prefix . $UNDERSCORE . q{unsplit},
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Get parameters
    # Exome/rapid analysis use combined reference for more power
    my $sex_check_min_F;
    if (   ( $consensus_analysis_type eq q{wes} ) || ( $consensus_analysis_type eq q{rapid} ) )
    {
        $sex_check_min_F = q{0.2 0.75};
    }
    my $extract_file;

    if ( scalar( @{ $active_parameter_href->{sample_ids} } ) > 1 ) {

        $extract_file = $binary_fileset_prefix . $DOT . q{prune.in};
    }

    plink_sex_check(
        {
            sex_check_min_f => $sex_check_min_F,
            extract_file    => $extract_file,
            read_freqfile_path => $binary_fileset_prefix . $DOT . q{frqx},
            binary_fileset_prefix => $binary_fileset_prefix . $UNDERSCORE . q{unsplit},
            outfile_prefix => $outfile_prefix,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    if ( $mip_program_mode == 1 ) {

        if ( scalar( @{ $active_parameter_href->{sample_ids} } ) > 1 )
        {    #Only perform if more than 1 sample

            ## Collect QC metadata info for later use
            add_program_outfile_to_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    program_name     => q{inbreeding_factor},
                    outdirectory     => $outfamily_directory,
                    outfile          => $family_id . $DOT . q{het},
                }
            );

            ## Collect QC metadata info for later use
            add_program_outfile_to_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    program_name     => q{relation_check},
                    outdirectory     => $outfamily_directory,
                    outfile          => $family_id . $DOT . q{mibs},
                }
            );
        }

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => q{plink_sexcheck},
                outdirectory     => $outfamily_directory,
                outfile          => $family_id . $DOT . q{sexcheck},
            }
        );

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => q{plink2},
                outdirectory     => $directory,
                outfile          => $stdout_file,
            }
        );

        slurm_submit_job_sample_id_dependency_family_dead_end(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $family_id,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );

    }

    close $FILEHANDLE;
    return;
}

1;
