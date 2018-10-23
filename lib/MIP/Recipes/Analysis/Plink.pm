package MIP::Recipes::Analysis::Plink;

use 5.026;
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
    our $VERSION = 1.04;

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
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
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
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
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

    use MIP::File::Format::Pedigree qw{ create_fam_file };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_module_parameters get_program_attributes };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Variantcalling::Bcftools qw(bcftools_view bcftools_annotate);
    use MIP::Program::Variantcalling::Plink
      qw{ plink_calculate_inbreeding plink_check_sex_chroms plink_create_mibs plink_fix_fam_ped_map_freq plink_sex_check plink_variant_pruning };
    use MIP::Program::Variantcalling::Vt qw(vt_uniq);
    use MIP::QC::Record
      qw{ add_program_outfile_to_sample_info add_program_metafile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constants
    Readonly my $CHR_X_NUMBER        => 23;
    Readonly my $CHR_Y_NUMBER        => 24;
    Readonly my $FEMALE_MAX_F        => 0.2;
    Readonly my $INDEP_WINDOW_SIZE   => 50;
    Readonly my $INDEP_STEP_SIZE     => 5;
    Readonly my $INDEP_VIF_THRESHOLD => 2;
    Readonly my $MALE_MIN_F          => 0.75;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $family_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            program_name   => $program_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path_prefix = $io{in}{file_path_prefix};
    my $infile_suffix      = $io{in}{file_suffix};
    my $infile_path        = $io{in}{file_path};

    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $human_genome_reference_version =
      $file_info_href->{human_genome_reference_version};
    my $human_genome_reference_source = $file_info_href->{human_genome_reference_source};
    my $job_id_chain                  = get_program_attributes(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
            attribute      => q{chain},
        }
    );
    my $program_mode = $active_parameter_href->{$program_name};
    my @sample_ids   = @{ $active_parameter_href->{sample_ids} };
    my ( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
    );

    ## Set outfiles depending on sample data
    my %plink_outfiles = (
        check_for_sex    => { plink_sexcheck => q{sexcheck}, },
        multiple_samples => {
            inbreeding_factor => q{het},
            relation_check    => q{mibs},
        },
    );

    my @plink_outfiles;
  MODE:
    while ( my ( $mode, $program_href ) = each %plink_outfiles ) {

      PLINK_PROGRAM:
        while ( my ( $file_name_prefix, $file_suffix ) = each %{$program_href} ) {

            if ( scalar @sample_ids > 1 and $mode eq q{multiple_samples} ) {

                push @plink_outfiles, $file_name_prefix . $DOT . $file_suffix;
                next;
            }
            if (    $active_parameter_href->{found_other_count} ne scalar @sample_ids
                and $mode eq q{check_for_sex} )
            {
                push @plink_outfiles, $file_name_prefix . $DOT . $file_suffix;
            }
        }
    }

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $family_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => \@plink_outfiles,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                program_name     => $program_name,
                temp_directory   => $temp_directory,
            }
        )
    );

    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my %outfile_path        = %{ $io{out}{file_path_href} };

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $recipe_file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $family_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            program_directory               => $program_name,
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
        }
    );

    # Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $program_info_file ) = splitpath($program_info_path);

    # To enable submission to %sample_info_qc later
    my $stdout_file_path = catfile( $directory, $program_info_file . q{.stdout.txt} );

    ### SHELL:

    my $plink_outfile_prefix =
      catfile( $outdir_path_prefix, $family_id . $UNDERSCORE . q{data} );
    my $family_file_path = catfile( $outdir_path_prefix, $family_id . $DOT . q{fam} );

    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            active_parameter_href => $active_parameter_href,
            fam_file_path         => $family_file_path,
            FILEHANDLE            => $FILEHANDLE,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Prepare input
    say {$FILEHANDLE} q{## Remove indels using bcftools};
    my $view_outfile_path =
      $infile_path_prefix . $UNDERSCORE . q{no_indels} . $infile_suffix;
    bcftools_view(
        {
            exclude_types_ref => [qw{ indels }],
            FILEHANDLE        => $FILEHANDLE,
            infile_path       => $infile_path,
            outfile_path      => $view_outfile_path,
            output_type       => q{v},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Create uniq IDs and remove duplicate variants};
    bcftools_annotate(
        {
            FILEHANDLE     => $FILEHANDLE,
            infile_path    => $view_outfile_path,
            output_type    => q{v},
            remove_ids_ref => [q{ID}],
            set_id         => q?+'%CHROM:%POS:%REF:%ALT'?,
        }
    );

    print {$FILEHANDLE} $PIPE . $SPACE;

    ## Drops duplicate variants that appear later in the the VCF file
    my $uniq_outfile_path =
      $infile_path_prefix . $UNDERSCORE . q{no_indels_ann_uniq} . $infile_suffix;
    vt_uniq(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $DASH,
            outfile_path => $uniq_outfile_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ### Plink variant pruning and creation of unique Ids
    say {$FILEHANDLE} q{## Create pruning set and uniq IDs};
    plink_variant_pruning(
        {
            const_fid           => $family_id,
            FILEHANDLE          => $FILEHANDLE,
            make_bed            => 1,
            outfile_prefix      => $plink_outfile_prefix,
            indep               => 1,
            indep_step_size     => $INDEP_STEP_SIZE,
            indep_vif_threshold => $INDEP_VIF_THRESHOLD,
            indep_window_size   => $INDEP_WINDOW_SIZE,
            set_missing_var_ids => q?@:#[?
              . $human_genome_reference_version
              . q?]\$1,\$2?,
            vcffile_path   => $uniq_outfile_path,
            vcf_half_call  => q{haploid},
            vcf_require_gt => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE}
      q{## Update Plink fam. Create ped and map file and frequency report};
    ## Get parameters
    my $allow_no_sex;

    ## If not all samples have a known sex
    if ( $active_parameter_href->{found_other_count} ) {

        $allow_no_sex = 1;
    }
    my $binary_fileset_prefix = $plink_outfile_prefix;

    plink_fix_fam_ped_map_freq(
        {
            allow_no_sex          => $allow_no_sex,
            binary_fileset_prefix => $binary_fileset_prefix,
            fam_file_path         => $family_file_path,
            FILEHANDLE            => $FILEHANDLE,
            freqx                 => 1,
            make_just_fam         => 1,
            outfile_prefix        => $plink_outfile_prefix,
            recode                => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    my $inbreeding_outfile_prefix = catfile( $outdir_path_prefix, $family_id );

    # Only perform if more than 1 sample
    if ( scalar @sample_ids > 1 ) {

        say {$FILEHANDLE} q{## Calculate inbreeding coefficients per family};
        plink_calculate_inbreeding(
            {
                binary_fileset_prefix   => $binary_fileset_prefix,
                extract_file            => $binary_fileset_prefix . $DOT . q{prune.in},
                FILEHANDLE              => $FILEHANDLE,
                inbreeding_coefficients => 1,
                het                     => 1,
                outfile_prefix          => $inbreeding_outfile_prefix,
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
                outfile_prefix => $inbreeding_outfile_prefix,
                ped_file_path  => $binary_fileset_prefix . $DOT . q{ped},
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    if ( $active_parameter_href->{found_other_count} ne scalar @sample_ids ) {

        ## Only if not all samples have unknown sex
        ### Plink sex-check
        ## Get parameters
        my $genome_build;

        ## Set correct build prefix
        if ( $human_genome_reference_source eq q{GRCh} ) {

            $genome_build = q{b} . $human_genome_reference_version;
        }
        else {

            $genome_build = q{hg} . $human_genome_reference_version;
        }

        plink_check_sex_chroms(
            {
                binary_fileset_prefix => $binary_fileset_prefix,
                FILEHANDLE            => $FILEHANDLE,
                no_fail               => 1,
                make_bed              => 1,
                outfile_prefix        => $plink_outfile_prefix . $UNDERSCORE . q{unsplit},
                regions_ref           => [ $CHR_X_NUMBER, $CHR_Y_NUMBER ],
                split_x               => $genome_build,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Get parameters
        my $sex_check_min_f;
        if ( $consensus_analysis_type eq q{wes} ) {

            $sex_check_min_f = $FEMALE_MAX_F . $SPACE . $MALE_MIN_F;
        }
        my $extract_file;
        my $read_freqfile_path;

        if ( scalar @sample_ids > 1 ) {

            $extract_file       = $binary_fileset_prefix . $DOT . q{prune.in};
            $read_freqfile_path = $binary_fileset_prefix . $DOT . q{frqx};
        }

        my $sex_check_outfile_prefix = catfile( $outdir_path_prefix, $family_id );

        plink_sex_check(
            {
                binary_fileset_prefix => $binary_fileset_prefix
                  . $UNDERSCORE
                  . q{unsplit},
                extract_file       => $extract_file,
                FILEHANDLE         => $FILEHANDLE,
                outfile_prefix     => $sex_check_outfile_prefix,
                read_freqfile_path => $read_freqfile_path,
                sex_check_min_f    => $sex_check_min_f,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $program_mode == 1 ) {

        while ( my ( $outfile_tag, $outfile_path ) = each %outfile_path ) {

            ## Collect QC metadata info for later use
            add_program_outfile_to_sample_info(
                {
                    path             => $outfile_path,
                    program_name     => $outfile_tag,
                    sample_info_href => $sample_info_href,
                }
            );

        }

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                path             => $stdout_file_path,
                program_name     => q{plink2},
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                dependency_method       => q{family_to_island},
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                job_id_chain            => $job_id_chain,
                recipe_file_path        => $recipe_file_path,
                sample_ids_ref          => \@{ $active_parameter_href->{sample_ids} },
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return;
}

1;
