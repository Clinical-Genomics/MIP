package MIP::Recipes::Analysis::Gatk_haplotypecaller;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use POSIX qw{ floor };
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
    our $VERSION = 1.15;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_gatk_haplotypecaller };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_gatk_haplotypecaller {

## Function : Gatk haplotypecaller analysis recipe
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
    my $temp_directory;
    my $xargs_file_counter;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id_ref => {
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
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
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
        xargs_file_counter => {
            allow       => qr{ \A\d+\z }xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw{ get_memory_constrained_core_number };
    use MIP::File::Format::Pedigree qw{ create_fam_file };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter
      qw{ get_gatk_intervals get_recipe_parameters get_recipe_attributes };
    use MIP::IO::Files qw{ xargs_migrate_contig_files };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Alignment::Gatk qw{ gatk_haplotypecaller };
    use MIP::Program::Variantcalling::Gatk qw{ gatk_gathervcfscloud };
    use MIP::QC::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Set::File qw{ set_io_files };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $JAVA_MEMORY_ALLOCATION        => 8;
    Readonly my $STANDARD_MIN_CONFIDENCE_THRSD => 10;

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my $indir_path_prefix       = $io{in}{dir_path_prefix};
    my $infile_suffix           = $io{in}{file_suffix};
    my $infile_name_prefix      = $io{in}{file_name_prefix};
    my $temp_infile_name_prefix = $io{temp}{file_name_prefix};
    my %temp_infile_path        = %{ $io{temp}{file_path_href} };

    my $job_id_chain = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $recipe_mode        = $active_parameter_href->{$recipe_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $analysis_type      = $active_parameter_href->{analysis_type}{$sample_id};
    my $xargs_file_path_prefix;

    ## Get module parameters
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    # Constrain parallelization to match available memory
    my $program_core_number = get_memory_constrained_core_number(
        {
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
            memory_allocation  => $JAVA_MEMORY_ALLOCATION,
            node_ram_memory    => $active_parameter_href->{node_ram_memory},
            recipe_core_number => $core_number,
        }
    );

    ## Outpaths
    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $sample_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => $file_info_href->{contigs_size_ordered},
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
                temp_directory   => $temp_directory,
            }
        )
    );
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my %outfile_path        = %{ $io{out}{file_path_href} };
    my %temp_outfile_path   = %{ $io{temp}{file_path_href} };

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $sample_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    # For ".fam" file
    my $outcase_file_directory =
      catdir( $active_parameter_href->{outdata_dir}, $case_id );

    ## Create .fam file to be used in variant calling analyses
    my $fam_file_path = catfile( $outcase_file_directory, $case_id . $DOT . q{fam} );
    create_fam_file(
        {
            active_parameter_href => $active_parameter_href,
            fam_file_path         => $fam_file_path,
            FILEHANDLE            => $FILEHANDLE,
            parameter_href        => $parameter_href,
            log                   => $log,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Generate gatk intervals. Chromosomes for WGS/WTS and paths to contig_bed_files for WES
    my %gatk_intervals = get_gatk_intervals(
        {
            analysis_type         => $analysis_type,
            contigs_ref           => \@{ $file_info_href->{contigs_size_ordered} },
            FILEHANDLE            => $FILEHANDLE,
            max_cores_per_node    => $core_number,
            outdirectory          => $temp_directory,
            reference_dir         => $active_parameter_href->{reference_dir},
            exome_target_bed_href => $active_parameter_href->{exome_target_bed},
            file_ending           => $file_info_href->{exome_target_bed}[1],
            log                   => $log,
            sample_id             => $sample_id,
        }
    );

    ## Set the PCR indel model for haplotypecaller
    my $pcr_indel_model = _get_pcr_indel_model(
        {
            analysis_type         => $analysis_type,
            active_parameter_href => $active_parameter_href,
        }
    );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    ($xargs_file_counter) = xargs_migrate_contig_files(
        {
            contigs_ref        => \@{ $file_info_href->{contigs_size_ordered} },
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_ending        => substr( $infile_suffix, 0, 2 ) . $ASTERISK,
            file_path          => $recipe_file_path,
            indirectory        => $indir_path_prefix,
            infile             => $infile_name_prefix,
            recipe_info_path   => $recipe_info_path,
            temp_directory     => $temp_directory,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## GATK HaplotypeCaller
    say {$FILEHANDLE} q{## GATK HaplotypeCaller};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $program_core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## GATK Haplotypecaller
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        gatk_haplotypecaller(
            {
                annotations_ref =>
                  \@{ $active_parameter_href->{gatk_haplotypecaller_annotation} },
                dbsnp_path =>
                  $active_parameter_href->{gatk_haplotypecaller_snp_known_set},
                dont_use_soft_clipped_bases =>
                  $active_parameter_href->{gatk_haplotypecaller_no_soft_clipped_bases},
                emit_ref_confidence =>
                  $active_parameter_href->{gatk_haplotypecaller_emit_ref_confidence},
                FILEHANDLE           => $XARGSFILEHANDLE,
                infile_path          => $temp_infile_path{$contig},
                intervals_ref        => $gatk_intervals{$contig},
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                memory_allocation    => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
                outfile_path         => $temp_outfile_path{$contig},
                pcr_indel_model      => $pcr_indel_model,
                pedigree             => $fam_file_path,
                referencefile_path   => $referencefile_path,
                standard_min_confidence_threshold_for_calling =>
                  $STANDARD_MIN_CONFIDENCE_THRSD,
                stderrfile_path => $stderrfile_path,
                temp_directory  => $temp_directory,
                verbosity       => $active_parameter_href->{gatk_logging_level},
                use_new_qual_calculator =>
                  $active_parameter_href->{gatk_use_new_qual_calculator},
                xargs_mode => 1,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## Concatenate contig VCFs
    my $concat_vcf_path = $outfile_path_prefix . $outfile_suffix;

    ## GATK GatherVcfsCloud
    my @contig_vcf_paths;
  CONTIG:
    foreach my $contig ( @{ $file_info_href->{contigs} } ) {

        push @contig_vcf_paths, $temp_outfile_path{$contig};
    }

    gatk_gathervcfscloud(
        {
            FILEHANDLE           => $FILEHANDLE,
            ignore_safety_checks => 0,
            infile_paths_ref     => \@contig_vcf_paths,
            memory_allocation    => q{Xmx4G},
            outfile_path         => $concat_vcf_path,
            temp_directory       => $temp_directory,
            verbosity            => $active_parameter_href->{gatk_logging_level},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    close $FILEHANDLE;
    close $XARGSFILEHANDLE;

    ## Set input files for next module
    set_io_files(
        {
            chain_id       => $job_id_chain,
            id             => $sample_id,
            file_info_href => $file_info_href,
            file_paths_ref => [$concat_vcf_path],
            recipe_name    => $recipe_name,
            stream         => q{out},
        }
    );

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => basename($concat_vcf_path),
                path             => $concat_vcf_path,
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command            => $profile_base_command,
                case_id                 => $case_id,
                dependency_method       => q{sample_to_sample},
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                job_id_chain            => $job_id_chain,
                recipe_file_path        => $recipe_file_path,
                sample_id               => $sample_id,
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub _get_pcr_indel_model {

## Function : Set the pcr indel model for GATK HaplotypeCaller
## Returns  : $pcr_indel_model
## Argumetns: active_parameter_href => Active parameter hash {REF}
##          : $analysis_type        => Analysis type

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $analysis_type;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$analysis_type,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $pcr_indel_model;

    ## Leave $pcr_indel_model as undef for WES
    if ( not $analysis_type eq q{wes} ) {
        $pcr_indel_model = $active_parameter_href->{gatk_haplotypecaller_pcr_indel_model};
    }

    return $pcr_indel_model;
}

1;
