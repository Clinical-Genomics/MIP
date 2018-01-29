package MIP::Recipes::Analysis::Gatk_genotypegvcfs;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
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
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_gatk_genotypegvcfs };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_gatk_genotypegvcfs {

## Function : GATK GenoTypeGVCFs.
## Returns  :
## Arguments: $active_parameter_href    => Active parameters for this analysis hash {REF}
##          : $call_type                => Call type
##          : $family_id                => Family id
##          : $file_info_href           => File info hash {REF}
##          : $infile_lane_prefix_href  => Infile(s) without the ".ending"
##          : $job_id_href              => Job id hash {REF}
##          : $outaligner_dir           => Outaligner_dir used in the analysis
##          : $outfamily_directory      => Outfamily directory
##          : $outfamily_file_directory => Outfamily file directory
##          : $parameter_href           => Parameter hash {REF}
##          : $program_name             => Program name
##          : $sample_info_href         => Info on samples and family hash {REF}
##          : $temp_directory           => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $outfamily_directory;
    my $outfamily_file_directory;
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
        outfamily_directory => {
            defined     => 1,
            required    => 1,
            store       => \$outfamily_directory,
            strict_type => 1,
        },
        outfamily_file_directory => {
            defined     => 1,
            required    => 1,
            store       => \$outfamily_file_directory,
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

    use MIP::File::Format::Pedigree qw{ create_fam_file gatk_pedigree_flag };
    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Language::Java qw{ java_core };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_step_in_parallel_to_family };
    use MIP::Program::Variantcalling::Gatk qw{ gatk_genotypegvcfs };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Set::File qw{ set_file_suffix };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};

    ## Gatk genotype is most safely processed in single thread mode, , but we need some java heap allocation
    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    ## Constants
    Readonly my $INCLUDE_NONVARIANT_SITES_TIME => 50;

    # If all sites should be included
    if ( $active_parameter_href->{gatk_genotypegvcfs_all_sites} == 1 ) {

        # Including all sites requires longer processing time
        $time = $INCLUDE_NONVARIANT_SITES_TIME;
    }

    my $sbatch_script_tracker = 0;

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    # Used downstream
    $parameter_href->{$mip_program_name}{$family_id}{indirectory} =
      $outfamily_directory;

    ## Tags
    my $outfile_tag =
      $file_info_href->{$family_id}{$mip_program_name}{file_tag};

    ## Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            program_name   => q{pgatk_haplotypecaller},
            suffix_key     => q{outfile_suffix},
        }
    );

    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            file_suffix => $parameter_href->{$mip_program_name}{outfile_suffix},
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

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

  CONTIG:
    foreach my $contig ( @{ $file_info_href->{contigs} } ) {

        ## Assign file_tags
        my $outfile_prefix =
          $family_id . $outfile_tag . $call_type . $UNDERSCORE . $contig;
        my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        my ($file_path) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                call_type             => $call_type,
                core_number           => $core_number,
                directory_id          => $family_id,
                FILEHANDLE            => $FILEHANDLE,
                job_id_href           => $job_id_href,
                process_time          => $time,
                program_directory     => catfile( $outaligner_dir, q{gatk} ),
                program_name          => $program_name,
                sleep                 => 1,
                source_environment_commands_ref => [$source_environment_cmd],
                temp_directory                  => $temp_directory,
            }
        );

        ## Collect infiles for all sample_ids to enable migration to temporary directory
        my @file_paths;

        ## Add merged infile name after merging all BAM files per sample_id
      SAMPLE:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Add merged infile name prefix after merging all BAM files per sample_id
            my $merged_infile_prefix = get_merged_infile_prefix(
                {
                    file_info_href => $file_info_href,
                    sample_id      => $sample_id,
                }
            );

            ## Assign directories
            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $outaligner_dir, q{gatk} );

            ## Assign file_tags
            my $infile_tag =
              $file_info_href->{$sample_id}{pgatk_haplotypecaller}{file_tag};
            my $infile_prefix =
              $merged_infile_prefix . $infile_tag . $UNDERSCORE . $contig;

            ## Collect for downstream use
            push
              @file_paths,
              catfile( $temp_directory, $infile_prefix . $infile_suffix );

            ## Copy file(s) to temporary directory
            say {$FILEHANDLE} q{## Copy file(s) to temporary directory};

            $insample_directory = catfile( $insample_directory,
                $infile_prefix . $infile_suffix . $ASTERISK );

            migrate_file(
                {
                    FILEHANDLE   => $FILEHANDLE,
                    infile_path  => $insample_directory,
                    outfile_path => $temp_directory,
                }
            );
            say {$FILEHANDLE} q{wait} . $NEWLINE;

        }

        ## GATK GenoTypeGVCFs
        say {$FILEHANDLE} q{## GATK GenoTypeGVCFs};

        ## Get parameters
        ## Check if "--pedigree" and "--pedigreeValidationType" should be included in analysis
        my %commands = gatk_pedigree_flag(
            {
                fam_file_path => catfile(
                    $outfamily_file_directory, $family_id . $DOT . q{fam}
                ),
                program_name => $program_name,
            }
        );
        ## Files for genotypegvcfs
        if ( $consensus_analysis_type eq q{wes} ) {

            push @file_paths,
              $active_parameter_href->{gatk_genotypegvcfs_ref_gvcf};
        }

        gatk_genotypegvcfs(
            {
                dbsnp_file_path =>
                  $active_parameter_href->{gatk_haplotypecaller_snp_known_set},
                FILEHANDLE       => $FILEHANDLE,
                intervals_ref    => [$contig],
                infile_paths_ref => \@file_paths,
                include_nonvariant_sites =>
                  $active_parameter_href->{gatk_genotypegvcfs_all_sites},
                java_jar => catfile(
                    $active_parameter_href->{gatk_path},
                    q{GenomeAnalysisTK.jar}
                ),
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                logging_level => $active_parameter_href->{gatk_logging_level},
                memory_allocation => q{Xmx24g},
                outfile_path      => $outfile_path_prefix . $outfile_suffix,
                pedigree_validation_type => $commands{pedigree_validation_type},
                pedigree                 => $commands{pedigree},
                referencefile_path =>
                  $active_parameter_href->{human_genome_reference},
                temp_directory => $temp_directory,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $outfile_path_prefix
                  . $outfile_suffix
                  . $ASTERISK,
                outfile_path => $outfamily_directory,
            }
        );
        say {$FILEHANDLE} q{wait} . $NEWLINE;

        close $FILEHANDLE;

        if ( $mip_program_mode == 1 ) {

            slurm_submit_job_sample_id_dependency_step_in_parallel_to_family(
                {
                    family_id               => $family_id,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    log                     => $log,
                    path                    => $job_id_chain,
                    sample_ids_ref =>
                      \@{ $active_parameter_href->{sample_ids} },
                    sbatch_file_name      => $file_path,
                    sbatch_script_tracker => $sbatch_script_tracker,
                }
            );
        }
        $sbatch_script_tracker++;    # Tracks nr of sbatch scripts
    }
    return;
}

1;
