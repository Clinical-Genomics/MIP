package MIP::Recipes::Gatk_genotypegvcfs;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use autodie qw{ :all };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Spec::Functions qw{ catdir catfile };

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_gatk_genotypegvcfs };

}

## Constants
Readonly my $NEWLINE     => qq{\n};
Readonly my $UNDERSCORE  => q{_};
Readonly my $DOT         => q{.};
Readonly my $ASTERISK    => q{*};
Readonly my $LONGER_TIME => 50;

sub analysis_gatk_genotypegvcfs {

## Function : GATK GenoTypeGVCFs.
## Returns  :
## Arguments: $parameter_href           => Parameter hash {REF}
##          : $active_parameter_href    => Active parameters for this analysis hash {REF}
##          : $sample_info_href         => Info on samples and family hash {REF}
##          : $file_info_href           => File info hash {REF}
##          : $infile_lane_prefix_href  => Infile(s) without the ".ending"
##          : $job_id_href              => Job id hash {REF}
##          : $family_id                => Family id
##          : $temp_directory           => Temporary directory
##          : $outfamily_directory      => Outfamily directory
##          : $outfamily_file_directory => Outfamily file directory
##          : $outaligner_dir           => Outaligner_dir used in the analysis {REF}
##          : $program_name             => Program name
##          : $call_type                => Call type

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id;
    my $temp_directory;
    my $outaligner_dir;
    my $call_type;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $outfamily_directory;
    my $outfamily_file_directory;
    my $program_name;

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
        outfamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfamily_directory,
        },
        outfamily_file_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfamily_file_directory,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
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
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Language::Java qw{ java_core };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_step_in_parallel_to_family };
    use MIP::Program::Variantcalling::Gatk qw{ gatk_genotypegvcfs };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $jobid_chain = $parameter_href->{$mip_program_name}{chain};
    say {$FILEHANDLE} "##Jobid_chain:$jobid_chain##";

    my $time        = $active_parameter_href->{module_time}{$mip_program_name};
    if ( $active_parameter_href->{gatk_genotypegvcfs_all_sites} == 1 ) {

        #Including all sites requires longer processing time
        $time = $LONGER_TIME;
    }
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};

#gatk genotype is most safely processed in single thread mode, , but we need some java heap allocation
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};

    my $sbatch_script_tracker = 0;

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    #Used downstream
    $parameter_href->{$mip_program_name}{$family_id}{indirectory} =
      $outfamily_directory;

    ## Files
    my $outfile_tag =
      $file_info_href->{$family_id}{$mip_program_name}{file_tag};

    ## Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{outfile_suffix},
            program_name   => q{pgatk_haplotypecaller},
        }
    );

    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
            jobid_chain    => $jobid_chain,
            file_suffix => $parameter_href->{$mip_program_name}{outfile_suffix},
        }
    );

    ## Create .fam file to be used in variant calling analyses
    my $fam_file_path =
      catfile( $outfamily_file_directory, $family_id . $DOT . q{fam} );
    create_fam_file(
        {
            parameter_href        => $parameter_href,
            active_parameter_href => $active_parameter_href,
            sample_info_href      => $sample_info_href,
            FILEHANDLE            => $FILEHANDLE,
            fam_file_path         => $fam_file_path,

        }
    );

    ## Split per contig
    foreach my $contig ( @{ $file_info_href->{contigs} } ) {

        ## Assign file_tags
        my $outfile_prefix =
          $family_id . $outfile_tag . $call_type . $UNDERSCORE . $contig;
        my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        my ($file_name) = program_prerequisites(
            {
                active_parameter_href => $active_parameter_href,
                job_id_href           => $job_id_href,
                FILEHANDLE            => $FILEHANDLE,
                directory_id          => $family_id,
                program_name          => $program_name,
                program_directory => catfile( lc($outaligner_dir), q{gatk} ),
                call_type         => $call_type,
                core_number       => $core_number,
                process_time      => $time,
                temp_directory    => $temp_directory,
                sleep             => 1
                , #Let process sleep randomly for 0-60 seconds to avoid race condition
            }
        );

        ## Collect infiles for all sample_ids to enable migration to temporary directory
        my @file_paths;
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Add merged infile name after merging all BAM files per sample_id
            my $infile = $file_info_href->{$sample_id}{merge_infile};    #Alias

            ## Assign directories
            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $outaligner_dir, q{gatk} );

            ## Assign file_tags
            my $infile_tag =
              $file_info_href->{$sample_id}{pgatk_haplotypecaller}{file_tag};
            my $infile_prefix = $infile . $infile_tag . $UNDERSCORE . $contig;

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
                active_parameter_href => $active_parameter_href,
                fam_file_path         => catfile(
                    $outfamily_file_directory, $family_id . $DOT . q{fam}
                ),
                program_name => $program_name,
            }
        );
        ## Files for genotypegvcfs
        if (   ( $consensus_analysis_type eq q{wes} )
            || ( $consensus_analysis_type eq q{rapid} ) )
        {

            push @file_paths,
              $active_parameter_href->{gatk_genotypegvcfs_ref_gvcf};
        }

        ## Writes java core commands to filehandle.
        java_core(
            {
                FILEHANDLE        => $FILEHANDLE,
                memory_allocation => q{Xmx24g},
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                temp_directory => $temp_directory,
                java_jar       => catfile(
                    $active_parameter_href->{gatk_path},
                    q{GenomeAnalysisTK.jar}
                ),
            }
        );

        gatk_genotypegvcfs(
            {
                intervals_ref    => [$contig],
                infile_paths_ref => \@file_paths,
                outfile_path     => $outfile_path_prefix . $outfile_suffix,
                logging_level => $active_parameter_href->{gatk_logging_level},
                referencefile_path =>
                  $active_parameter_href->{human_genome_reference},
                dbsnp =>
                  $active_parameter_href->{gatk_haplotypecaller_snp_known_set},
                pedigree_validation_type => $commands{pedigree_validation_type},
                pedigree                 => $commands{pedigree},
                include_nonvariant_sites =>
                  $active_parameter_href->{gatk_genotypegvcfs_all_sites},
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                infile_path => $outfile_path_prefix
                  . $outfile_suffix
                  . $ASTERISK,
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} q{wait} . $NEWLINE;

        close $FILEHANDLE;

        if (   ( $mip_program_mode == 1 )
            && ( !$active_parameter_href->{dry_run_all} ) )
        {

            submit_job(
                {
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    job_id_href             => $job_id_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    dependencies => q{sample_id_dependency_step_in_parallel},
                    path         => $jobid_chain,
                    sbatch_file_name      => $file_name,
                    sbatch_script_tracker => $sbatch_script_tracker
                }
            );
        }
        $sbatch_script_tracker++;    #Tracks nr of sbatch scripts

    }

    return;
}

1;
