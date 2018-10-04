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
    our $VERSION = 1.05;

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
##          : $family_id                => Family id
##          : $file_info_href           => File info hash {REF}
##          : $infile_lane_prefix_href  => Infile(s) without the ".ending"
##          : $job_id_href              => Job id hash {REF}
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
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_step_in_parallel_to_family };
    use MIP::Program::Variantcalling::Gatk
      qw{ gatk_genomicsdbimport  gatk_genotypegvcfs };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $job_id_chain = get_program_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            program_name   => $program_name,
        }
    );
    my $program_mode = $active_parameter_href->{$program_name};
    my $sbatch_script_tracker = 0;

    ## Gatk genotype is most safely processed in single thread mode, , but we need some java heap allocation
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

        ## Set and get the io files per chain, id and stream
    my %io = parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $family_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $family_id,
                iterators_ref    => $file_info_href->{contigs_size_ordered},
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                program_name     => $program_name,
                temp_directory   => $temp_directory,
            }
    );
    my @outfile_paths       = @{ $io{out}{file_paths} };
    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_name_prefix = $io{out}{file_name_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my %outfile_path   = %{ $io{out}{file_path_href} };
    my $temp_outfile_path_prefix = $io{temp}{file_path_prefix};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Create .fam file to be used in variant calling analyses
    my $fam_file_path =
      catfile( $outdir_path_prefix, $family_id . $DOT . q{fam} );
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

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        my ($file_path) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                core_number           => $core_number,
                directory_id          => $family_id,
                FILEHANDLE            => $FILEHANDLE,
                job_id_href           => $job_id_href,
                log                   => $log,
                process_time          => $time,
                program_directory     => $program_name,
                program_name          => $program_name,
                sleep                 => 1,
                source_environment_commands_ref => \@source_environment_cmds,
                temp_directory                  => $temp_directory,
            }
        );

        ## Collect infiles for all sample_ids to enable migration to temporary directory
	my @genotype_temp_infile_paths;
	my $process_batches_count = 1;
    while ( my ( $sample_id_index, $sample_id ) =
        each @{ $active_parameter_href->{sample_ids} } )
    {

        $process_batches_count = print_wait(
            {
                FILEHANDLE            => $FILEHANDLE,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                process_counter       => $sample_id_index,
            }
        );

	       ## Get the io infiles per chain and id
        my %sample_io = get_io_files(
            {
                id             => $sample_id,
                file_info_href => $file_info_href,
                parameter_href => $parameter_href,
                program_name   => $program_name,
                stream         => q{in},
                temp_directory => $temp_directory,
            }
        );
        my $infile_path_prefix = $sample_io{in}{file_path_prefix};
        my $infile_suffix      = $sample_io{in}{file_suffix};
        my $infile_path =
          $infile_path_prefix . substr( $infile_suffix, 0, 2 ) . $ASTERISK;
        my $temp_infile_path_prefix = $sample_io{temp}{file_path_prefix};
        my $temp_infile_path        = $temp_infile_path_prefix . $infile_suffix;

        ## Store temp infile path for each sample_id
        push @genotype_temp_infile_paths, $temp_infile_path;

            ## Copy file(s) to temporary directory
            say {$FILEHANDLE} q{## Copy file(s) to temporary directory};

            migrate_file(
                {
                    FILEHANDLE   => $FILEHANDLE,
                    infile_path  => $infile_path,
                    outfile_path => $temp_directory,
                }
            );
            say {$FILEHANDLE} q{wait} . $NEWLINE;

      }

        ## GATK GenomicsDBImport
        say {$FILEHANDLE} q{## GATK GenomicsDBImport};

        ## Files to import into GenomicsDB
        if ( $consensus_analysis_type eq q{wes} ) {

            push @genotype_temp_infile_paths,
              $active_parameter_href->{gatk_genotypegvcfs_ref_gvcf};
        }

        gatk_genomicsdbimport(
            {
                FILEHANDLE       => $FILEHANDLE,
                intervals_ref    => [$contig],
                infile_paths_ref => \@genotype_temp_infile_paths,
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                verbosity => $active_parameter_href->{gatk_logging_level},
                genomicsdb_workspace_path => $temp_outfile_path_prefix
                  . $UNDERSCORE . q{DB},
                referencefile_path =>
                  $active_parameter_href->{human_genome_reference},
                temp_directory    => $temp_directory,
                memory_allocation => q{Xmx8g},
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## GATK GenoTypeGVCFs
        say {$FILEHANDLE} q{## GATK GenoTypeGVCFs};

        gatk_genotypegvcfs(
            {
                dbsnp_path =>
                  $active_parameter_href->{gatk_haplotypecaller_snp_known_set},
                FILEHANDLE    => $FILEHANDLE,
                intervals_ref => [$contig],
                infile_path   => q{gendb://}
                  . $temp_outfile_path_prefix
                  . $UNDERSCORE . q{DB},
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                memory_allocation => q{Xmx8g},
                outfile_path      => $outfile_path{$contig},
                pedigree          => $fam_file_path,
                referencefile_path =>
                  $active_parameter_href->{human_genome_reference},
                temp_directory => $temp_directory,
                verbosity      => $active_parameter_href->{gatk_logging_level},
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        close $FILEHANDLE;

        if ( $program_mode == 1 ) {

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
