package MIP::Recipes::Analysis::Gatk_gathervcfs;

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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_gatk_gathervcfs };

}

## Constants
Readonly my $ASTERISK    => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_gatk_gathervcfs {

## Function : Gather VCFs produced after gatk_genotypegvcfs done per contig.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}

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

    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Gnu::Coreutils qw(gnu_mv);
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Processmanagement::Processes qw{ print_wait };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_view_and_index_vcf };
    use MIP::Program::Variantcalling::Gatk
      qw{ gatk_gathervcfscloud gatk_selectvariants };
    use MIP::QC::Record
      qw{ add_processing_metafile_to_sample_info add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set program mode
    my $program_mode = $active_parameter_href->{$program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$program_name}{chain};
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};

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
    my $infamily_directory = my $outfamily_directory =
      catfile( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir, q{gatk}, );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$family_id}{$program_name}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$program_name}{file_tag};

    ## Files
    my $infile_prefix  = $family_id . $infile_tag . $call_type;
    my $outfile_prefix = $family_id . $outfile_tag . $call_type;

    ## Tags
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ### Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
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
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            call_type             => $call_type,
            core_number           => $core_number,
            directory_id          => $family_id,
            FILEHANDLE            => $FILEHANDLE,
            job_id_href           => $job_id_href,
            log                   => $log,
            process_time          => $time,
            program_directory     => catfile( $outaligner_dir, q{gatk} ),
            program_name          => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    my $process_batches_count = 1;

    ## Gather vcf files
    my @vcffile_paths;

  CONTIG:
    while ( my ( $contig_index, $contig ) =
        each @{ $file_info_href->{contigs} } )
    {

        $process_batches_count = print_wait(
            {
                FILEHANDLE            => $FILEHANDLE,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                process_counter       => $contig_index,
            }
        );

        ## Store infile
        push @vcffile_paths,
          $file_path_prefix . $UNDERSCORE . $contig . $infile_suffix;

        ## Copy file(s) to temporary directory
        migrate_file(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => catfile(
                    $infamily_directory,
                    $infile_prefix
                      . $UNDERSCORE
                      . $contig
                      . $infile_suffix
                      . $ASTERISK
                ),
                outfile_path => $temp_directory
            }
        );
        ## Store infile
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## GATK GatherVcfsCloud
    gatk_gathervcfscloud(
        {
            FILEHANDLE           => $FILEHANDLE,
            ignore_safety_checks => 0,
            infile_paths_ref     => \@vcffile_paths,
            memory_allocation    => q{Xmx4G},
            outfile_path         => $outfile_path_prefix . $outfile_suffix,
            temp_directory       => $temp_directory,
            verbosity => $active_parameter_href->{gatk_logging_level},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Produce a bcf compressed and index from vcf
    if ( $active_parameter_href->{gatk_gathervcfs_bcf_file} ) {

        # Exome analysis
        if ( $consensus_analysis_type eq q{wes} ) {

            say {$FILEHANDLE} q{### Remove extra reference samples};
            say {$FILEHANDLE} q{## GATK SelectVariants};
            gatk_selectvariants(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => $outfile_path_prefix . $outfile_suffix,
                    java_use_large_pages =>
                      $active_parameter_href->{java_use_large_pages},
                    memory_allocation => q{Xmx2g},
                    outfile_path      => $outfile_path_prefix
                      . $UNDERSCORE
                      . q{incnonvariantloci}
                      . $outfile_suffix,
                    referencefile_path => $referencefile_path,
                    sample_names_ref =>
                      \@{ $active_parameter_href->{sample_ids} },
                    temp_directory => $temp_directory,
                    verbosity => $active_parameter_href->{gatk_logging_level},
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            ## Move to original filename
            gnu_mv(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => $outfile_path_prefix
                      . $UNDERSCORE
                      . q{incnonvariantloci}
                      . $outfile_suffix,
                    outfile_path => $outfile_path_prefix . $outfile_suffix,
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }

        ## Reformat variant calling file and index
        bcftools_view_and_index_vcf(
            {
                FILEHANDLE          => $FILEHANDLE,
                infile_path         => $outfile_path_prefix . $outfile_suffix,
                outfile_path_prefix => $outfile_path_prefix,
                output_type         => q{b},
            }
        );

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $outfile_path_prefix . $DOT . q{bcf} . $ASTERISK,
                outfile_path => $outfamily_directory,
            }
        );
        say {$FILEHANDLE} q{wait}, $NEWLINE;
    }

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $outfile_path_prefix . $outfile_suffix . $ASTERISK,
            outfile_path => $outfamily_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        if ( $active_parameter_href->{gatk_gathervcfs_bcf_file} == 1 ) {

            my $program_gbcf_file_path =
              catfile( $outfamily_directory, $outfile_prefix . $DOT . q{bcf} );
            add_processing_metafile_to_sample_info(
                {
                    metafile_tag     => q{gbcf_file},
                    path             => $program_gbcf_file_path,
                    sample_info_href => $sample_info_href,
                }
            );
        }

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

1;
