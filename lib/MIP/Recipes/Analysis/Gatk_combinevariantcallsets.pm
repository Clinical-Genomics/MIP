package MIP::Recipes::Analysis::Gatk_combinevariantcallsets;

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
    our @EXPORT_OK = qw{ analysis_gatk_combinevariantcallsets };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_gatk_combinevariantcallsets {

## Function : GATK CombineVariants to combine all variants call from different callers.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File info hash {REF}
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

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter
      qw{ get_module_parameters get_program_parameters get_program_attributes };
    use MIP::Language::Java qw{ java_core };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_view_and_index_vcf };
    use MIP::Program::Variantcalling::Gatk qw{ gatk_combinevariants };
    use MIP::QC::Record
      qw{ add_program_outfile_to_sample_info add_processing_metafile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Stores callers that have been executed
    my @variant_callers;

    ## Only process active callers
    foreach my $variant_caller (
        @{ $parameter_href->{dynamic_parameter}{variant_callers} } )
    {
        if ( $active_parameter_href->{$variant_caller} ) {

            push @variant_callers, $variant_caller;
        }
    }

    ## Stores the parallel chains that job ids should be inherited from
    my @parallel_chains;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    my $gatk_jar =
      catfile( $active_parameter_href->{gatk_path}, q{GenomeAnalysisTK.jar} );
    my $job_id_chain = get_program_attributes(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
            attribute      => q{chain},
        }
    );
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $program_mode       = $active_parameter_href->{$program_name};
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
            chain_id               => $job_id_chain,
            id                     => $family_id,
            file_info_href         => $file_info_href,
            file_name_prefixes_ref => [$family_id],
            outdata_dir            => $active_parameter_href->{outdata_dir},
            parameter_href         => $parameter_href,
            program_name           => $program_name,
            temp_directory         => $temp_directory,
        }
    );

    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my $outfile_path        = $outfile_path_prefix . $outfile_suffix;

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
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
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    ## Collect infiles for all sample_ids for joint calling program
    # Paths for variant callers to be merged
    my %file_path;
    my $stream = q{out};

    ## Collect file info and migrate files
  VARIANT_CALLER:
    foreach my $variant_caller (@variant_callers) {

        ## Get the io infiles per chain and id
        my %sample_io = get_io_files(
            {
                id             => $family_id,
                file_info_href => $file_info_href,
                parameter_href => $parameter_href,
                program_name   => $variant_caller,
                stream         => $stream,
                temp_directory => $temp_directory,
            }
        );
        my $infile_path_prefix = $sample_io{$stream}{file_path_prefix};
        my $infile_suffix      = $sample_io{$stream}{file_suffix};
        my $infile_path        = $infile_path_prefix . $infile_suffix;

        ## Only use first part of name
        my ($variant_caller_prio_tag) = split /_/sxm, $variant_caller;
        ## Collect both tag and path in the same string
        $file_path{$variant_caller} =
          $variant_caller_prio_tag . $SPACE . $infile_path;

        push @parallel_chains, $parameter_href->{$variant_caller}{chain};
    }

    ## GATK CombineVariants
    say {$FILEHANDLE} q{## GATK CombineVariants};

    my @combine_infile_paths = map { $file_path{$_} } @variant_callers;
    gatk_combinevariants(
        {
            exclude_nonvariants   => 1,
            FILEHANDLE            => $FILEHANDLE,
            genotype_merge_option => $active_parameter_href
              ->{gatk_combinevariants_genotype_merge_option},
            infile_paths_ref => \@combine_infile_paths,
            java_jar         => $gatk_jar,
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            logging_level     => $active_parameter_href->{gatk_logging_level},
            memory_allocation => q{Xmx2g},
            outfile_path      => $outfile_path,
            prioritize_caller =>
              $active_parameter_href->{gatk_combinevariants_prioritize_caller},
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    if ( $active_parameter_href->{gatk_combinevariantcallsets_bcf_file} ) {

        ## Reformat variant calling file and index
        bcftools_view_and_index_vcf(
            {
                FILEHANDLE          => $FILEHANDLE,
                infile_path         => $outfile_path,
                outfile_path_prefix => $outfile_path_prefix,
                output_type         => q{b},
            }
        );
    }

    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                path             => $outfile_path,
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );

        my $most_complete_format_key =
          q{most_complete} . $UNDERSCORE . substr $outfile_suffix, 1;
        add_processing_metafile_to_sample_info(
            {
                metafile_tag     => $most_complete_format_key,
                path             => $outfile_path,
                sample_info_href => $sample_info_href,
            }
        );

        if ( $active_parameter_href->{gatk_combinevariantcallsets_bcf_file} ) {

            my $bcf_suffix = $DOT . q{bcf};
            my $most_complete_bcf_key =
              q{most_complete} . $UNDERSCORE . substr $bcf_suffix, 1;
            my $bcf_file_path = $outfile_path_prefix . $bcf_suffix;
            add_processing_metafile_to_sample_info(
                {
                    metafile_tag     => $most_complete_bcf_key,
                    path             => $bcf_file_path,
                    sample_info_href => $sample_info_href,
                }
            );
        }

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                parallel_chains_ref     => \@parallel_chains,
                path                    => $job_id_chain,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                sbatch_file_name => $file_path,
            }
        );
    }
    return;
}

1;
