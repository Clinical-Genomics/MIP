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
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_gatk_combinevariantcallsets };

}

## Constants
Readonly my $ASTERIX    => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_gatk_combinevariantcallsets {

## Function : GATK CombineVariants to combine all variants call from different callers.
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

    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Language::Java qw{ java_core };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_view_and_index_vcf };
    use MIP::Program::Variantcalling::Gatk qw{ gatk_combinevariants };
    use MIP::QC::Record
      qw{ add_program_outfile_to_sample_info add_processing_metafile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };

    ## Stores callers that have been executed
    my @variant_callers;

    ## Stores the parallel chains that jobIds should be inherited from
    my @parallel_chains;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set program mode
    my $program_mode = $active_parameter_href->{$program_name};

    ## Alias
    my $job_id_chain       = $parameter_href->{$program_name}{chain};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $gatk_jar =
      catfile( $active_parameter_href->{gatk_path}, q{GenomeAnalysisTK.jar} );
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

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            call_type                       => $call_type,
            core_number                     => $core_number,
            directory_id                    => $family_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            program_directory               => catfile($outaligner_dir),
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ## Assign directories
    my $outfamily_directory = catfile( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir );

    # Used downstream
    $parameter_href->{$program_name}{indirectory} = $outfamily_directory;

    ## Assign file_tags
    my $outfile_tag =
      $file_info_href->{$family_id}{gatk_combinevariantcallsets}{file_tag};

    ## Will be set downstream
    my @infile_tags_and_paths;

    ## Files
    my $outfile_prefix = $family_id . $outfile_tag . $call_type;

    ## Paths
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ### Assign suffix
    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
            job_id_chain   => $job_id_chain,
            file_suffix    => $parameter_href->{$program_name}{outfile_suffix},
        }
    );

    ## Collect file info and migrate files
  VARIANT_CALLER:
    foreach my $variant_caller (
        @{ $parameter_href->{dynamic_parameter}{variant_callers} } )
    {

        next VARIANT_CALLER
          if ( not $active_parameter_href->{$variant_caller} );

        ### Expect vcf
        ## Assign directories
        my $program_outdirectory_name =
          $parameter_href->{$variant_caller}{outdir_name};
        my $infamily_directory = catfile( $active_parameter_href->{outdata_dir},
            $family_id, $outaligner_dir, $program_outdirectory_name );

        ## Assign file_tags
        my $infile_tag =
          $file_info_href->{$family_id}{$variant_caller}{file_tag};
        my $infile_prefix = $family_id . $infile_tag . $call_type;

        ## Assign suffix
        my $infile_suffix = get_file_suffix(
            {
                parameter_href => $parameter_href,
                program_name   => $variant_caller,
                suffix_key     => q{outfile_suffix},
            }
        );

        ## Collect both tag and path in the same string
        my $path = catfile( $temp_directory, $infile_prefix . $infile_suffix );
        push @infile_tags_and_paths,
          $program_outdirectory_name . $SPACE . $path;

        ## To prioritize downstream - 1. gatk 2. samtools etc.
        ## Determined by order_parameters order
        unshift @variant_callers, $program_outdirectory_name;

        ## Do not add MAIN chains
        if ( not $parameter_href->{$variant_caller}{chain} eq q{MAIN} ) {

            push @parallel_chains, $parameter_href->{$variant_caller}{chain};
        }

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        migrate_file(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => catfile(
                    $infamily_directory,
                    $infile_prefix . $infile_suffix . $ASTERIX
                ),
                outfile_path => $temp_directory
            }
        );
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## GATK CombineVariants
    say {$FILEHANDLE} q{## GATK CombineVariants};

    gatk_combinevariants(
        {
            exclude_nonvariants   => 1,
            FILEHANDLE            => $FILEHANDLE,
            genotype_merge_option => $active_parameter_href
              ->{gatk_combinevariants_genotype_merge_option},
            infile_paths_ref => \@infile_tags_and_paths,
            java_jar         => $gatk_jar,
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            logging_level     => $active_parameter_href->{gatk_logging_level},
            memory_allocation => q{Xmx2g},
            outfile_path      => $outfile_path_prefix . $outfile_suffix,
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
                infile_path  => $outfile_path_prefix . $DOT . q{bcf} . $ASTERIX,
                outfile_path => $outfamily_directory,
            }
        );
    }

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $outfile_path_prefix . $outfile_suffix,
            outfile_path => $outfamily_directory,
        }
    );
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

        my $most_complete_format_key =
          q{most_complete} . $UNDERSCORE . substr $outfile_suffix, 1;
        add_processing_metafile_to_sample_info(
            {
                metafile_tag     => $most_complete_format_key,
                path             => $program_outfile_path,
                sample_info_href => $sample_info_href,
            }
        );

        if ( $active_parameter_href->{gatk_combinevariantcallsets_bcf_file} ) {

            my $bcf_suffix = $DOT . q{bcf};
            my $most_complete_bcf_key =
              q{most_complete} . $UNDERSCORE . substr $bcf_suffix, 1;
            my $program_bcf_file_path =
              catfile( $outfamily_directory, $outfile_prefix . $bcf_suffix );
            add_processing_metafile_to_sample_info(
                {
                    metafile_tag     => $most_complete_bcf_key,
                    path             => $program_bcf_file_path,
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
