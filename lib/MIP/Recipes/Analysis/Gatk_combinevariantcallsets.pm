package MIP::Recipes::Analysis::Gatk_combinevariantcallsets;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $ASTERISK $DOT $LOG_NAME $NEWLINE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.14;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_gatk_combinevariantcallsets };

}

sub analysis_gatk_combinevariantcallsets {

## Function : GATK CombineVariants to combine all variants call from different callers.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id => {
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
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Program::Gnu::Coreutils qw{ gnu_cp };
    use MIP::Language::Java qw{ java_core };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bcftools qw{ bcftools_view_and_index_vcf };
    use MIP::Program::Gatk qw{ gatk_combinevariants };
    use MIP::Sample_info qw{ set_file_path_to_store set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Stores callers that have been executed
    my @variant_callers;

    ## Only process active callers
  CALLER:
    foreach my $variant_caller ( @{ $parameter_href->{cache}{variant_callers} } ) {
        if ( $active_parameter_href->{$variant_caller} ) {

            push @variant_callers, $variant_caller;
        }
    }

    ## Stores the parallel chains that job ids should be inherited from
    my @parallel_chains;

    ## Unpack parameters
    my $gatk_jar =
      catfile( $active_parameter_href->{gatk_path}, q{GenomeAnalysisTK.jar} );
    my $job_id_chain = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $recipe_mode        = $active_parameter_href->{$recipe_name};
    my %recipe_resource    = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set and get the io files per chain, id and stream
    my %io = parse_io_outfiles(
        {
            chain_id               => $job_id_chain,
            id                     => $case_id,
            file_info_href         => $file_info_href,
            file_name_prefixes_ref => [$case_id],
            outdata_dir            => $active_parameter_href->{outdata_dir},
            parameter_href         => $parameter_href,
            recipe_name            => $recipe_name,
        }
    );

    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my $outfile_path        = $outfile_path_prefix . $outfile_suffix;

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ($recipe_file_path) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $recipe_resource{core_number},
            directory_id                    => $case_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            log                             => $log,
            memory_allocation               => $recipe_resource{memory},
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
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
                id             => $case_id,
                file_info_href => $file_info_href,
                parameter_href => $parameter_href,
                recipe_name    => $variant_caller,
                stream         => $stream,
            }
        );
        my $infile_path_prefix = $sample_io{$stream}{file_path_prefix};
        my $infile_suffix      = $sample_io{$stream}{file_suffix};
        my $infile_path        = $infile_path_prefix . $infile_suffix;

        ## Only use first part of name
        my ($variant_caller_prio_tag) = split /_/sxm, $variant_caller;

        ## Collect both tag and path in the same string
        $file_path{$variant_caller} = $variant_caller_prio_tag . $SPACE . $infile_path;
        ## For single caller use - collect infile path without prio tag
        $file_path{infile_path} = $infile_path;

        push @parallel_chains, $parameter_href->{$variant_caller}{chain};
    }

    my @combine_infile_paths = map { $file_path{$_} } @variant_callers;

    ## Check that we have something to combine
    if ( scalar @variant_callers > 1 ) {

        ## GATK CombineVariants
        say {$filehandle} q{## GATK CombineVariants};

        gatk_combinevariants(
            {
                exclude_nonvariants => 1,
                filehandle          => $filehandle,
                genotype_merge_option =>
                  $active_parameter_href->{gatk_combinevariants_genotype_merge_option},
                infile_paths_ref     => \@combine_infile_paths,
                java_jar             => $gatk_jar,
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                logging_level        => $active_parameter_href->{gatk_logging_level},
                memory_allocation    => q{Xmx2g},
                outfile_path         => $outfile_path,
                prioritize_caller =>
                  $active_parameter_href->{gatk_combinevariants_prioritize_caller},
                referencefile_path => $referencefile_path,
                temp_directory     => $temp_directory,
            }
        );
        say {$filehandle} $NEWLINE;
    }
    else {

        say {$filehandle} q{## Renaming case to facilitate downstream processing};

        gnu_cp(
            {
                filehandle   => $filehandle,
                infile_path  => $file_path{infile_path},
                outfile_path => $outfile_path,
            }
        );
        say {$filehandle} $NEWLINE;
    }

    if ( $active_parameter_href->{gatk_combinevariantcallsets_bcf_file} ) {

        ## Reformat variant calling file and index
        bcftools_view_and_index_vcf(
            {
                filehandle          => $filehandle,
                infile_path         => $outfile_path,
                outfile_path_prefix => $outfile_path_prefix,
                output_type         => q{b},
            }
        );
    }

    close $filehandle;

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        if ( $active_parameter_href->{gatk_combinevariantcallsets_bcf_file} ) {

            my $bcf_suffix    = $DOT . q{bcf};
            my $bcf_file_path = $outfile_path_prefix . $bcf_suffix;
            set_file_path_to_store(
                {
                    format           => q{bcf},
                    id               => $case_id,
                    path             => $bcf_file_path,
                    path_index       => $bcf_file_path . $DOT . q{csi},
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                }
            );
            set_file_path_to_store(
                {
                    format           => q{bcf},
                    id               => $case_id,
                    path             => $bcf_file_path,
                    path_index       => $bcf_file_path . $DOT . q{csi},
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                }
            );
        }

        submit_recipe(
            {
                base_command         => $profile_base_command,
                case_id              => $case_id,
                dependency_method    => q{sample_to_case},
                job_id_chain         => $job_id_chain,
                job_id_href          => $job_id_href,
                job_reservation_name => $active_parameter_href->{job_reservation_name},
                log                  => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                parallel_chains_ref => \@parallel_chains,
                recipe_file_path    => $recipe_file_path,
                sample_ids_ref      => \@{ $active_parameter_href->{sample_ids} },
                submission_profile  => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

1;
