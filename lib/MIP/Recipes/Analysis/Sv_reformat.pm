package MIP::Recipes::Analysis::Sv_reformat;

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

# MIPs lib/
use MIP::Constants qw{ $DOT $EMPTY_STR $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.14;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_reformat_sv };

}

## Constants
Readonly my $JAVA_MEMORY_ALLOCATION => 20;

sub analysis_reformat_sv {

## Function : Concatenate and sort contig files. Optionally remove variants from genelist.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $filehandle              => Sbatch filehandle to write to
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $reference_dir           => MIP reference directory
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $sample_info_href;
    my $parameter_href;
    my $recipe_name;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
    my $reference_dir;
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
        reference_dir => {
            default     => $arg_href->{active_parameter_href}{reference_dir},
            store       => \$reference_dir,
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

    use MIP::Get::Analysis qw{ get_vcf_parser_analysis_suffix };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Gnu::Software::Gnu_grep qw{ gnu_grep };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bcftools qw{ bcftools_view_and_index_vcf };
    use MIP::Program::Picardtools qw{ picardtools_sortvcf };
    use MIP::Sample_info qw{ set_file_path_to_store
      set_most_complete_vcf
      set_recipe_metafile_in_sample_info
      set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my %io = get_io_files(
        {
            id             => $case_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
        }
    );

    my $infile_name_prefix = $io{in}{file_name_prefix};
    my @infile_paths       = @{ $io{in}{file_paths} };

    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my $job_id_chain            = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    my @vcfparser_analysis_types = get_vcf_parser_analysis_suffix(
        {
            vcfparser_outfile_count =>
              $active_parameter_href->{sv_vcfparser_outfile_count},
        }
    );

    ## Set and get the io files per chain, id and stream
    my @set_outfile_name_prefixes =
      map { $infile_name_prefix . $_ } @vcfparser_analysis_types;

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $case_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => \@vcfparser_analysis_types,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
            }
        )
    );

    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my @outfile_paths       = @{ $io{out}{file_paths} };
    my @outfile_suffixes    = @{ $io{out}{file_suffixes} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
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

    ## Sort vcf
  INFILE:
    while ( my ( $infile_index, $infile_path ) = each @infile_paths ) {

        ## Prepare for downstream compression and logging
        my $bcftools_suffix = $EMPTY_STR;
        my $metafile_tag    = q{research};
        if ( $infile_index == 1 ) {

            $bcftools_suffix = $DOT . q{selected};
            $metafile_tag    = q{clinical};
        }

        ## Get parameters for sort
        my $sequence_dict_file = catfile( $reference_dir,
            $file_info_href->{human_genome_reference_name_prefix} . $DOT . q{dict} );

        ## Sort variants in vcf format
        picardtools_sortvcf(
            {
                filehandle       => $filehandle,
                infile_paths_ref => [$infile_path],
                java_jar =>
                  catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                memory_allocation    => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
                outfile_path         => $outfile_paths[$infile_index],
                referencefile_path   => $active_parameter_href->{human_genome_reference},
                sequence_dictionary  => $sequence_dict_file,
                temp_directory       => $active_parameter_href->{temp_directory},
            }
        );
        say {$filehandle} $NEWLINE;

        ## Remove variants in hgnc_id list from vcf
        if ( $active_parameter_href->{sv_reformat_remove_genes_file} ) {

            my $filter_metafile_tag = q{sv_reformat_remove_genes_file_research};

            ## Update metafile_tag depending on select or research
            if ( $infile_index == 1 ) {

                $filter_metafile_tag = q{sv_reformat_remove_genes_file_clinical};
            }

            my $filter_file_path = catfile( $reference_dir,
                $active_parameter_href->{sv_reformat_remove_genes_file} );
            my $filter_outfile_path =
                $outfile_path_prefix
              . $UNDERSCORE
              . q{filtered}
              . $outfile_suffixes[$infile_index];
            ## Removes contig_names from contigs array if no male or other found
            gnu_grep(
                {
                    filehandle       => $filehandle,
                    filter_file_path => $filter_file_path,
                    infile_path      => $outfile_paths[$infile_index],
                    invert_match     => 1,
                    stdoutfile_path  => $filter_outfile_path,
                }
            );
            say {$filehandle} $NEWLINE;

            if ( $recipe_mode == 1 ) {

                ## Save filtered file
                set_recipe_metafile_in_sample_info(
                    {
                        recipe_name      => $recipe_name,
                        metafile_tag     => $filter_metafile_tag,
                        path             => $filter_outfile_path,
                        sample_info_href => $sample_info_href,
                    }
                );
            }
        }

        say {$filehandle} q{## Compress};
        ## Reformat variant calling file and index
        bcftools_view_and_index_vcf(
            {
                filehandle          => $filehandle,
                infile_path         => $outfile_paths[$infile_index],
                outfile_path_prefix => $outfile_path_prefix . $bcftools_suffix,
                output_type         => q{z},
            }
        );

        if ( $recipe_mode == 1 ) {

            set_most_complete_vcf(
                {
                    active_parameter_href => $active_parameter_href,
                    path                  => $outfile_paths[$infile_index],
                    recipe_name           => $recipe_name,
                    sample_info_href      => $sample_info_href,
                    vcf_file_key          => q{sv}
                      . $UNDERSCORE
                      . substr( $outfile_suffixes[0], 1 )
                      . $UNDERSCORE . q{file},
                    vcfparser_outfile_counter => $infile_index,
                }
            );
            set_most_complete_vcf(
                {
                    active_parameter_href => $active_parameter_href,
                    path                  => $outfile_paths[$infile_index] . $DOT . q{gz},
                    recipe_name           => $recipe_name,
                    sample_info_href      => $sample_info_href,
                    vcf_file_key          => q{sv}
                      . $UNDERSCORE
                      . substr( $outfile_suffixes[0], 1 )
                      . $UNDERSCORE
                      . q{binary_file},
                    vcfparser_outfile_counter => $infile_index,
                }
            );

            # Save clinical candidate list path
            set_recipe_metafile_in_sample_info(
                {
                    metafile_tag     => $metafile_tag,
                    path             => $outfile_paths[$infile_index] . $DOT . q{gz},
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                }
            );
            set_file_path_to_store(
                {
                    file_tag         => $metafile_tag,
                    file_type        => q{vcf},
                    path             => $outfile_paths[$infile_index] . $DOT . q{gz},
                    sample_info_href => $sample_info_href,
                }
            );
        }
    }
    say {$filehandle} $NEWLINE;

    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe_mode == 1 ) {

        submit_recipe(
            {
                base_command            => $profile_base_command,
                case_id                 => $case_id,
                dependency_method       => q{sample_to_case},
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_chain            => $job_id_chain,
                job_id_href             => $job_id_href,
                job_reservation_name    => $active_parameter_href->{job_reservation_name},
                log                     => $log,
                recipe_file_path        => $recipe_file_path,
                sample_ids_ref          => \@{ $active_parameter_href->{sample_ids} },
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

1;
