package MIP::Recipes::Analysis::Gatk_gathervcfs;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
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
use MIP::Constants qw{ $ASTERISK $DOT $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.09;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_gatk_gathervcfs };

}

sub analysis_gatk_gathervcfs {

## Function : Gather VCFs produced after gatk_genotypegvcfs done per contig.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}

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
    use MIP::Program::Gnu::Coreutils qw(gnu_mv);
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bcftools qw{ bcftools_view_and_index_vcf };
    use MIP::Program::Gatk qw{ gatk_gathervcfscloud gatk_selectvariants };
    use MIP::Sample_info
      qw{ set_processing_metafile_in_sample_info set_recipe_outfile_in_sample_info };
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
    my %infile_path        = %{ $io{in}{file_path_href} };

    my $job_id_chain = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my $referencefile_path      = $active_parameter_href->{human_genome_reference};
    my $recipe_mode             = $active_parameter_href->{$recipe_name};
    my %recipe_resource         = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $case_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => [$infile_name_prefix],
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my $outfile_path        = $io{out}{file_path};

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
            memory_allocation               => $recipe_resource{memory},
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
            temp_directory                  => $temp_directory,
        }
    );

    ## Gather vcf files
    my @gatk_infile_paths =
      map { $infile_path{$_} } @{ $file_info_href->{contigs} };

    ## GATK GatherVcfsCloud
    gatk_gathervcfscloud(
        {
            filehandle           => $filehandle,
            ignore_safety_checks => 0,
            infile_paths_ref     => \@gatk_infile_paths,
            memory_allocation    => q{Xmx4G},
            outfile_path         => $outfile_path,
            temp_directory       => $temp_directory,
            verbosity            => $active_parameter_href->{gatk_logging_level},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Produce a bcf compressed and index from vcf
    if ( $active_parameter_href->{gatk_gathervcfs_bcf_file} ) {

        # Exome or panel analysis
        if ( $consensus_analysis_type =~ /wes|panel/xms ) {

            say {$filehandle} q{### Remove extra reference samples};
            say {$filehandle} q{## GATK SelectVariants};
            gatk_selectvariants(
                {
                    filehandle  => $filehandle,
                    infile_path => $outfile_path,
                    java_use_large_pages =>
                      $active_parameter_href->{java_use_large_pages},
                    memory_allocation => q{Xmx2g},
                    outfile_path      => $outfile_path_prefix
                      . $UNDERSCORE
                      . q{incnonvariantloci}
                      . $outfile_suffix,
                    referencefile_path => $referencefile_path,
                    sample_names_ref   => \@{ $active_parameter_href->{sample_ids} },
                    temp_directory     => $temp_directory,
                    verbosity          => $active_parameter_href->{gatk_logging_level},
                }
            );
            say {$filehandle} $NEWLINE;

            ## Move to original filename
            gnu_mv(
                {
                    filehandle  => $filehandle,
                    infile_path => $outfile_path_prefix
                      . $UNDERSCORE
                      . q{incnonvariantloci}
                      . $outfile_suffix,
                    outfile_path => $outfile_path,
                }
            );
            say {$filehandle} $NEWLINE;
        }

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

        if ( $active_parameter_href->{gatk_gathervcfs_bcf_file} ) {

            my $gbcf_file_path = $outfile_path_prefix . $DOT . q{bcf};

            set_processing_metafile_in_sample_info(
                {
                    metafile_tag     => q{gbcf_file},
                    path             => $gbcf_file_path,
                    sample_info_href => $sample_info_href,
                }
            );
        }

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

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
                recipe_file_path   => $recipe_file_path,
                sample_ids_ref     => \@{ $active_parameter_href->{sample_ids} },
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

1;
