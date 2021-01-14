package MIP::Recipes::Analysis::Glnexus;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile devnull };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $NEWLINE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_glnexus };

}

sub analysis_glnexus {

## Function : Merges gvcfs from DeepVariant to generate a cohort vcf.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Recipe name
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
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_io_files };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Program::Bcftools qw{ bcftools_view_and_index_vcf };
    use MIP::Program::Gnu::Coreutils qw{ gnu_cp };
    use MIP::Program::Glnexus qw{ glnexus_merge };
    use MIP::Program::Htslib qw{ htslib_bgzip };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my %consensus_analysis_type_map =
      ( MIXED => q{WGS}, PANEL => q{WES}, WGS => q{WGS}, WES => q{WES} );

    my $consensus_analysis_type =
      $consensus_analysis_type_map{ uc $parameter_href->{cache}{consensus_analysis_type} };

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number = $recipe{core_number};
    my $time        = $recipe{time};

    ## Get the io infiles per chain and id
    my @genotype_infile_paths;
    my @genotype_infile_path_prefixes;

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Get the io infiles per chain and id
        my %sample_io = get_io_files(
            {
                id             => $sample_id,
                file_info_href => $file_info_href,
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
                stream         => q{in},
            }
        );
        push @genotype_infile_paths,         $sample_io{in}{file_path};
        push @genotype_infile_path_prefixes, $sample_io{in}{file_path_prefix};
    }

    my %io = parse_io_outfiles(
        {
            chain_id               => $recipe{job_id_chain},
            id                     => $case_id,
            file_info_href         => $file_info_href,
            file_name_prefixes_ref => [$case_id],
            outdata_dir            => $active_parameter_href->{outdata_dir},
            parameter_href         => $parameter_href,
            recipe_name            => $recipe_name,
        }
    );

    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_path =
      catdir( $active_parameter_href->{temp_directory}, $io{out}{file_name_prefix} . q{.vcf} );

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe{core_number},
            directory_id          => $case_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe{memory},
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
        }
    );

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    my $config_type = q{DeepVariant} . $consensus_analysis_type;

    if ( scalar @{ $active_parameter_href->{sample_ids} } > 1 ) {

        ## Glnexus
        say {$filehandle} q{## Glnexus};

        glnexus_merge(
            {
                config           => $config_type,
                dir              => catdir( $active_parameter_href->{temp_directory}, q{glnexus} ),
                filehandle       => $filehandle,
                infile_paths_ref => \@genotype_infile_paths,
                stdoutfile_path  => $outfile_path,
            }
        );
        say {$filehandle} $NEWLINE;

        say {$filehandle} q{## View};

        bcftools_view_and_index_vcf(
            {
                filehandle          => $filehandle,
                index_type          => q{tbi},
                infile_path         => $outfile_path,
                outfile_path_prefix => $outfile_path_prefix,
                output_type         => q{z},
                threads             => $core_number,
            }
        );
    }
    else {

        say {$filehandle} q{## Single sample - copy deepvariant vcf output and index};

      FILE_SUFFIX:
        foreach my $dv_file_name_suffix (qw { .vcf.gz .vcf.gz.tbi }) {
            gnu_cp {
                filehandle   => $filehandle,
                infile_path  => $genotype_infile_path_prefixes[0] . $dv_file_name_suffix,
                outfile_path => $outfile_path_prefix . $dv_file_name_suffix,
            };
            say {$filehandle} $NEWLINE;
        }
    }
    ## Close filehandle
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_case},
                log                               => $log,
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
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
