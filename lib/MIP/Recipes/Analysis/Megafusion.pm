package MIP::Recipes::Analysis::Megafusion;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $DOT $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_megafusion };

}

sub analysis_megafusion {

## Function : Convert tsv files from fusion callers to vcf format and combine into one
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Recipe name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_id;
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
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
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

    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bcftools qw{bcftools_view_and_index_vcf };
    use MIP::Program::Megafusion qw{ megafusion };
    use MIP::Program::Svdb qw{ svdb_merge };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info
      qw{ set_file_path_to_store set_recipe_metafile_in_sample_info set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %fusion_file_path;

  FUSION_CALLER_RECIPE:
    foreach my $fusion_caller_recipe ( @{ $active_parameter_href->{megafusion_callers} } ) {

        next FUSION_CALLER_RECIPE if not $active_parameter_href->{$fusion_caller_recipe};

        my %io = get_io_files(
            {
                id             => $sample_id,
                file_info_href => $file_info_href,
                parameter_href => $parameter_href,
                recipe_name    => $fusion_caller_recipe,
                stream         => q{out},
            }
        );
        $fusion_file_path{$fusion_caller_recipe}{infile_path} = $io{out}{file_path};
        $fusion_file_path{$fusion_caller_recipe}{file_name_prefix} =
          $io{out}{file_name_prefixes}[0];
    }

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    my %io = (
        parse_io_outfiles(
            {
                chain_id               => $recipe{job_id_chain},
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => [$sample_id],
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );

    my $outfile_name_prefix = $io{out}{file_name_prefix};
    my $outfile_path        = $io{out}{file_path};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outdir_path         = $io{out}{dir_path};
    my $outfile_suffix      = $io{out}{file_suffix};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe{core_number},
            directory_id          => $sample_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe{memory},
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
        }
    );

    my %config_path = (
        arriba_ar   => $active_parameter_href->{megafusion_arriba_config},
        star_fusion => $active_parameter_href->{megafusion_star_fusion_config},
    );

    my @fusion_vcfs;

    ### SHELL:
    say {$filehandle} q{## } . $recipe_name;

  FUSION_RECIPE:
    foreach my $fusion_recipe ( keys %fusion_file_path ) {

        my $megafusion_outfile_path =
          $outdir_path . $fusion_file_path{$fusion_recipe}{file_name_prefix} . $outfile_suffix;
        megafusion(
            {
                config_file_path => $config_path{$fusion_recipe},
                filehandle       => $filehandle,
                infile_path      => $fusion_file_path{$fusion_recipe}{infile_path},
                sample_id        => $sample_id,
                stdoutfile_path  => $megafusion_outfile_path,
            }
        );
        say {$filehandle} $NEWLINE;

        push @fusion_vcfs, $megafusion_outfile_path;
    }

    say {$filehandle} q{## Merge fusion calls};
    svdb_merge(
        {
            filehandle       => $filehandle,
            infile_paths_ref => \@fusion_vcfs,
            stdoutfile_path  => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Compress vcf for storing};
    bcftools_view_and_index_vcf(
        {
            filehandle          => $filehandle,
            index               => 1,
            index_type          => q{tbi},
            infile_path         => $outfile_path,
            outfile_path_prefix => $outfile_path_prefix,
            output_type         => q{z},
            threads             => $recipe{core_number},
        }
    );

    ## Close filehandles
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => $outfile_name_prefix,
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        set_file_path_to_store(
            {
                format           => q{vcf},
                id               => $sample_id,
                path             => $outfile_path_prefix . $DOT . q{vcf.gz},
                path_index       => $outfile_path_prefix . $DOT . q{vcf.gz.tbi},
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_sample},
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
                log                               => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_id          => $sample_id,
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

1;
