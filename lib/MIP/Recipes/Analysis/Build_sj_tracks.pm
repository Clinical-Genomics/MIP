package MIP::Recipes::Analysis::Build_sj_tracks;

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
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $DOT $LOG_NAME $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_build_sj_tracks };

}

sub analysis_build_sj_tracks {

## Function : Build splice junction tracks for IGV
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $case_id               => Family id
##          : $file_info_href        => File_info hash {REF}
##          : $job_id_href           => Job id hash {REF}
##          : $parameter_href        => Parameter hash {REF}
##          : $profile_base_command  => Submission profile base command
##          : $recipe_name           => Recipe name
##          : $sample_id             => Sample id
##          : $sample_info_href      => Info on samples and case hash {REF}

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
    use MIP::Language::Perl qw{ perl_nae_oneliners };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Program::Ucsc qw{ ucsc_wig_to_big_wig };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Reference qw{ write_contigs_size_file };
    use MIP::Sample_info
      qw{ set_file_path_to_store  set_recipe_metafile_in_sample_info set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my $star_sj_file_suffix  = q{.SJ.out.tab};
    my $star_wig_file_suffix = q{.Signal.UniqueMultiple.str1.out.wig};

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path_prefix = $io{in}{file_path_prefix};

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $recipe{job_id_chain},
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => [$infile_name_prefix],
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
    my $bed_outfile_path    = $outfile_path_prefix . q{.bed.gz};

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

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    say {$filehandle} q{## Create chromosome name and size file};
    my $contigs_size_file_path = catfile( $outdir_path, q{contigs_size_file} . $DOT . q{tsv} );
    write_contigs_size_file(
        {
            fai_file_path => $active_parameter_href->{human_genome_reference} . $DOT . q{fai},
            outfile_path  => $contigs_size_file_path,
        }
    );

    say {$filehandle} q{## Create big wig file};
    ucsc_wig_to_big_wig(
        {
            clip                   => 1,
            contigs_size_file_path => $contigs_size_file_path,
            filehandle             => $filehandle,
            infile_path            => $infile_path_prefix . $star_wig_file_suffix,
            outfile_path           => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Convert SJ file to bed};
    perl_nae_oneliners(
        {
            filehandle     => $filehandle,
            oneliner_name  => q{star_sj_tab_to_bed},
            print_newline  => 1,
            stdinfile_path => $infile_path_prefix . $star_sj_file_suffix,
            use_container  => 1,
        }

    );
    print {$filehandle} $PIPE . $SPACE;

    htslib_bgzip(
        {
            filehandle      => $filehandle,
            infile_path     => q{-},
            stdoutfile_path => $bed_outfile_path,
            write_to_stdout => 1,
        }
    );
    say {$filehandle} $NEWLINE;

    htslib_tabix(
        {
            filehandle  => $filehandle,
            infile_path => $bed_outfile_path,
            preset      => q{bed},
        }
    );
    say {$filehandle} $NEWLINE;

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
                format           => q{bigwig},
                id               => $sample_id,
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
                tag              => q{coverage},
            }
        );

        set_file_path_to_store(
            {
                format           => q{bed},
                id               => $sample_id,
                path             => $bed_outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
                tag              => q{junction},
            }
        );

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_island},
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
