package MIP::Recipes::Analysis::Fusion_report;

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
use MIP::Constants qw{ $DOT $LOG_NAME $NEWLINE $PIPE $SINGLE_QUOTE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_fusion_report };

}

sub analysis_fusion_report {

## Function : Generate clinical and research fusion reports from arriba fusion calls
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $job_id_href             => Job id hash {REF}    my ($arg_href) = @_;
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command    ## Flatten argument(s)
##          : $recipe_name             => Recipe name    my $active_parameter_href;
##          : $sample_id               => Sample id    my $file_info_href;
##          : $sample_info_href        => Info on samples and case hash {REF    my $job_id_href;

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

    use MIP::Analysis qw{ get_vcf_parser_analysis_suffix };
    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Language::Perl qw{ perl_nae_oneliners };
    use MIP::Program::Gnu::Software::Gnu_grep qw{ gnu_grep };
    use MIP::Program::Gnu::Coreutils qw{ gnu_uniq };
    use MIP::Program::Arriba qw{ draw_fusions };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{
      get_pedigree_sample_id_attributes
      set_recipe_metafile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => q{arriba_ar},
            stream         => q{out},
        }
    );
    my $infile_name_prefix = $io{out}{file_name_prefix};
    my $infile_path        = $io{out}{file_path};
    my $infile_path_prefix = $io{out}{file_path_prefix};
    my $infile_suffix      = $io{out}{file_suffix};

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    my $sample_display_name = get_pedigree_sample_id_attributes(
        {
            attribute        => q{sample_display_name},
            sample_id        => $sample_id,
            sample_info_href => $sample_info_href,
        }
    );

    if ( $sample_display_name
        and not $active_parameter_href->{fusion_use_sample_id_as_display_name} )
    {

        $infile_name_prefix = $sample_display_name;
    }

    my @report_types = get_vcf_parser_analysis_suffix(
        {
            vcfparser_outfile_count => $active_parameter_href->{fusion_outfile_count},
        }
    );

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $recipe{job_id_chain},
                id               => $sample_id,
                file_info_href   => $file_info_href,
                iterators_ref    => \@report_types,
                file_name_prefix => $infile_name_prefix,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
            }
        )
    );
    my $outdir_path = $io{out}{dir_path};
    my %outfile     = (
        research => {
            file_path   => $io{out}{file_paths}->[0],
            file_suffix => $io{out}{file_suffixes}->[0],
        },
    );
    if ( @report_types > 1 ) {

        $outfile{clinical}{file_path}   = $io{out}{file_paths}->[1];
        $outfile{clinical}{file_suffix} = $io{out}{file_suffixes}->[1];
    }

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

  TAG:
    foreach my $tag ( keys %outfile ) {

        my $fusion_file_path = $infile_path;

        if ( $outfile{$tag}{file_suffix} eq q{.selected.pdf} ) {

            perl_nae_oneliners(
                {
                    filehandle     => $filehandle,
                    oneliner_name  => q{get_gene_panel_hgnc_symbols},
                    print_newline  => 1,
                    stdinfile_path => $active_parameter_href->{fusion_select_file},
                    use_container  => 1,
                }
            );
            print {$filehandle} $PIPE . $SPACE;

            my $gene_list_path = catfile( $outdir_path, q{select_genes.txt} );
            gnu_uniq(
                {
                    filehandle      => $filehandle,
                    infile_path     => q{-},
                    stdoutfile_path => $gene_list_path,
                }
            );
            say {$filehandle} $NEWLINE;

            $fusion_file_path =
              catfile( $outdir_path, $infile_name_prefix . $DOT . q{selected} . $infile_suffix );
            gnu_grep(
                {
                    filehandle       => $filehandle,
                    filter_file_path => $gene_list_path,
                    infile_path      => $infile_path,
                    pattern          => $SINGLE_QUOTE . q{^#gene1} . $SINGLE_QUOTE,
                    stdoutfile_path  => $fusion_file_path,
                    word_regexp      => 1,
                }
            );
            say {$filehandle} $NEWLINE;
        }

        my $bam_file_path = $infile_path_prefix . $DOT . q{bam};
        draw_fusions(
            {
                alignment_file_path      => $bam_file_path,
                annotation_file_path     => $active_parameter_href->{transcript_annotation},
                cytoband_file_path       => $active_parameter_href->{fusion_cytoband_path},
                filehandle               => $filehandle,
                fusion_file_path         => $fusion_file_path,
                outfile_path             => $outfile{$tag}{file_path},
                protein_domain_file_path => $active_parameter_href->{fusion_protein_domain_path},
            }
        );
        say {$filehandle} $NEWLINE;
    }

    ## Close filehandle
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

      TAG:
        foreach my $tag ( keys %outfile ) {

            ## Collect QC metadata info for later use
            set_recipe_metafile_in_sample_info(
                {
                    metafile_tag     => $tag,
                    path             => $outfile{$tag}{file_path},
                    recipe_name      => $recipe_name,
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );
        }

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
