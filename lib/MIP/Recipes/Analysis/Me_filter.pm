package MIP::Recipes::Analysis::Me_filter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile devnull };
use List::Util qw{ first};
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ %ANALYSIS $DOT $DOUBLE_QUOTE $LOG_NAME $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_me_filter };

}

## Constants
Readonly my $ANNOTATION_DISTANCE    => $ANALYSIS{ANNOTATION_DISTANCE};
Readonly my $ANNOTATION_DISTANCE_MT => $ANALYSIS{ANNOTATION_DISTANCE_MT};
Readonly my $MT_CONTIG_ID_REGEXP    => q{MT | M | chrM};

sub analysis_me_filter {

## Function : Vcfparser performs parsing of varianteffectpredictor annotated ws SV variants
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $filehandle              => Sbatch filehandle to write to
##          : $file_info_href          => File info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

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
    my $xargs_file_counter;

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
        xargs_file_counter => {
            allow       => qr{ \A\d+\z }xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Analysis qw{ get_vcf_parser_analysis_suffix };
    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::List qw{ check_element_exist_hash_of_array };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bcftools qw{ bcftools_concat bcftools_sort bcftools_view };
    use MIP::Program::Htslib qw{ htslib_tabix };
    use MIP::Program::Mip qw{ mip_vcfparser };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info
      qw{ set_file_path_to_store set_gene_panel set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script};

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
            temp_directory => $temp_directory,
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path        = $io{in}{file_path};

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    my @vcfparser_analysis_types = get_vcf_parser_analysis_suffix(
        {
            vcfparser_outfile_count => $active_parameter_href->{sv_vcfparser_outfile_count},
        }
    );

    ## Set and get the io files per chain, id and stream
    my @set_outfile_name_prefixes =
      map { $infile_name_prefix . $_ } @vcfparser_analysis_types;
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $recipe{job_id_chain},
                id               => $case_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => \@vcfparser_analysis_types,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
                temp_directory   => $temp_directory,
            }
        )
    );
    my $outdir_path         = $io{out}{dir_path};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my @outfile_paths       = @{ $io{out}{file_paths} };
    my %outfile_path        = %{ $io{out}{file_path_href} };
    my @outfile_suffixes    = @{ $io{out}{file_suffixes} };

    ## Filehandles
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
            temp_directory        => $temp_directory,
        }
    );

    my $mt_contig = first { / $MT_CONTIG_ID_REGEXP /xms } @{ $file_info_href->{contigs} };
    my @select_feature_annotation_columns;
    my $select_file;
    my $select_file_matching_column;
    my $non_mt_outfile_path        = $outfile_path_prefix . q{_non_MT} . $outfile_suffixes[0];
    my $select_non_mt_outfile_path = $outfile_path_prefix . q{_non_MT} . $outfile_suffixes[1];
    my $mt_outfile_path            = $outfile_path_prefix . q{_MT} . $outfile_suffixes[0];
    my $select_mt_outfile_path     = $outfile_path_prefix . q{_MT} . $outfile_suffixes[1];

    if ( $active_parameter_href->{sv_vcfparser_select_file} ) {

        ## List of genes to analyse separately
        $select_file = catfile( $active_parameter_href->{sv_vcfparser_select_file} );

        ## Column of HGNC Symbol in select file ("-sf")
        $select_file_matching_column =
          $active_parameter_href->{sv_vcfparser_select_file_matching_column};

        if ( exists $active_parameter_href->{sv_vcfparser_select_feature_annotation_columns} ) {

            @select_feature_annotation_columns =
              @{ $active_parameter_href->{sv_vcfparser_select_feature_annotation_columns} };
        }
    }

    my $exclude_filter = _build_bcftools_frequency_filter(
        {
            frequency_threshold         => $active_parameter_href->{me_filter_frequency_threshold},
            me_annotate_query_info_href => $active_parameter_href->{me_annotate_query_files},
        }
    );

    say {$filehandle} q{## Vcfparser non MT contigs};
    bcftools_view(
        {
            exclude     => $exclude_filter,
            filehandle  => $filehandle,
            infile_path => $infile_path,
            targets     => q{^} . $mt_contig,
        }
    );
    print {$filehandle} $PIPE . $SPACE;

    mip_vcfparser(
        {
            filehandle           => $filehandle,
            infile_path          => catfile( dirname( devnull() ), q{stdin} ),
            log_file_path        => catfile( $outdir_path,         q{vcfparser} . q{.log} ),
            padding              => $ANNOTATION_DISTANCE,
            parse_vep            => $active_parameter_href->{sv_varianteffectpredictor},
            per_gene             => $active_parameter_href->{sv_vcfparser_per_gene},
            pli_values_file_path => $active_parameter_href->{vcfparser_pli_score_file},
            range_feature_annotation_columns_ref =>
              \@{ $active_parameter_href->{sv_vcfparser_range_feature_annotation_columns} },
            range_feature_file_path => $active_parameter_href->{sv_vcfparser_range_feature_file},
            select_feature_annotation_columns_ref => \@select_feature_annotation_columns,
            stdoutfile_path                       => $non_mt_outfile_path,
            select_feature_file_path              => $select_file,
            select_feature_matching_column        => $select_file_matching_column,
            select_outfile                        => $select_non_mt_outfile_path,
            variant_type                          => q{sv},
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Vcfparser MT};
    bcftools_view(
        {
            exclude     => $exclude_filter,
            filehandle  => $filehandle,
            infile_path => $infile_path,
            regions_ref => [$mt_contig],
        }
    );
    print {$filehandle} $PIPE . $SPACE;

    mip_vcfparser(
        {
            filehandle           => $filehandle,
            infile_path          => catfile( dirname( devnull() ), q{stdin} ),
            log_file_path        => catfile( $outdir_path,         q{vcfparser} . q{.log} ),
            padding              => $ANNOTATION_DISTANCE_MT,
            parse_vep            => $active_parameter_href->{me_varianteffectpredictor},
            per_gene             => $active_parameter_href->{sv_vcfparser_per_gene},
            pli_values_file_path => $active_parameter_href->{vcfparser_pli_score_file},
            range_feature_annotation_columns_ref =>
              \@{ $active_parameter_href->{sv_vcfparser_range_feature_annotation_columns} },
            range_feature_file_path => $active_parameter_href->{sv_vcfparser_range_feature_file},
            select_feature_annotation_columns_ref => \@select_feature_annotation_columns,
            stdoutfile_path                       => $mt_outfile_path,
            select_feature_file_path              => $select_mt_outfile_path,
            select_feature_matching_column        => $select_file_matching_column,
            select_outfile                        => $select_mt_outfile_path,
            variant_type                          => q{sv},
        }
    );
    say {$filehandle} $NEWLINE;

    ### Special case: replace all clinical mitochondrial variants with research mitochondrial variants
    $select_mt_outfile_path =
        $active_parameter_href->{sv_vcfparser_add_all_mt_var}
      ? $mt_outfile_path
      : $select_mt_outfile_path;

    ## Concatenate MT variants with the rest
    my @file_sets = (
        {
            files_to_concat_ref => [ $non_mt_outfile_path, $mt_outfile_path ],
            outfile             => $outfile_paths[0],
        },
        {
            files_to_concat_ref => [ $select_non_mt_outfile_path, $select_mt_outfile_path ],
            outfile             => $outfile_paths[1],
        },
    );

    foreach my $outfile_set (@file_sets) {

        bcftools_concat(
            {
                allow_overlaps   => 1,
                filehandle       => $filehandle,
                infile_paths_ref => $outfile_set->{files_to_concat_ref},
                output_type      => q{u},
                threads          => $recipe{core_number},
            }
        );
        print {$filehandle} $PIPE . $SPACE;

        bcftools_sort(
            {
                filehandle     => $filehandle,
                output_type    => q{z},
                outfile_path   => $outfile_set->{outfile},
                temp_directory => $temp_directory,
            }
        );
        say {$filehandle} $NEWLINE;

        htslib_tabix(
            {
                filehandle  => $filehandle,
                force       => 1,
                infile_path => $outfile_set->{outfile},
            }
        );
        say {$filehandle} $NEWLINE;
    }

    close $filehandle;

    if ( $recipe{mode} == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
                path             => $outfile_paths[0],
            }
        );

        my %gene_panels = (
            range_file  => q{sv_vcfparser_range_feature_file},
            select_file => q{sv_vcfparser_select_file},
        );
      GENE_PANEL:
        while ( my ( $gene_panel_key, $gene_panel_file ) = each %gene_panels ) {

            ## Collect databases(s) from a potentially merged gene panel file and adds them to sample_info
            set_gene_panel(
                {
                    aggregate_gene_panel_file => $active_parameter_href->{$gene_panel_file},
                    aggregate_gene_panels_key => $gene_panel_key,
                    recipe_name               => $recipe_name,
                    sample_info_href          => $sample_info_href,
                }
            );
        }

      OUTFILE:
        while ( my ( $file_key, $outfile ) = each %outfile_path ) {

            my $metafile_tag = $file_key =~ m/selected/xms ? q{clinical} : q{research};

            set_file_path_to_store(
                {
                    format           => q{vcf},
                    id               => $case_id,
                    path             => $outfile,
                    path_index       => $outfile . $DOT . q{tbi},
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                    tag              => $metafile_tag,
                }
            );
        }

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_case},
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
                log                               => $log,
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

sub _build_bcftools_frequency_filter {

## Function : Build the frequency filter command
## Returns  :
## Arguments: $frequency_threshold         => Frequency annotation to use in filtering
##          : $me_annotate_query_info_href => SVDB query files used to annotate ME {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $frequency_threshold;
    my $me_annotate_query_info_href;

    my $tmpl = {
        frequency_threshold => {
            required    => 1,
            store       => \$frequency_threshold,
            strict_type => 1,
        },
        me_annotate_query_info_href => {
            default  => {},
            required => 1,
            store    => \$me_annotate_query_info_href,
        }
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Skip if no threshold is set
    return if ( not $frequency_threshold );

    ## Find annotations to use for filtering
    my @frequency_filters;
    my $frequency_expression = $SPACE . q{>} . $SPACE . $frequency_threshold;

  ANNOTATION:
    while ( my ( $annotation_file, $annotation_tags ) = each %{$me_annotate_query_info_href} ) {

        my ( $tag_prefix, $fq_suffix, $cnt_suffix, $file_fq_tag, $file_cnt_tag, $use_in_filter ) =
          split /[|]/sxm, $annotation_tags;

        next ANNOTATION if ( not $use_in_filter );

        push @frequency_filters, q{INFO/} . $tag_prefix . $file_fq_tag . $frequency_expression;
    }

    ## build filter
    my $exclude_filter =
      $DOUBLE_QUOTE . join( $SPACE . $PIPE . $SPACE, @frequency_filters ) . $DOUBLE_QUOTE;

    return $exclude_filter;
}

1;
