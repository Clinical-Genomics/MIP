package MIP::Recipes::Analysis::Mip_vcfparser;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use List::Util qw{ any };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ %ANALYSIS $DASH $DOT $EMPTY_STR $LOG_NAME $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ analysis_mip_vcfparser analysis_mip_vcfparser_panel analysis_mip_vcfparser_sv_wes analysis_mip_vcfparser_sv_wgs };

}

## Constants
Readonly my $ANNOTATION_DISTANCE    => $ANALYSIS{ANNOTATION_DISTANCE};
Readonly my $ANNOTATION_DISTANCE_MT => $ANALYSIS{ANNOTATION_DISTANCE_MT};
Readonly my $MT_CONTIG_ID_REG_EXP   => q{MT | M | chrM};

sub analysis_mip_vcfparser {

## Function : Vcfparser performs parsing of varianteffectpredictor annotated variants
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $file_path               => File path
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
    my $file_path;
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
        file_path   => { store => \$file_path, strict_type => 1, },
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
    use MIP::Cluster qw{ get_core_number update_memory_allocation };
    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Mip qw{ mip_vcfparser };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Sample_info qw{ set_gene_panel set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    ## Get the io infiles per chain and id
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
    my %infile_path        = %{ $io{in}{file_path_href} };

    my @contigs_size_ordered = @{ $file_info_href->{contigs_size_ordered} };
    my %recipe               = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    my @vcfparser_analysis_types = get_vcf_parser_analysis_suffix(
        {
            vcfparser_outfile_count => $active_parameter_href->{vcfparser_outfile_count},
        }
    );

    ## Add "research" suffixes
    my @set_outfile_name_suffixes = ( keys %infile_path );

    ## Add "select" suffixes
    my @select_file_contigs = @{ $file_info_href->{sorted_select_file_contigs} };
  SUFFIX:
    foreach my $infile_suffix ( keys %infile_path ) {

        if ( any { $_ eq $infile_suffix } @select_file_contigs ) {

            ## Add "selected" suffixes
            push @set_outfile_name_suffixes,
              $infile_suffix . $UNDERSCORE . $vcfparser_analysis_types[1];
        }
    }

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $recipe{job_id_chain},
                id               => $case_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => \@set_outfile_name_suffixes,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
                temp_directory   => $temp_directory,
            }
        )
    );

    my $outdir_path   = $io{out}{dir_path};
    my %outfile_path  = %{ $io{out}{file_path_href} };
    my @outfile_paths = @{ $io{out}{file_paths} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle      = IO::Handle->new();
    my $xargsfilehandle = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    my $core_number = get_core_number(
        {
            max_cores_per_node   => $active_parameter_href->{max_cores_per_node},
            modifier_core_number => scalar @{ $file_info_href->{contigs} },
            recipe_core_number   => $recipe{core_number},
        }
    );
    ## Update memory depending on how many cores that are being used
    my $memory_allocation = update_memory_allocation(
        {
            node_ram_memory           => $active_parameter_href->{node_ram_memory},
            parallel_processes        => $core_number,
            process_memory_allocation => $recipe{memory},
        }
    );

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $case_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $memory_allocation,
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            temp_directory        => $temp_directory,
        }
    );

    ### SHELL:

    ## vcfparser
    say {$filehandle} q{## vcfparser};

    my $xargs_file_path_prefix;

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $core_number,
            filehandle         => $filehandle,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            xargsfilehandle    => $xargsfilehandle,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig (@contigs_size_ordered) {

        ## Get parameters
        my $padding;
        if ( $contig =~ / $MT_CONTIG_ID_REG_EXP /xms ) {

            # Special case for mitochondrial contig annotation
            $padding = $ANNOTATION_DISTANCE_MT;
        }

        my $log_file_path = catfile( $outdir_path, q{vcfparser} . $UNDERSCORE . $contig . q{.log} );
        my $vcfparser_xargs_file_path_prefix = $xargs_file_path_prefix . $DOT . $contig;
        my @select_feature_annotation_columns;
        my $select_file;
        my $select_file_matching_column;
        my $select_outfile;
        if ( $active_parameter_href->{vcfparser_select_file} ) {

            if (
                check_element_exist_hash_of_array(
                    {
                        element  => $contig,
                        hash_ref => $file_info_href,
                        key      => q{select_file_contigs},

                    }
                )
              )
            {
                # List of genes to analyse separately
                $select_file = catfile( $active_parameter_href->{vcfparser_select_file} );

                # Column of HGNC Symbol in SelectFile (-sf)
                $select_file_matching_column =
                  $active_parameter_href->{vcfparser_select_file_matching_column};

                if ( exists $active_parameter_href->{vcfparser_select_feature_annotation_columns} )
                {
                    @select_feature_annotation_columns =
                      @{ $active_parameter_href->{vcfparser_select_feature_annotation_columns} };
                }
                my $select_outfile_suffix_key =
                  $contig . $UNDERSCORE . $vcfparser_analysis_types[1];
                $select_outfile = $outfile_path{$select_outfile_suffix_key};
            }
        }

        mip_vcfparser(
            {
                filehandle    => $xargsfilehandle,
                infile_path   => $infile_path{$contig},
                log_file_path => $log_file_path,
                padding       => $padding,
                parse_vep     => $active_parameter_href->{varianteffectpredictor},
                range_feature_annotation_columns_ref =>
                  \@{ $active_parameter_href->{vcfparser_range_feature_annotation_columns} },
                range_feature_file_path => $active_parameter_href->{vcfparser_range_feature_file},
                select_feature_annotation_columns_ref => \@select_feature_annotation_columns,
                select_feature_file_path              => $select_file,
                select_feature_matching_column        => $select_file_matching_column,
                select_outfile                        => $select_outfile,
                stderrfile_path => $vcfparser_xargs_file_path_prefix . $DOT . q{stderr.txt},
                stdoutfile_path => $outfile_path{$contig},
            }
        );

        say {$xargsfilehandle} $NEWLINE;

        ### Special case: replace all clinical mitochondrial variants with research mitochondrial variants
        _add_all_mt_var_from_research_to_clinical(
            {
                add_all_mt_var => $active_parameter_href->{vcfparser_add_all_mt_var},
                contig         => $contig,
                filehandle     => $filehandle,
                infile_path    => $outfile_path{$contig},
                outfile_path   => $select_outfile,
            }
        );
    }

    if ( $recipe{mode} == 1 ) {

        my %gene_panels = (
            range_file  => q{vcfparser_range_feature_file},
            select_file => q{vcfparser_select_file},
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

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_paths[0],
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );
    }

    close $filehandle or $log->logcroak(q{Could not close $filehandle});
    close $xargsfilehandle
      or $log->logcroak(q{Could not close $xargsfilehandle});

    if ( $recipe{mode} == 1 ) {

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

sub analysis_mip_vcfparser_panel {

## Function : Vcfparser performs parsing of varianteffectpredictor annotated variants
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $file_path               => File path
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
    my $file_path;
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
        file_path   => { store => \$file_path, strict_type => 1, },
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

    use MIP::Analysis qw{ get_vcf_parser_analysis_suffix };
    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::List qw{ get_splitted_lists };
    use MIP::Parameter qw{ get_cache };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bcftools qw{ bcftools_concat bcftools_sort bcftools_view };
    use MIP::Program::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Program::Mip qw{ mip_vcfparser };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{
      set_file_path_to_store
      set_gene_panel
      set_recipe_metafile_in_sample_info
      set_recipe_outfile_in_sample_info
    };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    ## Get the io infiles per chain and id
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

    my $analysis_type = get_cache(
        {
            parameter_href => $parameter_href,
            parameter_name => q{consensus_analysis_type},
        }
    );
    my @vcfparser_analysis_suffixes = get_vcf_parser_analysis_suffix(
        {
            analysis_type           => $analysis_type,
            vcfparser_outfile_count => $active_parameter_href->{vcfparser_outfile_count},
        }
    );

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $recipe{job_id_chain},
                id               => $case_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => \@vcfparser_analysis_suffixes,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
                temp_directory   => $temp_directory,
            }
        )
    );
    my $outdir_path           = $io{out}{dir_path};
    my %outfile_path          = %{ $io{out}{file_path_href} };
    my @outfile_paths         = @{ $io{out}{file_paths} };
    my @outfile_suffixes      = @{ $io{out}{file_suffixes} };
    my $temp_file_path_prefix = $io{temp}{file_path_prefix};
    my @temp_file_suffixes    = map { s/[.]gz\z//xmsr } @outfile_suffixes;

    ## Construct vcfparser program outfiles
    my @temp_ids          = qw{ chr mt };
    my %temp_outfile_path = (
        all => [ map { $temp_file_path_prefix . $DOT . $_ . $temp_file_suffixes[0] } @temp_ids ],
        selected =>
          [ map { $temp_file_path_prefix . $DOT . $_ . $temp_file_suffixes[1] } @temp_ids ],
    );

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ($recipe_file_path) = setup_script(
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

    ### SHELL:

    ## Get parameters
    my $log_file_path = catfile( $outdir_path, q{vcfparser} . q{.log} );
    my @paddings      = ( $ANNOTATION_DISTANCE, $ANNOTATION_DISTANCE_MT );

    # List of genes to analyse separately
    my $select_file_path = $active_parameter_href->{vcfparser_select_file};

    # Column of HGNC Symbol in SelectFile (-sf)
    my $select_file_matching_column =
      $active_parameter_href->{vcfparser_select_file_matching_column};

    my @select_feature_annotation_columns;
    if ( exists $active_parameter_href->{vcfparser_select_feature_annotation_columns} ) {
        @select_feature_annotation_columns =
          @{ $active_parameter_href->{vcfparser_select_feature_annotation_columns} };
    }

    say {$filehandle} q{## Prepare file for vcfparser};
    my $bcftools_outfile_path = catfile( $outdir_path, $infile_name_prefix . q{.vcf.gz} );
    bcftools_sort(
        {
            filehandle     => $filehandle,
            infile_path    => $infile_path,
            outfile_path   => $bcftools_outfile_path,
            output_type    => q{z},
            temp_directory => $temp_directory,
        }
    );
    say {$filehandle} $NEWLINE;

    htslib_tabix(
        {
            filehandle  => $filehandle,
            infile_path => $bcftools_outfile_path,
            preset      => q{vcf},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Get contigs
    my ( $mt_contig_ref, $contigs_ref ) = get_splitted_lists(
        {
            regexp   => qr/M/,
            list_ref => $file_info_href->{primary_contigs},
        }
    );

    ## vcfparser
    say {$filehandle} q{## vcfparser};

    my @all_contigs = ( $contigs_ref, $mt_contig_ref );
  REGIONS_REF:
    while ( my ( $index, $regions_ref ) = each @all_contigs ) {

        say {$filehandle} qq{## parse @{ $regions_ref }};
        bcftools_view(
            {
                filehandle  => $filehandle,
                regions_ref => $regions_ref,
                infile_path => $bcftools_outfile_path,
                output_type => q{v},
            }
        );
        print {$filehandle} $PIPE . $SPACE;

        mip_vcfparser(
            {
                filehandle    => $filehandle,
                infile_path   => $DASH,
                log_file_path => $log_file_path,
                padding       => $paddings[$index],
                parse_vep     => $active_parameter_href->{varianteffectpredictor},
                range_feature_annotation_columns_ref =>
                  \@{ $active_parameter_href->{vcfparser_range_feature_annotation_columns} },
                range_feature_file_path => $active_parameter_href->{vcfparser_range_feature_file},
                select_feature_annotation_columns_ref => \@select_feature_annotation_columns,
                select_feature_file_path              => $select_file_path,
                select_feature_matching_column        => $select_file_matching_column,
                select_outfile                        => $temp_outfile_path{selected}[$index],
                stdoutfile_path                       => $temp_outfile_path{all}[$index],
            }
        );
        say {$filehandle} $NEWLINE;

      TEMP_OUTFILE_TYPE:
        foreach my $temp_outfile_type ( keys %temp_outfile_path ) {

            say {$filehandle} q{## Compress and index};
            htslib_bgzip(
                {
                    filehandle  => $filehandle,
                    infile_path => $temp_outfile_path{$temp_outfile_type}[$index],
                    threads     => $recipe{core_number},
                }
            );
            say {$filehandle} $NEWLINE;

            htslib_tabix(
                {
                    filehandle  => $filehandle,
                    infile_path => $temp_outfile_path{$temp_outfile_type}[$index] . q{.gz},
                    preset      => q{vcf},
                }
            );
            say {$filehandle} $NEWLINE;
        }
    }

    ## Concatenate chr outfiles with mt outfiles
  ANALYSIS_SUFFIXES:
    foreach my $analysis_suffix (@vcfparser_analysis_suffixes) {

        say {$filehandle} qq{## Concatenate and index $analysis_suffix};
        bcftools_concat(
            {
                filehandle       => $filehandle,
                infile_paths_ref =>
                  [ map { $_ . q{.gz} } @{ $temp_outfile_path{$analysis_suffix} } ],
                outfile_path => $outfile_path{$analysis_suffix},
                output_type  => q{z},
                rm_dups      => 0,
                threads      => $recipe{core_number},
            }
        );
        say {$filehandle} $NEWLINE;

        htslib_tabix(
            {
                filehandle  => $filehandle,
                infile_path => $outfile_path{$analysis_suffix},
                preset      => q{vcf},
            }
        );
        say {$filehandle} $NEWLINE;
    }

    close $filehandle or $log->logcroak(q{Could not close $filehandle});

    if ( $recipe{mode} == 1 ) {

        my %gene_panels = (
            range_file  => q{vcfparser_range_feature_file},
            select_file => q{vcfparser_select_file},
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
        if ( $analysis_type eq q{wts} ) {

            my @tags = qw{ research clinical };

          TAG:
            while ( my ( $index, $tag ) = each @tags ) {

                set_file_path_to_store(
                    {
                        format           => q{vcf},
                        id               => $case_id,
                        path             => $outfile_paths[$index],
                        path_index       => $outfile_paths[$index] . q{.tbi},
                        recipe_name      => $recipe_name,
                        sample_info_href => $sample_info_href,
                        tag              => $tag,
                    }
                );
            }
        }

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_paths[0],
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

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

sub analysis_mip_vcfparser_sv_wes {

## Function : Vcfparser performs parsing of varianteffectpredictor annotated wes SV variants
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
    use MIP::Cluster qw{ get_core_number update_memory_allocation };
    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Gatk qw{ gatk_concatenate_variants };
    use MIP::Program::Mip qw{ mip_vcfparser };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
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
    my $infile_path_prefix = $io{in}{file_path_prefix};
    my $infile_suffix      = $io{in}{file_suffix};
    my $infile_path        = $infile_path_prefix . $infile_suffix;

    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my %recipe                  = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    my $xargs_file_path_prefix;

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Always one file for wes
    Readonly my $CORE_NUMBER       => 1;
    Readonly my $MEMORY_ALLOCATION => 1;

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
    my @outfile_suffixes    = @{ $io{out}{file_suffixes} };

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ($recipe_file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $CORE_NUMBER,
            directory_id          => $case_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $MEMORY_ALLOCATION,
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            temp_directory        => $temp_directory,
        }
    );

    ### SHELL:

    ## vcfparser
    say {$filehandle} q{## vcfparser};

    my $log_file_path = catfile( $outdir_path, q{vcfparser.log} );
    my @select_feature_annotation_columns;
    my $select_file;
    my $select_file_matching_column;
    my $select_outfile;
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

        ## Select outfile
        $select_outfile = $outfile_path_prefix . $outfile_suffixes[1];
    }

    mip_vcfparser(
        {
            filehandle           => $filehandle,
            infile_path          => $infile_path,
            log_file_path        => $log_file_path,
            parse_vep            => $active_parameter_href->{sv_varianteffectpredictor},
            per_gene             => $active_parameter_href->{sv_vcfparser_per_gene},
            pli_values_file_path => $active_parameter_href->{vcfparser_pli_score_file},
            range_feature_annotation_columns_ref =>
              \@{ $active_parameter_href->{sv_vcfparser_range_feature_annotation_columns} },
            range_feature_file_path => $active_parameter_href->{sv_vcfparser_range_feature_file},
            select_feature_annotation_columns_ref => \@select_feature_annotation_columns,
            stderrfile_path                       => $recipe_file_path . $DOT . q{stderr.txt},
            stdoutfile_path                       => $outfile_path_prefix . $outfile_suffixes[0],
            select_feature_file_path              => $select_file,
            select_feature_matching_column        => $select_file_matching_column,
            select_outfile                        => $select_outfile,
        }
    );
    say {$filehandle} $NEWLINE;

    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
                path             => $select_outfile,
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

sub analysis_mip_vcfparser_sv_wgs {

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
    use MIP::Program::Gatk qw{ gatk_concatenate_variants };
    use MIP::Program::Mip qw{ mip_vcfparser };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
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

    my $indir_path_prefix  = $io{in}{dir_path_prefix};
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my %infile_path        = %{ $io{in}{file_path_href} };

    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my %recipe                  = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    my $xargs_file_path_prefix;

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
    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my @outfile_suffixes    = @{ $io{out}{file_suffixes} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle      = IO::Handle->new();
    my $xargsfilehandle = IO::Handle->new();

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

    my @contigs = @{ $file_info_href->{contigs} };

    ## vcfparser
    say {$filehandle} q{## vcfparser};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $recipe{core_number},
            filehandle         => $filehandle,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            xargsfilehandle    => $xargsfilehandle,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig (@contigs) {

        my $padding;
        if ( $contig =~ / $MT_CONTIG_ID_REG_EXP /sxm ) {

            ## Special case for mitochondrial contig annotation
            $padding = $ANNOTATION_DISTANCE_MT;
        }

        my $log_file_path = catfile( $outdir_path, q{vcfparser} . $UNDERSCORE . $contig . q{.log} );
        my $vcfparser_outfile_path = $outfile_path_prefix . $DOT . $contig . $outfile_suffix;
        my $vcfparser_xargs_file_path_prefix = $xargs_file_path_prefix . $DOT . $contig;
        my @select_feature_annotation_columns;
        my $select_file;
        my $select_file_matching_column;
        my $select_outfile;

        if ( $active_parameter_href->{sv_vcfparser_select_file} ) {

            if (
                check_element_exist_hash_of_array(
                    {
                        element  => $contig,
                        hash_ref => $file_info_href,
                        key      => q{select_file_contigs},
                    }
                )
              )
            {

                ## List of genes to analyse separately
                $select_file = catfile( $active_parameter_href->{sv_vcfparser_select_file} );

                ## Column of HGNC Symbol in select file ("-sf")
                $select_file_matching_column =
                  $active_parameter_href->{sv_vcfparser_select_file_matching_column};

                if (
                    exists $active_parameter_href->{sv_vcfparser_select_feature_annotation_columns}
                  )
                {

                    @select_feature_annotation_columns =
                      @{ $active_parameter_href->{sv_vcfparser_select_feature_annotation_columns} };
                }

                ## Select outfile
                $select_outfile = $outfile_path_prefix . $DOT . $contig . $outfile_suffixes[1];
            }
        }

        mip_vcfparser(
            {
                filehandle           => $xargsfilehandle,
                infile_path          => $infile_path{$contig},
                log_file_path        => $log_file_path,
                padding              => $padding,
                parse_vep            => $active_parameter_href->{sv_varianteffectpredictor},
                per_gene             => $active_parameter_href->{sv_vcfparser_per_gene},
                pli_values_file_path => $active_parameter_href->{vcfparser_pli_score_file},
                range_feature_annotation_columns_ref =>
                  \@{ $active_parameter_href->{sv_vcfparser_range_feature_annotation_columns} },
                range_feature_file_path =>
                  $active_parameter_href->{sv_vcfparser_range_feature_file},
                select_feature_annotation_columns_ref => \@select_feature_annotation_columns,
                stderrfile_path => $vcfparser_xargs_file_path_prefix . $DOT . q{stderr.txt},
                stdoutfile_path => $vcfparser_outfile_path,
                select_feature_file_path       => $select_file,
                select_feature_matching_column => $select_file_matching_column,
                select_outfile                 => $select_outfile,
            }
        );

        say {$xargsfilehandle} $NEWLINE;

        ### Special case: replace all clinical mitochondrial variants with research mitochondrial variants
        _add_all_mt_var_from_research_to_clinical(
            {
                add_all_mt_var => $active_parameter_href->{sv_vcfparser_add_all_mt_var},
                contig         => $contig,
                filehandle     => $filehandle,
                infile_path    => $vcfparser_outfile_path,
                outfile_path   => $select_outfile,
            }
        );
    }

  ANALYSIS_SUFFIXES:
    foreach my $analysis_suffix (@outfile_suffixes) {

        my @concat_contigs = @contigs;

        ## Update contigs list using select file contigs
        if ( $analysis_suffix eq q{.selected.vcf} ) {

            @concat_contigs = @{ $file_info_href->{select_file_contigs} };
        }
        gatk_concatenate_variants(
            {
                active_parameter_href => $active_parameter_href,
                elements_ref          => \@concat_contigs,
                continue              => 1,
                filehandle            => $filehandle,
                infile_postfix        => $analysis_suffix,
                infile_prefix         => $outfile_path_prefix,
                outfile_path_prefix   => $outfile_path_prefix,
                outfile_suffix        => $analysis_suffix,
            }
        );
    }
    say {$filehandle} q{wait} . $NEWLINE;

    close $filehandle or $log->logcroak(q{Could not close filehandle});
    close $xargsfilehandle
      or $log->logcroak(q{Could not close xargsfilehandle});

    if ( $recipe{mode} == 1 ) {

        my $outfile_sample_info_prefix =
          $outfile_path_prefix . $DOT . $contigs[0] . $outfile_suffix;

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
                path             => $outfile_sample_info_prefix,
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

sub _add_all_mt_var_from_research_to_clinical {

## Function : Replace all clinical mitochondrial variants with research mitochondrial variants
## Returns  :
## Arguments: $add_all_mt_var => Add all mt variants switch
##          : $contig         => Contig {REF}
##          : $filehandle     => Filehandle to write to
##          : $infile_path    => File to copy from
##          : $outfile_path   => File to replace

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contig;
    my $filehandle;
    my $infile_path;
    my $outfile_path;

    ## Default(s)
    my $add_all_mt_var;

    my $tmpl = {
        add_all_mt_var => {
            allow       => qr{ \A\d+\z }sxm,
            default     => 1,
            store       => \$add_all_mt_var,
            strict_type => 1,
        },
        contig => {
            defined     => 1,
            required    => 1,
            store       => \$contig,
            strict_type => 1,
        },
        filehandle => {
            defined  => 1,
            required => 1,
            store    => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path => {
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Gnu::Coreutils qw{ gnu_cp };

    return if ( not defined $outfile_path );

    if ( $add_all_mt_var and $contig =~ / $MT_CONTIG_ID_REG_EXP /sxm ) {

        say {$filehandle} q{## Replacing clinical MT variants with research MT variants};
        gnu_cp(
            {
                filehandle   => $filehandle,
                infile_path  => $infile_path,
                outfile_path => $outfile_path,
            }
        );
    }
    return;
}

1;
