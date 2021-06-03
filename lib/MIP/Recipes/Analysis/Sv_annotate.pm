package MIP::Recipes::Analysis::Sv_annotate;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catfile splitpath };
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
  qw{ $ASTERISK $COLON $DASH $DOT $DOUBLE_QUOTE $EMPTY_STR $LOG_NAME $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_sv_annotate };

}

sub analysis_sv_annotate {

## Function : Annotate structural variants
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $reference_dir           => MIP reference directory
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory {REF}

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
        reference_dir => {
            default     => $arg_href->{active_parameter_href}{reference_dir},
            store       => \$reference_dir,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Program::Gnu::Coreutils qw(gnu_mv);
    use MIP::Io::Read qw{ read_from_file };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bcftools
      qw{ bcftools_annotate bcftools_filter bcftools_view bcftools_view_and_index_vcf };
    use MIP::Program::Genmod qw{ genmod_annotate };
    use MIP::Program::Picardtools qw{ sort_vcf };
    use MIP::Program::Svdb qw{ svdb_query };
    use MIP::Program::Vcfanno qw{ vcfanno };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
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
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path_prefix = $io{in}{file_path_prefix};
    my $infile_suffix      = $io{in}{file_suffix};
    my $infile_path        = $infile_path_prefix . $infile_suffix;

    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my $sequence_dict_file      = catfile( $reference_dir,
        $file_info_href->{human_genome_reference_name_prefix} . $DOT . q{dict} );
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
                id                     => $case_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => [$infile_name_prefix],
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );

    my $outfile_path        = $io{out}{file_path};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};

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
            temp_directory        => $temp_directory,
        }
    );

    ## Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) =
      splitpath( $recipe_info_path . $DOT . q{stderr.txt} );
    my $stderrfile_path = $recipe_info_path . $DOT . q{stderr.txt};

    ### SHELL:

    ## Alternative file tag
    my $alt_file_tag = $EMPTY_STR;

    ## Store annotations to use in sv filtering
    my @svdb_query_annotations;

    if ( $active_parameter_href->{sv_svdb_query} ) {

        ## Set for first infile
        my $svdb_infile_path = $infile_path;

        ## Update alternative ending
        $alt_file_tag .= $UNDERSCORE . q{svdbq};

        ## Ensure correct infile
        my $annotation_file_counter = 0;

        ## Ensure correct outfiles
        my $outfile_tracker = 0;

      QUERIES:
        while ( my ( $query_db_file, $query_db_tag_info ) =
            each %{ $active_parameter_href->{sv_svdb_query_db_files} } )
        {

            if ($annotation_file_counter) {

                $svdb_infile_path =
                  $outfile_path_prefix . $alt_file_tag . $outfile_suffix . $DOT . $outfile_tracker;

                ## Increment now that infile has been set
                $outfile_tracker++;
            }
            ## Get parameters
# Split query_db_tag to decide svdb input query fields
# FORMAT: filename|OUT_FREQUENCY_INFO_KEY|OUT_ALLELE_COUNT_INFO_KEY|IN_FREQUENCY_INFO_KEY|IN_ALLELE_COUNT_INFO_KEY|USE_IN_FREQUENCY_FILTER
            my ( $query_db_tag, $out_frequency_tag_suffix, $out_allele_count_tag_suffix,
                $in_frequency_tag, $in_allele_count_tag, $is_frequency )
              = split /[|]/sxm, $query_db_tag_info;
            my $out_frequency_tag    = $out_frequency_tag_suffix    ||= $EMPTY_STR;
            my $out_allele_count_tag = $out_allele_count_tag_suffix ||= $EMPTY_STR;

            ## Add annotations to filter downstream
            if ( $is_frequency and $out_frequency_tag ) {

                push @svdb_query_annotations, $query_db_tag . $out_frequency_tag;
            }

            svdb_query(
                {
                    bnd_distance         => 25_000,
                    dbfile_path          => $query_db_file,
                    filehandle           => $filehandle,
                    infile_path          => $svdb_infile_path,
                    in_frequency_tag     => $in_frequency_tag,
                    in_allele_count_tag  => $in_allele_count_tag,
                    out_frequency_tag    => $query_db_tag . $out_frequency_tag,
                    out_allele_count_tag => $query_db_tag . $out_allele_count_tag,
                    stdoutfile_path      => $outfile_path_prefix
                      . $alt_file_tag
                      . $outfile_suffix
                      . $DOT
                      . $outfile_tracker,
                    overlap => 0.8,
                }
            );
            say {$filehandle} $NEWLINE;
            $annotation_file_counter++;
        }

        ## Rename to remove outfile_tracker
        gnu_mv(
            {
                filehandle  => $filehandle,
                infile_path => $outfile_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix
                  . $DOT
                  . $outfile_tracker,
                outfile_path => $outfile_path_prefix . $alt_file_tag . $outfile_suffix,
            }
        );
        say {$filehandle} $NEWLINE;
    }

    ## Alternative file tag
    my $outfile_alt_file_tag = $alt_file_tag . $UNDERSCORE . q{sorted};

    ## Writes sbatch code to supplied filehandle to sort variants in vcf format
    sort_vcf(
        {
            active_parameter_href => $active_parameter_href,
            filehandle            => $filehandle,
            infile_paths_ref      => [ $outfile_path_prefix . $alt_file_tag . $outfile_suffix ],
            outfile               => $outfile_path_prefix . $outfile_alt_file_tag . $outfile_suffix,
            sequence_dict_file    => $sequence_dict_file,
        }
    );
    say {$filehandle} $NEWLINE;

    $alt_file_tag = $outfile_alt_file_tag;

    ## Remove FILTER ne PASS and on frequency
    if ( $active_parameter_href->{sv_frequency_filter} ) {

        ## Build the exclude filter command
        my $exclude_filter = _build_bcftools_filter(
            {
                annotations_ref               => \@svdb_query_annotations,
                fqf_bcftools_filter_threshold =>
                  $active_parameter_href->{fqf_bcftools_filter_threshold},
            }
        );

        say {$filehandle} q{## Remove FILTER ne PASS and frequency over threshold};
        bcftools_view(
            {
                apply_filters_ref => [qw{ PASS }],
                exclude           => $exclude_filter,
                filehandle        => $filehandle,
                infile_path       => $outfile_path_prefix . $alt_file_tag . $outfile_suffix,
                outfile_path      => $outfile_path_prefix
                  . $alt_file_tag
                  . $UNDERSCORE . q{filt}
                  . $outfile_suffix,
                output_type => q{v},
            }
        );
        say {$filehandle} $NEWLINE;

        ## Update file tag
        $alt_file_tag .= $UNDERSCORE . q{filt};
    }

    ## Remove common variants
    if ( $active_parameter_href->{sv_frequency_filter} ) {

        say {$filehandle} q{## Remove common variants};
        vcfanno(
            {
                filehandle             => $filehandle,
                infile_path            => $outfile_path_prefix . $alt_file_tag . $outfile_suffix,
                luafile_path           => $active_parameter_href->{vcfanno_functions},
                stderrfile_path_append => $stderrfile_path,
                toml_configfile_path   => $active_parameter_href->{sv_vcfanno_config},
            }
        );
        print {$filehandle} $PIPE . $SPACE;

        ## Update file tag
        $alt_file_tag .= $UNDERSCORE . q{bcftools_filter};

        my %vcfanno_config = read_from_file(
            {
                format => q{toml},
                path   => $active_parameter_href->{sv_vcfanno_config},
            }
        );

        ## Store vcfanno annotations
        my @vcfanno_annotations;

      ANNOTATION:
        foreach my $annotation_href ( @{ $vcfanno_config{annotation} } ) {

            push @vcfanno_annotations, @{ $annotation_href->{names} };
        }
        ## Build the exclude filter command
        my $exclude_filter = _build_bcftools_filter(
            {
                annotations_ref               => \@vcfanno_annotations,
                fqf_bcftools_filter_threshold =>
                  $active_parameter_href->{fqf_bcftools_filter_threshold},
            }
        );

        bcftools_filter(
            {
                exclude                => $exclude_filter,
                filehandle             => $filehandle,
                infile_path            => $DASH,
                outfile_path           => $outfile_path_prefix . $alt_file_tag . $outfile_suffix,
                output_type            => q{v},
                stderrfile_path_append => $stderrfile_path,
            }
        );
        say {$filehandle} $NEWLINE;
    }

    ## Then we have something to rename
    if ( $alt_file_tag ne $EMPTY_STR ) {

        ## Writes sbatch code to supplied filehandle to sort variants in vcf format
        sort_vcf(
            {
                active_parameter_href => $active_parameter_href,
                filehandle            => $filehandle,
                infile_paths_ref      => [ $outfile_path_prefix . $alt_file_tag . $outfile_suffix ],
                outfile               => $outfile_path,
                sequence_dict_file    => $sequence_dict_file,
            }
        );
        say {$filehandle} $NEWLINE;
    }

    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => q{sv_annotate},
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

sub _build_bcftools_filter {

## Function : Build the exclude filter command
## Returns  :
## Arguments: $fqf_bcftools_filter_threshold => Exclude variants with frequency above filter threshold
##          : $annotations_ref               => Annotations to use in filtering

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fqf_bcftools_filter_threshold;
    my $annotations_ref;

    my $tmpl = {
        fqf_bcftools_filter_threshold => {
            defined     => 1,
            required    => 1,
            store       => \$fqf_bcftools_filter_threshold,
            strict_type => 1,
        },
        annotations_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$annotations_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $exclude_filter;
    my $threshold = $SPACE . q{>} . $SPACE . $fqf_bcftools_filter_threshold . $SPACE;

    $exclude_filter =
        $DOUBLE_QUOTE
      . q{INFO/}
      . join( $threshold . $PIPE . $SPACE . q{INFO/}, @{$annotations_ref} )
      . $threshold
      . $DOUBLE_QUOTE;
    return $exclude_filter;
}

1;
