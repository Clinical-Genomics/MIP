package MIP::Recipes::Analysis::Variant_annotation;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname fileparse };
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
use MIP::Constants qw{ $DASH $DOT $LOG_NAME $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      analysis_variant_annotation
      analysis_variant_annotation_panel
    };

}

sub analysis_variant_annotation {

## Function : Annotate vcf with allelle frequencies
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

    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Program::Bcftools qw{ bcftools_concat bcftools_index bcftools_view };
    use MIP::Program::Vcfanno qw{ vcfanno };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
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
    my %infile_path        = %{ $io{in}{file_path_href} };

    my @contigs_size_ordered = @{ $file_info_href->{contigs_size_ordered} };

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
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
                iterators_ref    => \@contigs_size_ordered,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
            }
        )
    );

    my @outfile_paths       = @{ $io{out}{file_paths} };
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my %outfile_path        = %{ $io{out}{file_path_href} };

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
        }
    );

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    my $loqusdb_header_path = _build_loqusdb_headers(
        {
            filehandle          => $filehandle,
            outfile_path_prefix => $outfile_path_prefix,
            vcfanno_config_name => $active_parameter_href->{vcfanno_config},
        }
    );

    ## Create file commands for xargs
    my ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number      => $recipe{core_number},
            filehandle       => $filehandle,
            file_path        => $recipe_file_path,
            recipe_info_path => $recipe_info_path,
            xargsfilehandle  => $xargsfilehandle,
        }
    );

  CONTIG:
    foreach my $contig (@contigs_size_ordered) {

        ## Get parameters
        my $stderrfile_path = $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};

        vcfanno(
            {
                filehandle           => $xargsfilehandle,
                infile_path          => $infile_path{$contig},
                luafile_path         => $active_parameter_href->{vcfanno_functions},
                stderrfile_path      => $stderrfile_path,
                toml_configfile_path => $active_parameter_href->{vcfanno_config},
            }
        );
        print {$xargsfilehandle} $PIPE . $SPACE;

        ## Add loqusdb headers if defined
        _compress_and_add_loqusdb_headers(
            {
                filehandle          => $xargsfilehandle,
                infile_path         => $DASH,
                loqusdb_header_path => $loqusdb_header_path,
                outfile_path        => $outfile_path{$contig},
                stderrfile_path     => $stderrfile_path,
            }
        );
        say {$xargsfilehandle} $NEWLINE;

    }

    say {$filehandle} q{## Index outfiles};

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
    foreach my $contig (@contigs_size_ordered) {
        ## Get parameters
        my $stderrfile_path = $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};

        bcftools_index(
            {
                filehandle      => $xargsfilehandle,
                infile_path     => $outfile_path{$contig},
                output_type     => q{tbi},
                stderrfile_path => $stderrfile_path,
            }
        );
        say {$xargsfilehandle} $NEWLINE;
    }

    close $xargsfilehandle or $log->logcroak(q{Could not close xargsfilehandle});

    say {$filehandle} q{## Concatenate outfiles for CHAIN_RHOVIZ and CHAIN_UPD};

    my $concat_outfile_path = $outfile_path_prefix . $DOT . q{vcf.gz};
    bcftools_concat(
        {
            filehandle       => $filehandle,
            infile_paths_ref => \@outfile_paths,
            output_type      => q{z},
            outfile_path     => $concat_outfile_path,
            rm_dups          => 0,
            threads          => $recipe{core_number} - 1,
        }
    );
    say {$filehandle} $NEWLINE;

    bcftools_index(
        {
            infile_path => $concat_outfile_path,
            filehandle  => $filehandle,
            output_type => q{tbi},
        }
    );

    ## Close filehandles
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

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

sub analysis_variant_annotation_panel {

## Function : Annotate vcf with allelle frequencies
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

    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Program::Bcftools qw{ bcftools_index bcftools_view };
    use MIP::Program::Vcfanno qw{ vcfanno };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
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
    my $infile_path        = $io{in}{file_path};

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set and get the io files per chain, id and stream
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

    vcfanno(
        {
            filehandle           => $filehandle,
            infile_path          => $infile_path,
            luafile_path         => $active_parameter_href->{vcfanno_functions},
            toml_configfile_path => $active_parameter_href->{vcfanno_config},
        }
    );
    print {$filehandle} $PIPE . $SPACE;

    bcftools_view(
        {
            filehandle   => $filehandle,
            infile_path  => $DASH,
            outfile_path => $outfile_path,
            output_type  => q{z},
        }
    );
    say {$filehandle} $NEWLINE;

    _add_loqusdb_headers(
        {
            filehandle          => $filehandle,
            infile_path         => $outfile_path,
            outfile_path_prefix => $outfile_path_prefix,
            outfile_suffix      => $DOT . q{vcf.gz},
            vcfanno_config_name => $active_parameter_href->{vcfanno_config},
        }
    );

    say {$filehandle} q{## Index outfiles};

    bcftools_index(
        {
            filehandle  => $filehandle,
            infile_path => $outfile_path,
            output_type => q{tbi},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Close filehandleS
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

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

sub _add_loqusdb_headers {

## Function : Add relevant loqusDB headers for downstream processing
## Returns  :
## Arguments: $filehandle          => Filehandle to write to
##          : $infile_path         => Infile path to read from
##          : $outfile_path_prefix => Outfile path
##          : $outfile_suffix      => Outfile suffix
##          : $vcfanno_config_name => Name of vcfanno config

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path_prefix;
    my $outfile_suffix;
    my $vcfanno_config_name;

    my $tmpl = {
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path_prefix,
            strict_type => 1,
        },
        outfile_suffix => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_suffix,
            strict_type => 1,
        },
        vcfanno_config_name => {
            defined     => 1,
            required    => 1,
            store       => \$vcfanno_config_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Bcftools qw{ bcftools_annotate };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mv};

    my $loqusdb_header_path = _build_loqusdb_headers(
        {
            filehandle          => $filehandle,
            outfile_path_prefix => $outfile_path_prefix,
            vcfanno_config_name => $vcfanno_config_name,
        }
    );

    return if ( not $loqusdb_header_path );

    my $annotate_outfile_path = $outfile_path_prefix . $UNDERSCORE . q{annotated.vcf.gz};
    bcftools_annotate(
        {
            filehandle      => $filehandle,
            headerfile_path => $loqusdb_header_path,
            infile_path     => $infile_path,
            outfile_path    => $annotate_outfile_path,
            output_type     => q{z},
        }
    );
    say {$filehandle} $NEWLINE;

    gnu_mv(
        {
            filehandle   => $filehandle,
            infile_path  => $annotate_outfile_path,
            outfile_path => $outfile_path_prefix . $outfile_suffix,
        }
    );
    say {$filehandle} $NEWLINE;
    return;
}

sub _build_loqusdb_headers {

## Function : Build relevant loqusDB headers for downstream processing
## Returns  : $loqusdb_header_path
## Arguments: $filehandle          => Filehandle to write to
##          : $outfile_path_prefix => Outfile path
##          : $vcfanno_config_name => Name of vcfanno config

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $outfile_path_prefix;
    my $vcfanno_config_name;

    my $tmpl = {
        filehandle          => { store => \$filehandle, },
        outfile_path_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path_prefix,
            strict_type => 1,
        },
        vcfanno_config_name => {
            defined     => 1,
            required    => 1,
            store       => \$vcfanno_config_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Io::Read qw{ read_from_file };
    use MIP::Language::Perl qw{ perl_nae_oneliners };
    use MIP::Program::Bcftools qw{ bcftools_view };

    my $loqusdb_reference_file;

    my %vcfanno_config = read_from_file(
        {
            format => q{toml},
            path   => $vcfanno_config_name,
        }
    );

  ANNOTATION:
    foreach my $annotation_href ( @{ $vcfanno_config{annotation} } ) {

        $loqusdb_reference_file =
          $annotation_href->{file} =~ /loqusdb_\w+_\w+-/xsm ? $annotation_href->{file} : undef;

        last ANNOTATION if ($loqusdb_reference_file);
    }

    ## Nothing to process - skip
    return if ( not $loqusdb_reference_file );

    say {$filehandle} q{## Build loqusdb headers};

    bcftools_view(
        {
            filehandle  => $filehandle,
            header_only => 1,
            infile_path => $loqusdb_reference_file,
        }
    );
    print {$filehandle} $PIPE . $SPACE;

    my $loqusdb_header_path = $outfile_path_prefix . $DOT . q{loqusdb_header};
    perl_nae_oneliners(
        {
            filehandle      => $filehandle,
            oneliner_name   => q{get_vcf_loqusdb_headers},
            stdoutfile_path => $loqusdb_header_path,
            use_container   => 1,
        }
    );
    say {$filehandle} $NEWLINE;

    return $loqusdb_header_path;
}

sub _compress_and_add_loqusdb_headers {

## Function : Compress and relevant loqusDB headers for downstream processing
## Returns  :
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $loqusdb_header_path    => Path to loqusdb header file
##          : $outfile_path           => Outfile path
##          : $stderrfile_path        => Name of vcfanno config

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $loqusdb_header_path;
    my $outfile_path;
    my $stderrfile_path;

    my $tmpl = {
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        loqusdb_header_path => {
            store       => \$loqusdb_header_path,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        stderrfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$stderrfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Bcftools qw{ bcftools_annotate bcftools_view };

    if ( not $loqusdb_header_path ) {

        bcftools_view(
            {
                filehandle             => $filehandle,
                infile_path            => $infile_path,
                outfile_path           => $outfile_path,
                output_type            => q{z},
                stderrfile_path_append => $stderrfile_path,
            }
        );
    }
    else {

        bcftools_annotate(
            {
                filehandle             => $filehandle,
                headerfile_path        => $loqusdb_header_path,
                infile_path            => $infile_path,
                outfile_path           => $outfile_path,
                output_type            => q{z},
                stderrfile_path_append => $stderrfile_path,
            }
        );
    }
    return;
}

1;
