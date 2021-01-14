package MIP::Recipes::Analysis::Chromograph;

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
use MIP::Constants qw{ $DOT $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      analysis_chromograph_cov
      analysis_chromograph_rhoviz
      analysis_chromograph_upd
    };

}

sub analysis_chromograph_cov {

## Function : Visualize chromosomes using chromograph with tiddit coverage data
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

    use MIP::Contigs qw{ delete_contig_elements };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Program::Chromograph qw{ chromograph };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Sample_info qw{ set_file_path_to_store set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

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
    my $infile_path        = $io{in}{file_path};

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Get all contigs excluding the mitochondria
    my @contigs = delete_contig_elements(
        {
            contigs_ref        => $file_info_href->{contigs},
            remove_contigs_ref => [qw{ MT }],
        }
    );
    my @outfile_name_prefixes = map { $infile_name_prefix . $UNDERSCORE . $_ } @contigs;
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $recipe{job_id_chain},
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => \@outfile_name_prefixes,
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );
    my @outfile_paths = @{ $io{out}{file_paths} };
    my $outdir_path   = $io{out}{dir_path};

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
    chromograph(
        {
            coverage_file_path => $infile_path,
            euploid            => 1,
            filehandle         => $filehandle,
            outdir_path        => $outdir_path,
            step               => $active_parameter_href->{tiddit_coverage_bin_size},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Close filehandle
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_paths[0],
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

      OUTFILE_PATH:
        foreach my $outfile_path (@outfile_paths) {

            set_file_path_to_store(
                {
                    format           => q{png},
                    id               => $sample_id,
                    path             => $outfile_path,
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                    tag              => q{tcov},
                }
            );
        }

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

sub analysis_chromograph_rhoviz {

## Function : Visualize chromosomes using chromograph with rhocall_viz data
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

    use MIP::File::Path qw{ remove_file_path_suffix };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Chromograph qw{ chromograph };
    use MIP::Program::Gnu::Coreutils qw{ gnu_cp };
    use MIP::Sample_info qw{ set_file_path_to_store set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

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
    my $infile_name_prefix = remove_file_path_suffix(
        {
            file_path         => $io{in}{file_names}[0],
            file_suffixes_ref => [ $io{in}{file_suffix} ],
        }
    );
    my %infile_path = _build_infile_path_hash(
        {
            file_path   => $io{in}{file_paths}[0],
            file_suffix => $io{in}{file_suffix},
        }
    );

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    my @outfile_name_prefixes = _build_outfile_name_prefixes(
        {
            contigs_ref        => $file_info_href->{contigs},
            infile_name_prefix => $infile_name_prefix,
            infile_path_href   => \%infile_path,
        }
    );
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $recipe{job_id_chain},
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => \@outfile_name_prefixes,
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );

    my $outdir_path   = $io{out}{dir_path};
    my @outfile_paths = @{ $io{out}{file_paths} };
    my %outfile_path  = (
        autozyg => [ grep { /autozyg/xms } @outfile_paths ],
        fracsnp => [ grep { /fracsnp/xms } @outfile_paths ],
    );

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $recipe{core_number},
            directory_id                    => $sample_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            memory_allocation               => $recipe{memory},
            process_time                    => $recipe{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe{load_env_ref},
        }
    );

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    my %chromograph_infile_path = (
        autozyg => catfile( $outdir_path, $infile_name_prefix . $DOT . q{autozyg.bed} ),
        fracsnp => catfile( $outdir_path, $infile_name_prefix . $DOT . q{fracsnp.wig} ),
    );

    ## Move and rename files so that the correct out files are generated
  FILE_TYPE:
    foreach my $file_type ( keys %infile_path ) {

        gnu_cp(
            {
                filehandle   => $filehandle,
                infile_path  => $infile_path{$file_type},
                outfile_path => $chromograph_infile_path{$file_type},
            }
        );
        print {$filehandle} $NEWLINE;
    }
    print {$filehandle} $NEWLINE;

    ## Process regions of autozygosity
    chromograph(
        {
            euploid           => 1,
            filehandle        => $filehandle,
            outdir_path       => $outdir_path,
            autozyg_file_path => $chromograph_infile_path{autozyg},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Process fraction of SNP
    chromograph(
        {
            euploid           => 1,
            filehandle        => $filehandle,
            outdir_path       => $outdir_path,
            fracsnp_file_path => $chromograph_infile_path{fracsnp},
        }
    );
    say {$filehandle} $NEWLINE;

    # Close filehandle
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_paths[0],
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        _set_chromograph_file_paths_to_store(
            {
                outfile_path_href => \%outfile_path,
                recipe_name       => $recipe_name,
                sample_id         => $sample_id,
                sample_info_href  => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_island},
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
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

sub analysis_chromograph_upd {

## Function : Visualize chromosomes using chromograph with upd data
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

    use MIP::Contigs qw{ delete_contig_elements };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Pedigree qw{ is_sample_proband_in_trio };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Chromograph qw{ chromograph };
    use MIP::Sample_info qw{ set_file_path_to_store set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Only run on proband in trio
    return
      if (
        not is_sample_proband_in_trio(
            {
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        )
      );

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
    my @infile_name_prefixes =
      map { s/$io{in}{file_suffix}//xmsr } @{ $io{in}{file_names} };
    my $infile_path_href = $io{in}{file_path_href};

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    my @outfile_name_prefixes;
    my @contigs = delete_contig_elements(
        {
            contigs_ref        => $file_info_href->{contigs},
            remove_contigs_ref => [q{MT}],
        }
    );
  INFILE_NAME_PREFIX:

    foreach my $infile_name_prefix (@infile_name_prefixes) {

        push @outfile_name_prefixes, map { $infile_name_prefix . $UNDERSCORE . $_ } @contigs;
    }
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $recipe{job_id_chain},
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => \@outfile_name_prefixes,
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );
    my $outdir_path   = $io{out}{dir_path};
    my @outfile_paths = @{ $io{out}{file_paths} };
    my %outfile_path  = %{ $io{out}{file_path_href} };

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

    ## Process regions file from UPD if wgs
    if ( $active_parameter_href->{analysis_type}{$sample_id} eq q{wgs} ) {

        chromograph(
            {
                euploid               => 1,
                filehandle            => $filehandle,
                outdir_path           => $outdir_path,
                upd_regions_file_path => $infile_path_href->{regions},
            }
        );
        say {$filehandle} $NEWLINE;
    }

    ## Process sites file from UPD
    chromograph(
        {
            euploid             => 1,
            filehandle          => $filehandle,
            outdir_path         => $outdir_path,
            upd_sites_file_path => $infile_path_href->{sites},
        }
    );
    say {$filehandle} $NEWLINE;

    # Close filehandles
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_paths[0],
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

      OUTFILE_PATH:
        while ( my ( $call_region, $outfile_path ) = each %outfile_path ) {

            my ($call_type) = $call_region =~ /(sites|regions)/xms;
            set_file_path_to_store(
                {
                    format           => q{png},
                    id               => $sample_id,
                    path             => $outfile_path,
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                    tag              => $call_type,
                }
            );
        }

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_island},
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
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

sub _build_infile_path_hash {

## Function : Build infile path hash
## Returns  : %infile_path
## Arguments: $file_path   => Infile path
##          : $file_suffix => File suffix

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $file_suffix;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        file_suffix => {
            defined     => 1,
            required    => 1,
            store       => \$file_suffix,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw{ remove_file_path_suffix };

    my $infile_path_prefix = remove_file_path_suffix(
        {
            file_path         => $file_path,
            file_suffixes_ref => [$file_suffix],
        }
    );
    my %infile_path = (
        autozyg => $infile_path_prefix . $DOT . q{bed},
        fracsnp => $infile_path_prefix . $DOT . q{wig},
    );

    return %infile_path;
}

sub _build_outfile_name_prefixes {

## Function : Build outfile name prefixes
## Returns  : @outfile_name_prefixes
## Arguments: $contigs_ref        => Active parameters for this analysis hash {REF}
##          : $infile_path_href   => Family id
##          : $infile_name_prefix => File_info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contigs_ref;
    my $infile_path_href;
    my $infile_name_prefix;

    my $tmpl = {
        contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$contigs_ref,
            strict_type => 1,
        },
        infile_path_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_path_href,
            strict_type => 1,
        },
        infile_name_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$infile_name_prefix,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Contigs qw{ delete_contig_elements };

    my @outfile_name_prefixes;
    my @outfile_contigs = delete_contig_elements(
        {
            contigs_ref        => $contigs_ref,
            remove_contigs_ref => [q{MT}],
        }
    );

  INFILE_TYPE:
    foreach my $infile_type ( keys %{$infile_path_href} ) {

        push @outfile_name_prefixes,
          map { $infile_name_prefix . $DOT . $infile_type . $UNDERSCORE . $_ } @outfile_contigs;
    }
    return @outfile_name_prefixes;
}

sub _set_chromograph_file_paths_to_store {

## Function : Set chromograph_rhoviz outfiles to store in sample_info
## Returns  :
## Arguments: $outfile_path_href => Path of file
##          : $recipe_name       => Recipe name that produced the file
##          : $sample_id         => Id associated with file (sample_id|case_id)
##          : $sample_info_href  => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $outfile_path_href;
    my $recipe_name;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        outfile_path_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$outfile_path_href,
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

  OUTFILE_TYPE:
    foreach my $outfile_type ( keys %{$outfile_path_href} ) {

      FILE_PATH:
        foreach my $outfile_path ( @{ $outfile_path_href->{$outfile_type} } ) {

            set_file_path_to_store(
                {
                    format           => q{png},
                    id               => $sample_id,
                    path             => $outfile_path,
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                    tag              => $outfile_type,
                }
            );
        }
    }
    return;
}

1;
