package MIP::Recipes::Download::Gnomad;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catfile devnull };
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
  qw{ $DASH $DOT $EMPTY_STR $ESCAPE $FORWARD_SLASH $NEWLINE $PIPE $SINGLE_QUOTE $SPACE $UNDERSCORE };

## Constants
Readonly my $DOWNLOAD_TRIES   => 12;
Readonly my $READ_TIMEOUT_SEC => 20;
Readonly my $TIMEOUT_SEC      => 20;
Readonly my $WAIT_RETRY_SEC   => 300;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ download_gnomad download download_gnomad_mt download_gnomad_xargs};

}

sub download_gnomad {

## Function : Download gnomad
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this download hash {REF}
##          : $genome_version        => Human genome version
##          : $job_id_href           => The job_id hash {REF}
##          : $profile_base_command  => Submission profile base command
##          : $recipe_name           => Recipe name
##          : $reference_href        => Reference hash {REF}
##          : $reference_version     => Reference version
##          : $quiet                 => Quiet (no output)
##          : $temp_directory        => Temporary directory for recipe
##          : $verbose               => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $genome_version;
    my $job_id_href;
    my $recipe_name;
    my $reference_href;
    my $reference_version;

    ## Default(s)
    my $profile_base_command;
    my $quiet;
    my $temp_directory;
    my $verbose;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        genome_version => {
            store       => \$genome_version,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
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
        reference_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$reference_href,
            strict_type => 1,
        },
        reference_version => {
            defined     => 1,
            required    => 1,
            store       => \$reference_version,
            strict_type => 1,
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$quiet,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Slurm_processes qw{ slurm_submit_job_no_dependency_dead_end };
    use MIP::Program::Bcftools qw{ bcftools_index };
    use MIP::Recipes::Download::Get_reference qw{ get_reference };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger( uc q{mip_download} );

    ## Unpack parameters
    my $reference_dir = $active_parameter_href->{reference_dir};

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Filehandle(s)
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $recipe{core_number},
            directory_id                    => q{mip_download},
            filehandle                      => $filehandle,
            info_file_id                    => $genome_version . $UNDERSCORE . $reference_version,
            job_id_href                     => $job_id_href,
            memory_allocation               => $recipe{memory},
            outdata_dir                     => $reference_dir,
            outscript_dir                   => $reference_dir,
            process_time                    => $recipe{time},
            recipe_data_directory_path      => $active_parameter_href->{reference_dir},
            recipe_directory                => $recipe_name . $UNDERSCORE . $reference_version,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe{load_env_ref},
        }
    );

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    get_reference(
        {
            filehandle     => $filehandle,
            recipe_name    => $recipe_name,
            reference_dir  => $reference_dir,
            reference_href => $reference_href,
            quiet          => $quiet,
            verbose        => $verbose,
        }
    );

    my $reformated_outfile = join $UNDERSCORE,
      ( $genome_version, $recipe_name, q{reformated}, $DASH . $reference_version . q{-.vcf.gz} );
    my $reformated_outfile_path = catfile( $reference_dir, $reformated_outfile );

    my %gnomad_post_processing = (
        q{r2.0.1} => {
            arg_href => {
                info_keys_ref => [qw{ INFO/AF INFO/AF_POPMAX }],
            },
            method => \&_annotate,
        },
        q{r2.1.1} => {
            arg_href => {
                info_keys_ref => [qw{ INFO/AF INFO/AF_popmax }],
            },
            method => \&_annotate,
        },
        q{r2.1.1_sv} => {
            arg_href => {
                info_keys_ref => [qw{ INFO/AC INFO/AF INFO/POPMAX_AF }],
            },
            method => \&_annotate,
        },
    );

    $gnomad_post_processing{$reference_version}{method}->(
        {
            %{ $gnomad_post_processing{$reference_version}{arg_href} },
            filehandle   => $filehandle,
            infile_path  => catfile( $reference_dir, $reference_href->{outfile} ),
            outfile_path => $reformated_outfile_path,
        }
    );

    bcftools_index(
        {
            filehandle  => $filehandle,
            infile_path => $reformated_outfile_path,
            output_type => q{tbi},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Create AF file for bcftools roh
    _build_af_file(
        {
            filehandle        => $filehandle,
            file_name         => $reformated_outfile,
            infile_path       => $reformated_outfile_path,
            reference_dir     => $reference_dir,
            reference_version => $reference_version,
        }
    );

    ## Close filehandle
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

        ## No upstream or downstream dependencies
        slurm_submit_job_no_dependency_dead_end(
            {
                base_command     => $profile_base_command,
                job_id_href      => $job_id_href,
                log              => $log,
                sbatch_file_name => $recipe_file_path,
            }
        );
    }
    return 1;
}

sub download_gnomad_mt {

## Function : Download gnomad mitochondrial
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this download hash {REF}
##          : $genome_version        => Human genome version
##          : $job_id_href           => The job_id hash {REF}
##          : $profile_base_command  => Submission profile base command
##          : $quiet                 => Quiet (no output)
##          : $recipe_name           => Recipe name
##          : $reference_href        => Reference hash {REF}
##          : $reference_version     => Reference version
##          : $temp_directory        => Temporary directory for recipe
##          : $verbose               => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $genome_version;
    my $job_id_href;
    my $recipe_name;
    my $reference_href;
    my $reference_version;

    ## Default(s)
    my $profile_base_command;
    my $quiet;
    my $temp_directory;
    my $verbose;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        genome_version => {
            store       => \$genome_version,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        profile_base_command => {
            default     => q{sbatch},
            store       => \$profile_base_command,
            strict_type => 1,
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$quiet,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        reference_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$reference_href,
            strict_type => 1,
        },
        reference_version => {
            defined     => 1,
            required    => 1,
            store       => \$reference_version,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Slurm_processes qw{ slurm_submit_job_no_dependency_dead_end };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Recipes::Download::Get_reference qw{ get_reference };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger( uc q{mip_download} );

    ## Unpack parameters
    my $reference_dir = $active_parameter_href->{reference_dir};

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Filehandle(s)
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $recipe{core_number},
            directory_id                    => q{mip_download},
            filehandle                      => $filehandle,
            info_file_id                    => $genome_version . $UNDERSCORE . $reference_version,
            job_id_href                     => $job_id_href,
            memory_allocation               => $recipe{memory},
            outdata_dir                     => $reference_dir,
            outscript_dir                   => $reference_dir,
            process_time                    => $recipe{time},
            recipe_data_directory_path      => $active_parameter_href->{reference_dir},
            recipe_directory                => $recipe_name . $UNDERSCORE . $reference_version,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe{load_env_ref},
        }
    );

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    get_reference(
        {
            filehandle     => $filehandle,
            recipe_name    => $recipe_name,
            reference_dir  => $reference_dir,
            reference_href => $reference_href,
            quiet          => $quiet,
            verbose        => $verbose,
        }
    );

    if ( $genome_version eq q{grch37} ) {

        _reformat_for_grch37(
            {
                filehandle  => $filehandle,
                infile_path => catfile( $reference_dir, $reference_href->{outfile} )
            }
        );
    }

    ## Close filehandleS
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

        ## No upstream or downstream dependencies
        slurm_submit_job_no_dependency_dead_end(
            {
                base_command     => $profile_base_command,
                job_id_href      => $job_id_href,
                log              => $log,
                sbatch_file_name => $recipe_file_path,
            }
        );
    }
    return 1;
}

sub download_gnomad_xargs {

## Function : Download gnomad using xargs
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this download hash {REF}
##          : $genome_version        => Human genome version
##          : $job_id_href           => The job_id hash {REF}
##          : $profile_base_command  => Submission profile base command
##          : $recipe_name           => Recipe name
##          : $reference_href        => Reference hash {REF}
##          : $reference_version     => Reference version
##          : $quiet                 => Quiet (no output)
##          : $temp_directory        => Temporary directory for recipe
##          : $verbose               => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $genome_version;
    my $job_id_href;
    my $recipe_name;
    my $reference_href;
    my $reference_version;

    ## Default(s)
    my $profile_base_command;
    my $quiet;
    my $temp_directory;
    my $verbose;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        genome_version => {
            store       => \$genome_version,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
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
        reference_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$reference_href,
            strict_type => 1,
        },
        reference_version => {
            defined     => 1,
            required    => 1,
            store       => \$reference_version,
            strict_type => 1,
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$quiet,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Slurm_processes qw{ slurm_submit_job_no_dependency_dead_end };
    use MIP::Program::Bcftools qw{ bcftools_concat };
    use MIP::Program::Gnu::Coreutils qw{ gnu_rm };
    use MIP::Program::Htslib qw{ htslib_tabix };
    use MIP::Program::Wget qw{ wget };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger( uc q{mip_download} );

    ## Unpack parameters
    my $reference_dir = $active_parameter_href->{reference_dir};

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Filehandle(s)
    # Create anonymous filehandle
    my $filehandle      = IO::Handle->new();
    my $xargsfilehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $recipe{core_number},
            directory_id                    => q{mip_download},
            filehandle                      => $filehandle,
            info_file_id                    => $genome_version . $UNDERSCORE . $reference_version,
            job_id_href                     => $job_id_href,
            memory_allocation               => $recipe{memory},
            outdata_dir                     => $reference_dir,
            outscript_dir                   => $reference_dir,
            process_time                    => $recipe{time},
            recipe_data_directory_path      => $active_parameter_href->{reference_dir},
            recipe_directory                => $recipe_name . $UNDERSCORE . $reference_version,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe{load_env_ref},
        }
    );

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

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

    ## Expand to list
    my @gnomad_files = glob $reference_href->{file};
    my @gnomad_file_paths;

  GNOMAD_FILE:
    foreach my $gnomad_file (@gnomad_files) {

        my $gnomad_url = $reference_href->{url_prefix} . $gnomad_file;

        wget(
            {
                filehandle        => $xargsfilehandle,
                outfile_path      => catfile( dirname( devnull() ), q{stdout} ),
                quiet             => 1,
                read_timeout      => $READ_TIMEOUT_SEC,
                retry_connrefused => 1,
                timeout           => $TIMEOUT_SEC,
                tries             => $DOWNLOAD_TRIES,
                url               => $gnomad_url,
                verbose           => 0,
                wait_retry        => $WAIT_RETRY_SEC,
            }
        );
        print {$xargsfilehandle} $PIPE . $SPACE;

        my $gnomad_file_path = catfile( $reference_dir, $gnomad_file );
        _annotate(
            {
                escape        => $ESCAPE,
                filehandle    => $xargsfilehandle,
                info_keys_ref => [qw{ INFO/AF INFO/AF_popmax }],
                infile_path   => $DASH,
                outfile_path  => $gnomad_file_path,
            }
        );
        push @gnomad_file_paths, $gnomad_file_path;
    }

    my $outfile_path = catfile( $reference_dir, $reference_href->{outfile} );
    bcftools_concat(
        {
            filehandle       => $filehandle,
            rm_dups          => 0,
            output_type      => q{z},
            threads          => $recipe{core_number},
            outfile_path     => $outfile_path,
            infile_paths_ref => \@gnomad_file_paths,
        }
    );
    say {$filehandle} $NEWLINE;

    htslib_tabix(
        {
            filehandle  => $filehandle,
            preset      => q{vcf},
            infile_path => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Create AF file for bcftools roh
    _build_af_file(
        {
            filehandle        => $filehandle,
            file_name         => $reference_href->{outfile},
            infile_path       => $outfile_path,
            reference_dir     => $reference_dir,
            reference_version => $reference_version,
        }
    );

    gnu_rm(
        {
            filehandle  => $filehandle,
            infile_path => catfile( $reference_dir, $reference_href->{file} ),
        }
    );

    ## Close filehandles
    close $filehandle      or $log->logcroak(q{Could not close filehandle});
    close $xargsfilehandle or $log->logcroak(q{Could not close xargsfilehandle});

    if ( $recipe{mode} == 1 ) {

        ## No upstream or downstream dependencies
        slurm_submit_job_no_dependency_dead_end(
            {
                base_command     => $profile_base_command,
                job_id_href      => $job_id_href,
                log              => $log,
                sbatch_file_name => $recipe_file_path,
            }
        );
    }
    return 1;
}

sub _build_af_file {

## Function : Build allele frequency file for bcftools roh
## Returns  :
## Arguments: $file_name         => Name of downloaded file
##          : $filehandle        => Filehandle
##          : $infile_path       => Path to reformatted file
##          : $reference_dir     => Reference directory
##          : $reference_version => Reference version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_name;
    my $filehandle;
    my $infile_path;
    my $reference_dir;
    my $reference_version;

    my $tmpl = {
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        reference_dir => {
            defined     => 1,
            required    => 1,
            store       => \$reference_dir,
            strict_type => 1,
        },
        reference_version => {
            allow       => [qw{ r2.0.1 r2.1.1 r2.1.1_sv r3.1.1 r3.1.2 }],
            required    => 1,
            store       => \$reference_version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw{ remove_file_path_suffix };
    use MIP::Program::Bcftools qw{ bcftools_query };
    use MIP::Program::Htslib qw{ htslib_bgzip htslib_tabix };

    ## Don't build file for SV:s
    return if ( $reference_version eq q{r2.1.1_sv} );

    my $outfile_no_suffix = remove_file_path_suffix(
        {
            file_path         => $file_name,
            file_suffixes_ref => [qw{ .vcf .vcf.gz }],
        }
    );
    my $allele_frq_file_path = catfile( $reference_dir, $outfile_no_suffix . q{.tab.gz} );

    bcftools_query(
        {
            filehandle       => $filehandle,
            format           => q{'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n'},
            infile_paths_ref => [$infile_path],
        }
    );
    print {$filehandle} $SPACE . $PIPE . $SPACE;

    htslib_bgzip(
        {
            filehandle      => $filehandle,
            stdoutfile_path => $allele_frq_file_path,
            write_to_stdout => 1,
        }
    );
    say {$filehandle} $NEWLINE;

    htslib_tabix(
        {
            begin       => 2,
            end         => 2,
            filehandle  => $filehandle,
            infile_path => $allele_frq_file_path,
            sequence    => 1,
        }
    );
    say {$filehandle} $NEWLINE;

    return;
}

sub _annotate {

## Function : Annotate gnomad vcf
## Returns  :
## Arguments: $filehandle       => Filehandle
##          : $infile_path      => Path to infile
##          : $info_keys_ref    => INFO keys
##          : $outfile_path     => Path to reformatted file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $info_keys_ref;
    my $outfile_path;

    ## Default
    my $escape;

    my $tmpl = {
        escape => {
            default => $EMPTY_STR,
            store   => \$escape,
        },
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        info_keys_ref => {
            default     => [],
            required    => 1,
            store       => \$info_keys_ref,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Bcftools qw{ bcftools_annotate };

    ## Annotate
    ## Only include sites for which at least one of the info keys are above zero
    my $include_record = join $SPACE . $PIPE x 2 . $SPACE, map { $_ . q{>0} } @{$info_keys_ref};
    bcftools_annotate(
        {
            filehandle     => $filehandle,
            include        => $escape . $SINGLE_QUOTE . $include_record . $escape . $SINGLE_QUOTE,
            infile_path    => $infile_path,
            outfile_path   => $outfile_path,
            output_type    => q{z},
            remove_ids_ref => [ map { q{^} . $_ } @{$info_keys_ref} ],
        }
    );
    say {$filehandle} $NEWLINE;

    return;
}

sub _reformat_for_grch37 {

## Function : Rename chrM to MT for grch37
## Returns  :
## Arguments: $filehandle  => Filehandle
##          : $infile_path => Path to reformatted file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;

    my $tmpl = {
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Perl qw{ perl_nae_oneliners };
    use MIP::Program::Bcftools qw{ bcftools_view };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Program::Htslib qw{ htslib_bgzip htslib_tabix };

    bcftools_view(
        {
            filehandle  => $filehandle,
            infile_path => $infile_path,
        }
    );
    print {$filehandle} $PIPE . $SPACE;

    perl_nae_oneliners(
        {
            filehandle    => $filehandle,
            oneliner_name => q{synonyms_grch38_to_grch37},
        }
    );
    print {$filehandle} $PIPE . $SPACE;

    my $temp_outfile_path = $infile_path . q{_temp.vcf.gz};
    htslib_bgzip(
        {
            filehandle      => $filehandle,
            stdoutfile_path => $temp_outfile_path,
            write_to_stdout => 1,
        }
    );
    say {$filehandle} $NEWLINE;

    gnu_mv(
        {
            filehandle   => $filehandle,
            infile_path  => $temp_outfile_path,
            outfile_path => $infile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    htslib_tabix(
        {
            filehandle  => $filehandle,
            infile_path => $infile_path,
            preset      => q{vcf},
        }
    );
    say {$filehandle} $NEWLINE;

    return;
}

1;
