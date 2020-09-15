package MIP::Recipes::Analysis::Chromograph;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
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

    # Set the version for version checking
    our $VERSION = 1.08;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_chromograph analysis_chromograph_proband };

}

sub analysis_chromograph {

## Function : Visualize chromosomes using chromograph using tiddit coverage data
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

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Program::Tar qw{ tar };
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
            recipe_name    => q{tiddit_coverage},
            stream         => q{out},
        }
    );
    my $infile_name_prefix = $io{out}{file_name_prefix};
    my $infile_path        = $io{out}{file_path};

    my $job_id_chain = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $sample_id,
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
            active_parameter_href           => $active_parameter_href,
            core_number                     => $recipe_resource{core_number},
            directory_id                    => $sample_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            log                             => $log,
            memory_allocation               => $recipe_resource{memory},
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
        }
    );

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    ## Process the wig file from tiddit_coverage
    chromograph(
        {
            coverage_file_path => $infile_path,
            filehandle         => $filehandle,
            outdir_path        => catdir( $outfile_path_prefix, q{coverage} ),
            step               => $active_parameter_href->{tiddit_coverage_bin_size},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Generate chromosome ideograms
    chromograph(
        {
            filehandle     => $filehandle,
            ideo_file_path => $active_parameter_href->{chromograph_cytoband_file},
            outdir_path    => catdir( $outfile_path_prefix, q{ideogram} ),
        }
    );
    say {$filehandle} $NEWLINE;

    tar(
        {
            create       => 1,
            filehandle   => $filehandle,
            file_path    => $outfile_path,
            filter_gzip  => 1,
            in_paths_ref => [$outfile_path_prefix],
        }
    );
    say {$filehandle} $NEWLINE;

    ## Close filehandle
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        set_file_path_to_store(
            {
                format           => q{tar},
                id               => $sample_id,
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command         => $profile_base_command,
                case_id              => $case_id,
                dependency_method    => q{case_to_sample},
                job_id_chain         => $job_id_chain,
                job_id_href          => $job_id_href,
                job_reservation_name => $active_parameter_href->{job_reservation_name},
                log                  => $log,
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

sub analysis_chromograph_proband {

## Function : Visualize chromosomes using chromograph with tiddit_coverage and upd data
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

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mkdir gnu_sort };
    use MIP::Program::Tar qw{ tar };
    use MIP::Program::Chromograph qw{ chromograph };
    use MIP::Program::Upd qw{ upd_call };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Reference qw{ write_contigs_size_file };
    use MIP::Sample_info
      qw{ get_family_member_id set_file_path_to_store set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Program::Ucsc qw{ ucsc_bed_to_big_bed };

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
    my $infile_path        = $infile_path_prefix . q{.vcf.gz};

    my $job_id_chain = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Switch from case id to sample id for the outfiles since this is done per sample
    $infile_name_prefix =~ s/$case_id/$sample_id/xms;
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => [$infile_name_prefix],
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );

    my $outdir_path         = $io{out}{dir_path};
    my $outfile_path        = $io{out}{file_path};
    my $outfile_path_prefix = $io{out}{file_path_prefix};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $recipe_resource{core_number},
            directory_id                    => $sample_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            log                             => $log,
            memory_allocation               => $recipe_resource{memory},
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
        }
    );

    ## Create chromosome name and size file
    my $contigs_size_file_path =
      catfile( $outdir_path, q{contigs_size_file} . $DOT . q{tsv} );
    write_contigs_size_file(
        {
            fai_file_path => $active_parameter_href->{human_genome_reference}
              . $DOT . q{fai},
            outfile_path => $contigs_size_file_path,
        }
    );
    say {$filehandle} $NEWLINE;

    %io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => q{tiddit_coverage},
            stream         => q{out},
        }
    );
    my $tiddit_cov_infile_path = $io{out}{file_path};
    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    ## UPD
    ## Get family hash
    my %family_member_id =
      get_family_member_id( { sample_info_href => $sample_info_href } );

    gnu_mkdir(
        {
            filehandle       => $filehandle,
            indirectory_path => $outfile_path_prefix,
            parents          => 1,
        }
    );
    say {$filehandle} $NEWLINE;

    my @call_types         = qw{ sites regions };
    my $upd_outfile_suffix = $DOT . q{bed};

  CALL_TYPE:
    foreach my $call_type (@call_types) {

        my $upd_oufile_prefix_path = $outfile_path_prefix . $UNDERSCORE . $call_type;
        my $upd_outfile_path       = $upd_oufile_prefix_path . $upd_outfile_suffix;
        my $sort_outfile_path =
          $upd_oufile_prefix_path . $UNDERSCORE . q{sorted} . $upd_outfile_suffix;
        my $ucsc_outfile_prefix_path = catfile( $outfile_path_prefix, $call_type );

        upd_call(
            {
                af_tag       => q{GNOMADAF},
                call_type    => $call_type,
                father_id    => $family_member_id{father},
                filehandle   => $filehandle,
                infile_path  => $infile_path,
                mother_id    => $family_member_id{mother},
                outfile_path => $upd_outfile_path,
                proband_id   => $sample_id,
            }
        );
        say {$filehandle} $NEWLINE;

        say {$filehandle} q{## Sort bed file};
        gnu_sort(
            {
                filehandle   => $filehandle,
                keys_ref     => [ q{1,1}, q{2,2n} ],
                infile_path  => $upd_outfile_path,
                outfile_path => $sort_outfile_path,
            }
        );
        say {$filehandle} $NEWLINE;

        say {$filehandle} q{## Create bed index files};
        ucsc_bed_to_big_bed(
            {
                contigs_size_file_path => $contigs_size_file_path,
                filehandle             => $filehandle,
                infile_path            => $sort_outfile_path,
                outfile_path           => $ucsc_outfile_prefix_path . $DOT . q{bb},
            }
        );
        say {$filehandle} $NEWLINE;
    }

    ## Process the wig file from tiddit_coverage
    chromograph(
        {
            coverage_file_path => $tiddit_cov_infile_path,
            filehandle         => $filehandle,
            outdir_path        => catdir( $outfile_path_prefix, q{coverage} ),
            step               => $active_parameter_href->{tiddit_coverage_bin_size},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Process regions file from UPD
    chromograph(
        {
            filehandle            => $filehandle,
            outdir_path           => catdir( $outfile_path_prefix, q{upd_regions} ),
            upd_regions_file_path => $outfile_path_prefix
              . $UNDERSCORE
              . q{regions}
              . $upd_outfile_suffix,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Process sites file from UPD
    chromograph(
        {
            filehandle          => $filehandle,
            outdir_path         => catdir( $outfile_path_prefix, q{upd_sites} ),
            upd_sites_file_path => $outfile_path_prefix
              . $UNDERSCORE
              . q{sites}
              . $upd_outfile_suffix,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Generate chromosome ideograms
    chromograph(
        {
            filehandle     => $filehandle,
            ideo_file_path => $active_parameter_href->{chromograph_cytoband_file},
            outdir_path    => catdir( $outfile_path_prefix, q{ideogram} ),
        }
    );
    say {$filehandle} $NEWLINE;

    tar(
        {
            create       => 1,
            filehandle   => $filehandle,
            file_path    => $outfile_path,
            filter_gzip  => 1,
            in_paths_ref => [$outfile_path_prefix],
        }
    );
    say {$filehandle} $NEWLINE;

    ## Close filehandleS
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        set_file_path_to_store(
            {
                format           => q{tar},
                id               => $sample_id,
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command      => $profile_base_command,
                case_id           => $case_id,
                dependency_method => q{case_to_sample},
                job_id_chain      => $job_id_chain,
                job_id_href       => $job_id_href,
                log               => $log,
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
