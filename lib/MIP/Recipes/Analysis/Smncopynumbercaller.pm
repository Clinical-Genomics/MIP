package MIP::Recipes::Analysis::Smncopynumbercaller;

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
use MIP::Constants qw{ $GENOME_VERSION $LOG_NAME $NEWLINE $SINGLE_QUOTE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_smncopynumbercaller };

}

sub analysis_smncopynumbercaller {

## Function : Call the copy number of full-length SMN1, full-length SMN2,
##          : as well as SMN2Δ7–8 (SMN2 with a deletion of Exon7-8).
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
    use MIP::Language::Perl qw{ perl_base };
    use MIP::Program::Gnu::Coreutils qw{ gnu_echo };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Smncopynumbercaller qw{ smn_caller };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{ set_file_path_to_store set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    my @infile_paths;
    my %sample_file_prefix;
  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Get the io infiles per chain and id
        my %sample_io = get_io_files(
            {
                id             => $sample_id,
                file_info_href => $file_info_href,
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
                stream         => q{in},
            }
        );
        my $file_name_prefix = $sample_io{in}{file_name_prefix};
        my $file_path_prefix = $sample_io{in}{file_path_prefix};
        my $file_suffix      = $sample_io{in}{file_suffix};
        push @infile_paths, $file_path_prefix . $file_suffix;
        $sample_file_prefix{$sample_id} = $file_name_prefix;
    }

    my %io = parse_io_outfiles(
        {
            chain_id               => $recipe{job_id_chain},
            id                     => $case_id,
            file_info_href         => $file_info_href,
            file_name_prefixes_ref => [$case_id],
            outdata_dir            => $active_parameter_href->{outdata_dir},
            parameter_href         => $parameter_href,
            recipe_name            => $recipe_name,
        }
    );
    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_name_prefix = $io{out}{file_name_prefix};
    my $outfile_path        = $io{out}{file_path};

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

    ## Create manifest file
    my $manifest_file_path = catfile( $outdir_path_prefix, q{manifest.txt} );

    gnu_echo(
        {
            enable_interpretation => 1,
            filehandle            => $filehandle,
            no_trailing_newline   => 1,
            outfile_path          => $manifest_file_path,
            strings_ref           => [ join q{\n}, @infile_paths ],
        }
    );
    say {$filehandle} $NEWLINE;

    smn_caller(
        {
            filehandle         => $filehandle,
            manifest_file_path => $manifest_file_path,
            genome_version     => $GENOME_VERSION,
            outfile_prefix     => $outfile_name_prefix,
            outdir_path        => $outdir_path_prefix,
            thread_number      => $recipe{core_number},
        }
    );
    say {$filehandle} $NEWLINE;

    _use_sample_id_in_output(
        {
            filehandle              => $filehandle,
            outfile_path            => $outfile_path,
            sample_file_prefix_href => \%sample_file_prefix,
        }
    );

    ## Close filehandle
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

        set_file_path_to_store(
            {
                format           => q{meta},
                id               => $case_id,
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

sub _use_sample_id_in_output {

## Function : Rename file_name_prefix to sample_id for sample column
## Returns  :
## Arguments: $filehandle              => Filehandle to write to
##          : $outfile_path            => Outfile path to use for search and replace
##          : $sample_file_prefix_href => Map of file_name_prefix and sample_id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $outfile_path;
    my $sample_file_prefix_href;

    my $tmpl = {
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        sample_file_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_file_prefix_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    while ( my ( $sample_id, $file_name_prefix ) = each %{$sample_file_prefix_href} ) {

        my @perl_commands = perl_base(
            {
                command_line  => 1,
                inplace       => 1,
                print         => 1,
                use_container => 1,
            }
        );
        push @perl_commands,
          ( $SINGLE_QUOTE, qq{s/$file_name_prefix/$sample_id/g}, $SINGLE_QUOTE, $outfile_path );
        say {$filehandle} join $SPACE, @perl_commands;
    }
    return;
}

1;
