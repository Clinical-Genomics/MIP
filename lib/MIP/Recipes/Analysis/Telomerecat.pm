package MIP::Recipes::Analysis::Telomerecat;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $FORWARD_SLASH $LOG_NAME $NEWLINE $SINGLE_QUOTE $SPACE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_telomerecat };

}

sub analysis_telomerecat {

## Function : Analyse telomere lengths using Telomerecat
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $case_id               => Family id
##          : $file_info_href        => File_info hash {REF}
##          : $job_id_href           => Job id hash {REF}
##          : $parameter_href        => Parameter hash {REF}
##          : $profile_base_command  => Submission profile base command
##          : $recipe_name           => Recipe name
##          : $sample_info_href      => Info on samples and case hash {REF}

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

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Program::Telomerecat qw{ telomerecat_bam2length };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Sample_info
      qw{ get_pedigree_sample_id_attributes set_file_path_to_store set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my $job_id_chain = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $recipe_mode = $active_parameter_href->{$recipe_name};

    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    my $use_sample_id_as_display_name =
      $active_parameter_href->{telomerecat_use_sample_id_as_display_name};

    my @infile_paths;
    my %sample_display;
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
        my $file_path_prefix = $sample_io{in}{file_path_prefix};
        my $file_name_prefix = $sample_io{in}{file_name_prefix};
        my $file_suffix      = $sample_io{in}{file_suffix};
        push @infile_paths, $file_path_prefix . $file_suffix;

        ## Collect display name
        my $sample_display_name = get_pedigree_sample_id_attributes(
            {
                attribute        => q{sample_display_name},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        $sample_display{ $file_name_prefix . $file_suffix } = $sample_display_name;
    }

    ## Set and get the io files per chain, id and stream
    my %io = parse_io_outfiles(
        {
            chain_id               => $job_id_chain,
            id                     => $case_id,
            file_info_href         => $file_info_href,
            file_name_prefixes_ref => [$case_id],
            outdata_dir            => $active_parameter_href->{outdata_dir},
            parameter_href         => $parameter_href,
            recipe_name            => $recipe_name,
        }
    );
    my $outfile_path = ${ $io{out}{file_paths} }[0];

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $recipe_resource{core_number},
            directory_id                    => $case_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            memory_allocation               => $recipe_resource{memory},
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
        }
    );

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    telomerecat_bam2length(
        {
            filehandle       => $filehandle,
            infile_paths_ref => \@infile_paths,
            outfile_path     => $outfile_path,
            processes        => $recipe_resource{core_number},
            temp_directory   => $active_parameter_href->{temp_directory},
        }
    );
    say {$filehandle} $NEWLINE;

    _rename_sample(
        {
            file_path                     => $outfile_path,
            filehandle                    => $filehandle,
            sample_display_href           => \%sample_display,
            use_sample_id_as_display_name => $use_sample_id_as_display_name,
        }
    );

    ## Close filehandle
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe_mode == 1 ) {

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
                base_command         => $profile_base_command,
                case_id              => $case_id,
                dependency_method    => q{sample_to_case},
                job_id_chain         => $job_id_chain,
                job_id_href          => $job_id_href,
                job_reservation_name => $active_parameter_href->{job_reservation_name},
                log                  => $log,
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

sub _rename_sample {

## Function : Change sample names in Telomerecat outfile
## Returns  :
## Arguments: $file_path                     => Path to Telomerecat outfile
##          : $filehandle                    => Filehandle
##          : $sample_display_href           => Sample name hash {REF}
##          : $use_sample_id_as_display_name => Use sample id as sample display name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $filehandle;
    my $sample_display_href;
    my $use_sample_id_as_display_name;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        sample_display_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_display_href,
            strict_type => 1,
        },
        use_sample_id_as_display_name => {
            allow       => [ undef, 0, 1 ],
            required    => 1,
            store       => \$use_sample_id_as_display_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Perl qw{ perl_base };
    use MIP::Unix::Write_to_file qw{ unix_write_to_file };

    return if $use_sample_id_as_display_name;

    say {$filehandle} q{## Rename file name to sample display name};
    while ( my ( $file_name, $sample_display_name ) = each %{$sample_display_href} ) {

        if ( $sample_display_name and not $use_sample_id_as_display_name ) {

            my @perl_cmds = perl_base(
                {
                    command_line => 1,
                    inplace      => 1,
                    print        => 1,
                }
            );

            my $regexp =
                $SINGLE_QUOTE . q?s/\A?
              . $file_name
              . $FORWARD_SLASH
              . $sample_display_name
              . $FORWARD_SLASH
              . $SINGLE_QUOTE;
            push @perl_cmds, $regexp;

            push @perl_cmds, $file_path;

            unix_write_to_file(
                {
                    commands_ref => \@perl_cmds,
                    filehandle   => $filehandle,
                    separator    => $SPACE,

                }
            );
            print {$filehandle} $NEWLINE;
        }
    }
    return;
}

1;
