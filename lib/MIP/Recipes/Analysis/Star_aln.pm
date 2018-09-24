package MIP::Recipes::Analysis::Star_aln;

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

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_star_aln };

}

## Constants
Readonly my $ASTERISK    => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_star_aln {

## Function : Alignment of fastq files using star aln
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
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
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Program::Alignment::Picardtools
      qw{ picardtools_addorreplacereadgroups };
    use MIP::Program::Alignment::Star qw{ star_aln };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_step_in_parallel };
    use MIP::QC::Record
      qw{ add_program_outfile_to_sample_info add_program_metafile_to_sample_info add_processing_metafile_to_sample_info };
    use MIP::Set::File qw{ set_file_suffix };
    use MIP::Script::Setup_script qw{ setup_script };

    ## PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program mode
    my $program_mode = $active_parameter_href->{$program_name};

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            program_name   => $program_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my @infile_paths         = @{ $io{in}{file_paths} };
    my @infile_names         = @{ $io{in}{file_names} };
    my @infile_name_prefixes = @{ $io{in}{file_name_prefixes} };
    my @temp_infile_paths    = @{ $io{temp}{file_paths} };

    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Assign job_id_chain
    my $job_id_chain = $parameter_href->{$program_name}{chain};

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id       => $job_id_chain,
                id             => $sample_id,
                file_info_href => $file_info_href,
                file_name_prefixes_ref =>
                  \@{ $infile_lane_prefix_href->{$sample_id} },
                outdata_dir    => $active_parameter_href->{outdata_dir},
                parameter_href => $parameter_href,
                program_name   => $program_name,
                temp_directory => $temp_directory,
            }
        )
    );

    my $outdir_path                = $io{out}{dir_path};
    my $outfile_suffix             = $io{out}{file_suffix};
    my @outfile_name_prefixes      = @{ $io{out}{file_name_prefixes} };
    my @outfile_paths              = @{ $io{out}{file_paths} };
    my @outfile_path_prefixes      = @{ $io{out}{file_path_prefixes} };
    my @temp_outfile_path_prefixes = @{ $io{temp}{file_path_prefixes} };
    my @temp_outfile_paths         = @{ $io{temp}{file_paths} };

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Assign file_tags
    my $outfile_tag =
      $file_info_href->{$sample_id}{$program_name}{file_tag};

    ### Assign suffix
    ## Set file suffix for next module within jobid chain
    set_file_suffix(
        {
            file_suffix    => $parameter_href->{$program_name}{outfile_suffix},
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
        }
    );

    # Too avoid adjusting infile_index in submitting to jobs
    my $paired_end_tracker = 0;

    ## Perform per single-end or read pair
  INFILE_PREFIX:
    while ( my ( $infile_index, $infile_prefix ) =
        each @{ $infile_lane_prefix_href->{$sample_id} } )
    {

        ## Assign file features
        my $file_path           = $temp_outfile_paths[$infile_index];
        my $file_path_prefix    = $temp_outfile_path_prefixes[$infile_index];
        my $outfile_name_prefix = $outfile_name_prefixes[$infile_index];
        my $outfile_path        = $outfile_paths[$infile_index];
        my $outfile_path_prefix = $outfile_path_prefixes[$infile_index];

        # Collect paired-end or single-end sequence run mode
        my $sequence_run_mode =
          $sample_info_href->{sample}{$sample_id}{file}{$infile_prefix}
          {sequence_run_type};

        # Collect interleaved info
        my $interleaved_fastq_file =
          $sample_info_href->{sample}{$sample_id}{file}{$infile_prefix}
          {interleaved};

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        my ( $file_name, $program_info_path ) = setup_script(
            {
                active_parameter_href           => $active_parameter_href,
                core_number                     => $core_number,
                directory_id                    => $sample_id,
                FILEHANDLE                      => $FILEHANDLE,
                job_id_href                     => $job_id_href,
                log                             => $log,
                program_directory               => $program_name,
                program_name                    => $program_name,
                process_time                    => $time,
                sleep                           => 1,
                source_environment_commands_ref => \@source_environment_cmds,
                temp_directory                  => $temp_directory,
            }
        );

        ## Copies file to temporary directory.
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};

        # Read 1
        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $infile_paths[$paired_end_tracker],
                outfile_path => $temp_directory,
            }
        );

        # If second read direction is present
        if ( $sequence_run_mode eq q{paired-end} ) {

            # Read 2
            migrate_file(
                {
                    FILEHANDLE   => $FILEHANDLE,
                    infile_path  => $infile_paths[ $paired_end_tracker + 1 ],
                    outfile_path => $temp_directory,
                }
            );
        }
        say {$FILEHANDLE} q{wait}, $NEWLINE;

        ## Star aln
        say {$FILEHANDLE} q{## Aligning reads with } . $program_name;

        ### Get parameters
        ## Infile(s)
        my @fastq_files = $temp_infile_paths[$paired_end_tracker];

        # If second read direction is present
        if ( $sequence_run_mode eq q{paired-end} ) {

            # Increment to collect correct read 2 from %infile
            $paired_end_tracker++;
            push @fastq_files, $temp_infile_paths[$paired_end_tracker];
        }
        my $referencefile_dir_path =
            $active_parameter_href->{star_aln_reference_genome}
          . $file_info_href->{star_aln_reference_genome}[0];

        star_aln(
            {
                FILEHANDLE       => $FILEHANDLE,
                align_intron_max => $active_parameter_href->{align_intron_max},
                align_mates_gap_max =>
                  $active_parameter_href->{align_mates_gap_max},
                align_sjdb_overhang_min =>
                  $active_parameter_href->{align_sjdb_overhang_min},
                chim_junction_overhang_min =>
                  $active_parameter_href->{chim_junction_overhang_min},
                chim_segment_min => $active_parameter_href->{chim_segment_min},
                genome_dir_path  => $referencefile_dir_path,
                infile_paths_ref => \@fastq_files,
                outfile_name_prefix => $file_path_prefix . $DOT,
                thread_number       => $core_number,
                two_pass_mode       => $active_parameter_href->{two_pass_mode},
            },
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Increment paired end tracker
        $paired_end_tracker++;

        picardtools_addorreplacereadgroups(
            {
                create_index => q{true},
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $file_path_prefix
                  . $DOT
                  . q{Aligned.sortedByCoord.out}
                  . $outfile_suffix,
                java_jar => catfile(
                    $active_parameter_href->{picardtools_path},
                    q{picard.jar}
                ),
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                memory_allocation       => q{Xmx1g},
                outfile_path            => $file_path,
                readgroup_id            => $infile_prefix,
                readgroup_library       => q{RNA},
                readgroup_platform      => $active_parameter_href->{platform},
                readgroup_platform_unit => q{0},
                readgroup_sample        => $sample_id,
            },
        );

        say {$FILEHANDLE} $NEWLINE;

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $file_path_prefix . $ASTERISK,
                outfile_path => $outdir_path,
                recursive    => 1,
            }
        );
        say {$FILEHANDLE} q{wait}, $NEWLINE;

        ## Close FILEHANDLES
        close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

        if ( $program_mode == 1 ) {

            ## Collect QC metadata info for later use
            add_program_outfile_to_sample_info(
                {
                    infile           => $outfile_name_prefix,
                    path             => $outfile_path,
                    program_name     => $program_name,
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );

            my $star_aln_log = $outfile_path_prefix . $DOT . q{Log.final.out};
            add_program_metafile_to_sample_info(
                {
                    infile           => $outfile_name_prefix,
                    path             => => $star_aln_log,
                    program_name     => $program_name,
                    metafile_tag     => q{log},
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );

            my $most_complete_format_key =
              q{most_complete} . $UNDERSCORE . substr $outfile_suffix, 1;
            add_processing_metafile_to_sample_info(
                {
                    metafile_tag     => $most_complete_format_key,
                    path             => $outfile_name_prefix,
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );

            slurm_submit_job_sample_id_dependency_step_in_parallel(
                {
                    family_id               => $family_id,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    log                     => $log,
                    path                    => $job_id_chain,
                    sample_id               => $sample_id,
                    sbatch_file_name        => $file_name,
                    sbatch_script_tracker   => $infile_index,
                }
            );
        }
    }
    return;
}

1;

