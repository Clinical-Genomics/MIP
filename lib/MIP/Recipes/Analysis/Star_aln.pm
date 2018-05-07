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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_star_aln };

}

## Constants
Readonly my $ASTERIX    => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_star_aln {

## Function : Alignment of fastq files using star aln
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infiles_ref             => Infiles hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $insample_directory      => In sample directory
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outsample_directory     => Out sample directory
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infiles_ref;
    my $infile_lane_prefix_href;
    my $insample_directory;
    my $job_id_href;
    my $outsample_directory;
    my $parameter_href;
    my $program_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
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
        infiles_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infiles_ref,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        insample_directory => {
            defined     => 1,
            required    => 1,
            store       => \$insample_directory,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir,
            strict_type => 1,
        },
        outsample_directory => {
            defined     => 1,
            required    => 1,
            store       => \$outsample_directory,
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

    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Program::Alignment::Picardtools
      qw{ picardtools_addorreplacereadgroups };
    use MIP::Program::Alignment::Star qw{ star_aln };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_step_in_parallel };
    use MIP::QC::Record
      qw{ add_program_outfile_to_sample_info add_processing_metafile_to_sample_info };
    use MIP::Set::File qw{ set_file_suffix };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
      );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Assign file_tags
    my $outfile_tag =
      $file_info_href->{$sample_id}{$mip_program_name}{file_tag};

    ### Assign suffix
    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            file_suffix => $parameter_href->{$mip_program_name}{outfile_suffix},
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

        ## Assign file tags
        my $file_path_prefix = catfile( $temp_directory, $infile_prefix );
        my $outfile_path_prefix = $file_path_prefix . $outfile_tag;

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
                program_directory               => lc $outaligner_dir,
                program_name                    => $program_name,
                process_time                    => $time,
                source_environment_commands_ref => \@source_environment_cmds,
                temp_directory                  => $temp_directory,
            }
        );

        ## Copies file to temporary directory.
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};

        # Read 1
        my $insample_dir_fastqc_path_read_one =
          catfile( $insample_directory, $infiles_ref->[$paired_end_tracker] );
        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $insample_dir_fastqc_path_read_one,
                outfile_path => $temp_directory,
            }
        );

        # If second read direction is present
        if ( $sequence_run_mode eq q{paired-end} ) {

            my $insample_dir_fastqc_path_read_two =
              catfile( $insample_directory,
                $infiles_ref->[ $paired_end_tracker + 1 ] );

            # Read 2
            migrate_file(
                {
                    FILEHANDLE   => $FILEHANDLE,
                    infile_path  => $insample_dir_fastqc_path_read_two,
                    outfile_path => $temp_directory,
                }
            );
        }
        say {$FILEHANDLE} q{wait}, $NEWLINE;

        ## Star aln
        say {$FILEHANDLE} q{## Aligning reads with } . $program_name;

        ### Get parameters

        ## Infile(s)
        my @fastq_files =
          ( catfile( $temp_directory, $infiles_ref->[$paired_end_tracker] ) );

        # If second read direction is present
        if ( $sequence_run_mode eq q{paired-end} ) {

            # Increment to collect correct read 2 from %infile
            $paired_end_tracker = $paired_end_tracker + 1;
            push @fastq_files,
              catfile( $temp_directory, $infiles_ref->[$paired_end_tracker] );
            catfile( $temp_directory, $infiles_ref->[$paired_end_tracker] );
        }
        my $referencefile_dir_path = $active_parameter_href->{reference_dir};

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
                outfile_name_prefix => $outfile_path_prefix . $DOT,
                thread_number       => $core_number,
                two_pass_mode       => $active_parameter_href->{two_pass_mode},
            },
        );
        say {$FILEHANDLE} $NEWLINE;

        picardtools_addorreplacereadgroups(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $outfile_path_prefix
                  . $DOT
                  . q{Aligned.sortedByCoord.out}
                  . $outfile_suffix,
                java_jar => catfile(
                    $active_parameter_href->{picardtools_path},
                    q{picard.jar}
                ),
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                memory_allocation  => q{Xmx1g},
                outfile_path       => $outfile_path_prefix . $outfile_suffix,
                readgroup_id       => $infile_prefix,
                readgroup_library  => q{RNA},
                readgroup_platform => $active_parameter_href->{platform},
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
                infile_path  => $outfile_path_prefix . $ASTERIX,
                outfile_path => $outsample_directory,
            }
        );
        say {$FILEHANDLE} q{wait}, $NEWLINE;

        ## Close FILEHANDLES
        close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

        if ( $mip_program_mode == 1 ) {

            my $program_outfile_path =
              $outfile_path_prefix . $UNDERSCORE . q{bam};

            ## Collect QC metadata info for later use
            add_program_outfile_to_sample_info(
                {
                    path             => $program_outfile_path,
                    program_name     => $program_name,
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );

            my $most_complete_format_key =
              q{most_complete} . $UNDERSCORE . substr $outfile_suffix, 1;
            my $qc_metafile_path =
              catfile( $outsample_directory, $infile_prefix . $outfile_suffix );
            add_processing_metafile_to_sample_info(
                {
                    metafile_tag     => $most_complete_format_key,
                    path             => $qc_metafile_path,
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

