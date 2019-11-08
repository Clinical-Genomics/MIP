package MIP::Recipes::Analysis::Star_aln;

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
use MIP::Constants qw{ $ASTERISK $DOT $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.13;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_star_aln };

}

sub analysis_star_aln {

## Function : Alignment of fastq files using star aln
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_id;
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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Program::Picardtools qw{ picardtools_addorreplacereadgroups };
    use MIP::Program::Star qw{ star_aln };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Sample_info
      qw{ get_read_group get_sequence_run_type set_recipe_outfile_in_sample_info set_recipe_metafile_in_sample_info set_processing_metafile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ## PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            file_info_href => $file_info_href,
            id             => $sample_id,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my @infile_paths         = @{ $io{in}{file_paths} };
    my @infile_names         = @{ $io{in}{file_names} };
    my @infile_name_prefixes = @{ $io{in}{file_name_prefixes} };
    my $recipe_mode          = $active_parameter_href->{$recipe_name};
    my $job_id_chain         = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
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
                file_name_prefixes_ref => \@{ $infile_lane_prefix_href->{$sample_id} },
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
                temp_directory         => $temp_directory,
            }
        )
    );
    my $outfile_suffix        = $io{out}{file_suffix};
    my @outfile_name_prefixes = @{ $io{out}{file_name_prefixes} };
    my @outfile_paths         = @{ $io{out}{file_paths} };
    my @outfile_path_prefixes = @{ $io{out}{file_path_prefixes} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    # Too avoid adjusting infile_index in submitting to jobs
    my $paired_end_tracker = 0;

    ## Perform per single-end or read pair
  INFILE_PREFIX:
    while ( my ( $infile_index, $infile_prefix ) =
        each @{ $infile_lane_prefix_href->{$sample_id} } )
    {

        ## Assign file features
        my $outfile_name_prefix = $outfile_name_prefixes[$infile_index];
        my $outfile_path        = $outfile_paths[$infile_index];
        my $outfile_path_prefix = $outfile_path_prefixes[$infile_index];

        # Collect paired-end or single-end sequence run mode
        my $sequence_run_type = get_sequence_run_type(
            {
                infile_lane_prefix => $infile_prefix,
                sample_id          => $sample_id,
                sample_info_href   => $sample_info_href,
            }
        );

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
                recipe_directory                => $recipe_name,
                recipe_name                     => $recipe_name,
                process_time                    => $recipe_resource{time},
                sleep                           => 1,
                source_environment_commands_ref => $recipe_resource{load_env_ref},
                temp_directory                  => $temp_directory,
            }
        );

        ### SHELL

        ## Star aln
        say {$filehandle} q{## Aligning reads with } . $recipe_name;

        ### Get parameters
        ## Infile(s)
        my @fastq_files = $infile_paths[$paired_end_tracker];

        # If second read direction is present
        if ( $sequence_run_type eq q{paired-end} ) {

            # Increment to collect correct read 2 from %infile
            $paired_end_tracker++;
            push @fastq_files, $infile_paths[$paired_end_tracker];
        }
        my $referencefile_dir_path =
            $active_parameter_href->{star_aln_reference_genome}
          . $file_info_href->{star_aln_reference_genome}[0];

        star_aln(
            {
                filehandle          => $filehandle,
                align_intron_max    => $active_parameter_href->{align_intron_max},
                align_mates_gap_max => $active_parameter_href->{align_mates_gap_max},
                align_sjdb_overhang_min =>
                  $active_parameter_href->{align_sjdb_overhang_min},
                chim_junction_overhang_min =>
                  $active_parameter_href->{chim_junction_overhang_min},
                chim_out_type       => $active_parameter_href->{chim_out_type},
                chim_segment_min    => $active_parameter_href->{chim_segment_min},
                genome_dir_path     => $referencefile_dir_path,
                infile_paths_ref    => \@fastq_files,
                outfile_name_prefix => $outfile_path_prefix . $DOT,
                thread_number       => $recipe_resource{core_number},
                two_pass_mode       => $active_parameter_href->{two_pass_mode},
            },
        );
        say {$filehandle} $NEWLINE;

        ## Increment paired end tracker
        $paired_end_tracker++;

        my %read_group = get_read_group(
            {
                infile_prefix    => $infile_prefix,
                platform         => $active_parameter_href->{platform},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        picardtools_addorreplacereadgroups(
            {
                create_index => q{true},
                filehandle   => $filehandle,
                infile_path  => $outfile_path_prefix
                  . $DOT
                  . q{Aligned.sortedByCoord.out}
                  . $outfile_suffix,
                java_jar =>
                  catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
                java_use_large_pages    => $active_parameter_href->{java_use_large_pages},
                memory_allocation       => q{Xmx1g},
                outfile_path            => $outfile_path,
                readgroup_id            => $read_group{id},
                readgroup_library       => $read_group{lb},
                readgroup_platform      => $read_group{pl},
                readgroup_platform_unit => $read_group{pu},
                readgroup_sample        => $read_group{sm},
            },
        );

        say {$filehandle} $NEWLINE;

        ## Close filehandleS
        close $filehandle or $log->logcroak(q{Could not close filehandle});

        if ( $recipe_mode == 1 ) {

            ## Collect QC metadata info for later use
            set_recipe_outfile_in_sample_info(
                {
                    infile           => $outfile_name_prefix,
                    path             => $outfile_path,
                    recipe_name      => $recipe_name,
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );

            my $star_aln_log = $outfile_path_prefix . $DOT . q{Log.final.out};
            set_recipe_metafile_in_sample_info(
                {
                    infile           => $outfile_name_prefix,
                    path             => $star_aln_log,
                    recipe_name      => $recipe_name,
                    metafile_tag     => q{log},
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );

            my $most_complete_format_key =
              q{most_complete} . $UNDERSCORE . substr $outfile_suffix, 1;
            set_processing_metafile_in_sample_info(
                {
                    metafile_tag     => $most_complete_format_key,
                    path             => $outfile_path,
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );

            submit_recipe(
                {
                    base_command            => $profile_base_command,
                    case_id                 => $case_id,
                    dependency_method       => q{sample_to_sample_parallel},
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_chain            => $job_id_chain,
                    job_id_href             => $job_id_href,
                    log                     => $log,
                    recipe_file_path        => $recipe_file_path,
                    recipe_files_tracker    => $infile_index,
                    sample_id               => $sample_id,
                    submission_profile => $active_parameter_href->{submission_profile},
                }
            );
        }
    }
    return 1;
}

1;
