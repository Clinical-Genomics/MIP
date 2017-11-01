package MIP::Recipes::Analysis::Bwa_rapid;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use autodie qw{:all};
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Basename qw(fileparse);
use File::Spec::Functions qw(catdir catfile);

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{Exporter};

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(analysis_bwa_rapid);

}

##Constants
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub bwa_mem_rapid {

##bwa_mem_rapid

##Function : Performs alignment of single and paired-end as well as interleaved fastq(.gz) files.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infiles_ref, $infile_lane_prefix_href, $job_id_href, $insample_directory, $outsample_directory, $sample_id, $program_name, $family_id, $outaligner_dir, $temp_directory
##         : $parameter_href          => Parameter hash {REF}
##         : $active_parameter_href   => Active parameters for this analysis hash {REF}
##         : $sample_info_href        => Info on samples and family hash {REF}
##         : $file_info_href          => File info hash {REF}
##         : $infiles_ref             => Infiles hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href             => Job id hash {REF}
##         : $insample_directory      => In sample directory
##         : $outsample_directory     => Out sample directory
##         : $sample_id               => Sample id
##         : $program_name            => Program name
##         : $family_id               => Family id
##         : $outaligner_dir          => Outaligner_dir used in the analysis
##         : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
    my $temp_directory;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infiles_ref;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $insample_directory;
    my $outsample_directory;
    my $sample_id;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href
        },
        infiles_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infiles_ref
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        insample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$insample_directory
        },
        outsample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outsample_directory
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use MIP::IO::Files qw(migrate_file);
    use MIP::Program::Alignment::Samtools qw(samtools_view samtools_stats);
    use MIP::Program::Alignment::Bwa qw(bwa_mem run_bwamem);
    use Program::Variantcalling::Bedtools qw (intersectbed);
    use MIP::Program::Alignment::Sambamba qw(sambamba_sort);
    use MIP::QC::Record qw(add_program_outfile_to_sample_info);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_step_in_parallel);

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};
    my $jobid_chain      = $parameter_href->{$mip_program_name}{chain};

    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};
    my $infile_size;
    my $total_sbatch_counter = 0;

    # Too avoid adjusting infile_index in submitting to jobs
    my $paired_end_tracker = 0;

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Assign directories
    # Used downstream
    $parameter_href->{$mip_program_name}{$sample_id}{indirectory} =
      $outsample_directory;

    ## Assign file tags
    my $outfile_tag =
      $file_info_href->{$sample_id}{$mip_program_name}{file_tag};

    ### Assign suffix
    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
            jobid_chain    => $jobid_chain,
            file_suffix => $parameter_href->{$mip_program_name}{outfile_suffix},
        }
    );

    my $uncompressed_bam_output;
    if ( $outfile_suffix eq q{.bam} ) {

        # Used downstream in samtools view
        $uncompressed_bam_output = 1;
    }

    ## Collect fastq file(s) size and interleaved info
    while ( my ( $infile_index, $infile_prefix ) =
        each( @{ $infile_lane_prefix_href->{$sample_id} } ) )
    {

        ## Assign file tags
        my $file_path_prefix = catfile( $temp_directory, $infile_prefix );
        my $outfile_path_prefix = $file_path_prefix . $outfile_tag;

        # Collect paired-end or single-end sequence run mode
        my $sequence_run_mode =
          $sample_info_href->{sample}{$sample_id}{file}{$infile_prefix}
          {sequence_run_type};

        my $interleaved_fastq_file =
          $sample_info_href->{sample}{$sample_id}{file}{$infile_prefix}
          {interleaved};
        my $fastq_file_first = $infiles_ref->[$infile_index];

        # Initiate for potential second read
        my $fastq_file_second;

        ## Fastq.gz
        if ( $fastq_file_first =~ /.fastq.gz$/ ) {

            ## Files are already gz and presently the scalar for compression has not been investigated. Therefore no automatic time allocation can be performed.

            # If second read direction is present
            if ( $sequence_run_mode eq q{paired-end} ) {

                $fastq_file_second =
                  $infiles_ref->[ $infile_index + $infile_index ];
                $infile_size =
                  -s catfile( $insample_directory, $fastq_file_second );
            }
            else {
                # Single end read

                $infile_size =
                  -s catfile( $insample_directory, $fastq_file_first );
            }
        }
        else {
            # Files are in fastq format and not compressed

            # If second read direction is present
            if ( $sequence_run_mode eq q{paired-end} ) {
                $fastq_file_second =
                  $infiles_ref->[ $infile_index + $infile_index ];

                ## Collect .fastq file size to enable estimation of time required for aligning, +1 for syncing multiple infiles per sample_id. Hence, filesize will be calculated on read2 (should not matter).
                $infile_size =
                  -s catfile( $insample_directory, $fastq_file_second );
            }
            else {
                #Single end read

                $infile_size =
                  -s catfile( $insample_directory, $fastq_file_first );
            }
        }

        ## Parallelize alignment by splitting of alignment processes as the files are read
        if ( $consensus_analysis_type eq q{rapid} ) {

            my $seq_length = $sample_info_href->{sample}{$sample_id}{file}
              {$infile_prefix}{sequence_length};
            my ( $number_nodes, $read_nr_of_lines ) =
              determine_nr_of_rapid_nodes(
                {
                    seq_length  => $seq_length,
                    infile_size => $infile_size,
                }
              );

            # Parallization for each file handled
            for (
                my $sbatch_counter = 0 ;
                $sbatch_counter < $number_nodes - 1 ;
                $sbatch_counter++
              )
            {

                ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
                my ( $file_name, $program_info_path ) = setup_script(
                    {
                        active_parameter_href => $active_parameter_href,
                        job_id_href           => $job_id_href,
                        FILEHANDLE            => $FILEHANDLE,
                        directory_id          => $sample_id,
                        program_name          => $program_name,
                        program_directory     => lc $outaligner_dir,
                        core_number =>
                          $active_parameter_href->{max_cores_per_node},
                        process_time => $time,
                        sleep        => 1,
                    }
                );

                # Split to enable submission to &sample_info_qc later
                my ( $volume, $directory, $stderr_file ) =
                  File::Spec->splitpath( $program_info_path . q{.stderr.txt} );

                # Constant for gz files
                my $read_start = $sbatch_counter * $read_nr_of_lines;

                # Constant for gz files
                my $read_stop = $read_start + ceil( $read_nr_of_lines + 1 );

                my $infile;

                # If second read direction is present
                if ( $sequence_run_mode eq 'paired-end' ) {

                    # For required .fastq file
                    $infile = $fastq_file_second;
                }
                else {
                    # Single end read

                    # For required .fastq file
                    $infile = $fastq_file_first;
                }

                ## BWA Mem for each read batch
                print $FILEHANDLE q{bwa mem };

              # Mark shorter split hits as secondary (for Picard compatibility).
                print $FILEHANDLE q{-M };

                # Number of threads
                print $FILEHANDLE q{-t }
                  . $active_parameter_href->{max_cores_per_node} . " ";

                if ($interleaved_fastq_file) {

                    print $FILEHANDLE "-p ";    #interleaved fastq mode
                }

                ## Read group header line
                print $FILEHANDLE q?-R "@RG\t?;
                print $FILEHANDLE q?ID:? . $infile_prefix . q?\t?;
                print $FILEHANDLE q?SM:? . $sample_id . q?\t?;
                print $FILEHANDLE q?PL:?
                  . $active_parameter_href->{platform} . q?" ?;

                print $FILEHANDLE $active_parameter_href
                  ->{human_genome_reference} . " ";    #Reference
                print $FILEHANDLE "<( ";               #Pipe to BWA Mem (Read 1)
                print $FILEHANDLE "zcat ";             #Decompress Read 1
                print $FILEHANDLE catfile( $insample_directory, $infile )
                  . " ";                               #Read 1
                print $FILEHANDLE "| ";                #Pipe
                print $FILEHANDLE q?perl -ne 'if ( ($.>?
                  . $read_start
                  . q?) && ($.<?
                  . $read_stop
                  . q?) ) {print $_;}' ?;    #Limit to sbatch script interval
                print $FILEHANDLE ") ";      #End Read 1

                if ( $sequence_run_mode eq 'paired-end' )
                {                            #Second read direction if present

                    print $FILEHANDLE "<( ";      #Pipe to BWA Mem (Read 2)
                    print $FILEHANDLE "zcat ";    #Decompress Read 2
                    print $FILEHANDLE catfile( $insample_directory,
                        $infiles_ref->[ $infile_index + $infile_index + 1 ] )
                      . " ";                      #Read 2
                    print $FILEHANDLE "| ";       #Pipe
                    print $FILEHANDLE q?perl -ne 'if ( ($.>?
                      . $read_start
                      . q?) && ($.<?
                      . $read_stop
                      . q?) ) {print $_;}' ?;   #Limit to sbatch script interval
                    print $FILEHANDLE ") ";     #End Read 2
                }

                print $FILEHANDLE
                  "| ";    #Pipe SAM to BAM conversion of aligned reads
                samtools_view(
                    {
                        infile_path => "-",
                        FILEHANDLE  => $FILEHANDLE,
                        thread_number =>
                          $active_parameter_href->{module_core_number}
                          {$mip_program_name},
                        auto_detect_input_format => 1,
                        with_header              => 1,
                        uncompressed_bam_output  => $uncompressed_bam_output,
                    }
                );
                print $FILEHANDLE "| ";    #Pipe
                intersectbed(
                    {
                        with_header => 1,
                        infile_path => "stdin",
                        intersectfile_path =>
                          $active_parameter_href->{bwa_mem_rapid_db},
                        outfile_path => catfile(
                            $outsample_directory,
                            $infile_prefix . "_"
                              . $sbatch_counter
                              . $outfile_suffix
                        ),
                        FILEHANDLE => $FILEHANDLE,
                    }
                );
                say $FILEHANDLE "\n";

                print $FILEHANDLE "samtools sort ";
                print $FILEHANDLE catfile( $outsample_directory,
                    $infile_prefix . "_" . $sbatch_counter . $outfile_suffix )
                  . " ";    #Infile
                say $FILEHANDLE catfile( $outsample_directory,
                    $infile_prefix . "_" . $sbatch_counter . $outfile_tag ),
                  "\n";     #OutFile

                print $FILEHANDLE "samtools index ";
                say $FILEHANDLE catfile(
                    $outsample_directory,
                    $infile_prefix . "_"
                      . $sbatch_counter
                      . $outfile_tag
                      . $outfile_suffix
                  ),
                  "\n";     #OutFile

                close($FILEHANDLE);

                if ( $mip_program_mode == 1 ) {

                    slurm_submit_job_sample_id_dependency_step_in_parallel(
                        {
                            job_id_href             => $job_id_href,
                            infile_lane_prefix_href => $infile_lane_prefix_href,
                            family_id               => $family_id,
                            sample_id               => $sample_id,
                            path                    => $jobid_chain,
                            log                     => $log,
                            sbatch_file_name        => $file_name,
                            sbatch_script_tracker   => $total_sbatch_counter
                        }
                    );
                }
                $total_sbatch_counter++;

                ## Save sbatch Counter to track how many read batch processes we have engaged
                $sample_info_href->{sample}{$sample_id}{$infile_prefix}
                  {pbwa_mem}{read_batch_process} =
                  $sbatch_counter + 1;    #Used to be  $sbatch_counter
                $sample_info_href->{sample}{$sample_id}{pbwa_mem}
                  {sbatch_batch_processes} = $total_sbatch_counter;
            }
        }
    }
    return;
}

1;
