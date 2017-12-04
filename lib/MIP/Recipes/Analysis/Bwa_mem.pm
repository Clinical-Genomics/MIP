package MIP::Recipes::Analysis::Bwa_mem;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use autodie qw{ :all };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Basename qw{ dirname fileparse };
use File::Spec::Functions qw{ catdir catfile devnull splitpath };

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_bwa_mem };

}

##Constants
Readonly my $ASTERIX      => q{*};
Readonly my $DOT          => q{.};
Readonly my $DOUBLE_QUOTE => q{"};
Readonly my $EMPTY_STR    => q{};
Readonly my $NEWLINE      => qq{\n};
Readonly my $PIPE         => q{|};
Readonly my $SPACE        => q{ };
Readonly my $UNDERSCORE   => q{_};

sub analysis_bwa_mem {

## Function : Performs alignment of single and paired-end as well as interleaved fastq(.gz) files.
## Returns  :
## Arguments: $parameter_href          => Parameter hash {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infiles_ref             => Infiles hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $insample_directory      => In sample directory
##          : $outsample_directory     => Out sample directory
##          : $sample_id               => Sample id
##          : $program_name            => Program name
##          : $family_id               => Family id
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

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

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
    my $temp_directory;

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

    use MIP::IO::Files qw{ migrate_file };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_step_in_parallel };
    use MIP::Program::Alignment::Bwa qw{ bwa_mem run_bwamem };
    use MIP::Program::Alignment::Samtools qw{ samtools_view samtools_stats };
    use MIP::Program::Alignment::Sambamba qw{ sambamba_sort };
    use MIP::QC::Record
      qw{ add_program_outfile_to_sample_info add_program_metafile_to_sample_info add_processing_metafile_to_sample_info };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};

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
            job_id_chain   => $job_id_chain,
            file_suffix => $parameter_href->{$mip_program_name}{outfile_suffix},
        }
    );

    my $uncompressed_bam_output;
    if ( $outfile_suffix eq q{.bam} ) {

        # Used downstream in samtools view
        $uncompressed_bam_output = 1;
    }

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
                active_parameter_href => $active_parameter_href,
                job_id_href           => $job_id_href,
                FILEHANDLE            => $FILEHANDLE,
                directory_id          => $sample_id,
                program_name          => $program_name,
                program_directory     => lc $outaligner_dir,
                core_number           => $core_number,
                process_time          => $time,
                temp_directory        => $temp_directory,
                sleep                 => 1,
            }
        );

        # Split to enable submission to %sample_info_qc later
        my ( $volume, $directory, $stderr_file ) =
          splitpath( $program_info_path . $DOT . q{stderr.txt} );

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

        ### BWA MEM
        say {$FILEHANDLE} q{## Aligning reads with }
          . $program_name
          . q{ and sorting via Sambamba};

        ## Detect version and source of the human_genome_reference: Source (hg19 or GRCh) and return the correct bwa_mem binary
        my $bwa_binary = _select_bwamem_binary(
            {
                human_genome_reference_source =>
                  $file_info_href->{human_genome_reference_source},
                human_genome_reference_version =>
                  $file_info_href->{human_genome_reference_version},
            }
        );

        ### Get parameters

        ## Infile(s)
        my $fastq_file_path =
          catfile( $temp_directory, $infiles_ref->[$paired_end_tracker] );
        my $second_fastq_file_path;

        # If second read direction is present
        if ( $sequence_run_mode eq q{paired-end} ) {

            # Increment to collect correct read 2 from %infile
            $paired_end_tracker = $paired_end_tracker + 1;
            $second_fastq_file_path =
              catfile( $temp_directory, $infiles_ref->[$paired_end_tracker] );
        }

        ## Read group header line
        my @read_group_headers = (
            $DOUBLE_QUOTE . q{@RG} . q{\t},
            q{ID:} . $infile_prefix . q{\t},
            q{SM:} . $sample_id . q{\t},
            q{PL:} . $active_parameter_href->{platform} . $DOUBLE_QUOTE,
        );

        ## Prepare for downstream processing
        # Can be either infile or instream
        my $sambamba_sort_infile;

        # Prior to ALTs in reference genome
        if ( $bwa_binary eq q{bwa mem} ) {

            bwa_mem(
                {
                    infile_path        => $fastq_file_path,
                    second_infile_path => $second_fastq_file_path,
                    idxbase            => $referencefile_path,
                    FILEHANDLE         => $FILEHANDLE,
                    thread_number      => $core_number,
                    read_group_header =>
                      join( $EMPTY_STR, @read_group_headers ),
                    interleaved_fastq_file  => $interleaved_fastq_file,
                    mark_split_as_secondary => 1,
                }
            );

            #Pipe SAM to BAM conversion of aligned reads
            print {$FILEHANDLE} $PIPE . $SPACE;

            samtools_view(
                {
                    infile_path              => q{-},
                    FILEHANDLE               => $FILEHANDLE,
                    thread_number            => $core_number,
                    auto_detect_input_format => 1,
                    with_header              => 1,
                    uncompressed_bam_output  => $uncompressed_bam_output,
                }
            );
            print {$FILEHANDLE} $PIPE . $SPACE;

            ## Set sambamba sort input; Pipe from samtools view
            $sambamba_sort_infile =
              catfile( dirname( devnull() ), q{stdin} );
        }

        # If post to ALTs in reference genome
        if ( $bwa_binary eq q{run-bwamem} ) {

            run_bwamem(
                {
                    infile_path          => $fastq_file_path,
                    second_infile_path   => $second_fastq_file_path,
                    idxbase              => $referencefile_path,
                    outfiles_prefix_path => $file_path_prefix,
                    FILEHANDLE           => $FILEHANDLE,
                    thread_number        => $core_number,
                    read_group_header =>
                      join( $EMPTY_STR, @read_group_headers ),
                    hla_typing => $active_parameter_href->{bwa_mem_hla},
                }
            );
            print {$FILEHANDLE} $PIPE . $SPACE;
            print {$FILEHANDLE} q{sh} . $SPACE;
            say   {$FILEHANDLE} $NEWLINE;

            ## Set sambamba sort input; Sort directly from run-bwakit
            $sambamba_sort_infile =
              $file_path_prefix . $DOT . q{aln} . $outfile_suffix;
        }
        ## Increment paired end tracker
        $paired_end_tracker++;

        ## Sort the output from bwa mem|run-bwamem
        sambamba_sort(
            {
                infile_path   => $sambamba_sort_infile,
                outfile_path  => $outfile_path_prefix . $outfile_suffix,
                FILEHANDLE    => $FILEHANDLE,
                show_progress => 1,
                memory_limit =>
                  $active_parameter_href->{bwa_sambamba_sort_memory_limit},
                temp_directory => $temp_directory,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        if ( $bwa_binary eq q{bwa mem} ) {

            ## BAMS, bwa_mem logs etc.
            migrate_file(
                {
                    infile_path => $outfile_path_prefix . $DOT . $ASTERIX,
                    ,
                    outfile_path => $outsample_directory,
                    FILEHANDLE   => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} q{wait}, $NEWLINE;
        }
        if ( $bwa_binary eq q{run-bwamem} ) {

            my @outfiles = (
                $outfile_path_prefix
                  . substr( $outfile_suffix, 0, 2 )
                  . $ASTERIX,
                $file_path_prefix . $DOT . q{log} . $ASTERIX,
                $file_path_prefix . $DOT . q{hla} . $ASTERIX,
            );
            foreach my $outfile (@outfiles) {

                migrate_file(
                    {
                        infile_path  => $outfile,
                        outfile_path => $outsample_directory,
                        FILEHANDLE   => $FILEHANDLE,
                    }
                );
            }
            say {$FILEHANDLE} q{wait}, $NEWLINE;
        }

        if ( $active_parameter_href->{bwa_mem_bamstats} ) {

            samtools_stats(
                {
                    infile_path => $outfile_path_prefix . $outfile_suffix,
                    FILEHANDLE  => $FILEHANDLE,
                    auto_detect_input_format => 1,
                }
            );
            print {$FILEHANDLE} $PIPE . $SPACE;

            ## Collect raw total sequences and reads mapped from samtools stats and calculate the percentage. Write it to stdout
            _add_percentage_mapped_reads_from_samtools(
                {
                    outfile_path => $outfile_path_prefix . $DOT . q{stats},
                    FILEHANDLE   => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            ## Copies file from temporary directory.
            say {$FILEHANDLE} q{## Copy file from temporary directory};
            migrate_file(
                {
                    infile_path  => $outfile_path_prefix . $DOT . q{stats},
                    outfile_path => $outsample_directory,
                    FILEHANDLE   => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} q{wait}, $NEWLINE;
        }

        if (   ( $active_parameter_href->{bwa_mem_cram} )
            && ( $outfile_suffix ne q{.cram} ) )
        {

            say {$FILEHANDLE} q{## Create CRAM file from SAM|BAM};
            samtools_view(
                {
                    infile_path  => $outfile_path_prefix . $outfile_suffix,
                    outfile_path => $outfile_path_prefix . $DOT . q{cram},
                    referencefile_path => $referencefile_path,
                    output_format      => q{cram},
                    FILEHANDLE         => $FILEHANDLE,
                    with_header        => 1,
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            ## Copies file from temporary directory.
            say {$FILEHANDLE} q{## Copy file from temporary directory};
            migrate_file(
                {
                    infile_path  => $outfile_path_prefix . $DOT . q{cram},
                    outfile_path => $outsample_directory,
                    FILEHANDLE   => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} q{wait}, $NEWLINE;
        }

        close $FILEHANDLE;

        if ( $mip_program_mode == 1 ) {

            my $most_complete_format_key =
              q{most_complete} . $UNDERSCORE . substr $outfile_suffix, 1;
            my $qc_metafile_path =
              catfile( $outsample_directory, $infile_prefix . $outfile_suffix );
            add_processing_metafile_to_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    sample_id        => $sample_id,
                    metafile_tag     => $most_complete_format_key,
                    path             => $qc_metafile_path,
                }
            );

            if (   ( $active_parameter_href->{bwa_mem_cram} )
                && ( $outfile_suffix ne q{.cram} ) )
            {

                # Required for analysisRunStatus check downstream
                my $qc_cram_path = catfile( $outsample_directory,
                    $infile_prefix . $outfile_tag . $DOT . q{cram} );
                add_program_metafile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        sample_id        => $sample_id,
                        program_name     => $program_name,
                        infile           => $infile_prefix,
                        metafile_tag     => q{cram},
                        path             => $qc_cram_path,
                    }
                );

            }
            if ( $active_parameter_href->{bwa_mem_bamstats} ) {

                ## Collect QC metadata info for later use
                my $qc_stats_outfile =
                  $infile_prefix . $outfile_tag . $DOT . q{stats};
                add_program_outfile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        sample_id        => $sample_id,
                        program_name     => q{bamstats},
                        infile           => $infile_prefix,
                        outdirectory     => $outsample_directory,
                        outfile          => $qc_stats_outfile,
                    }
                );
            }

            if ( $bwa_binary eq q{bwa mem} ) {

                add_program_outfile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        sample_id        => $sample_id,
                        program_name     => $program_name,
                        infile           => $infile_prefix,
                        outdirectory     => $directory,
                        outfile          => $stderr_file,
                    }
                );
            }
            if ( $bwa_binary eq q{run-bwamem} ) {

                my $qc_bwa_log = $infile_prefix . $DOT . q{log.bwamem};
                add_program_outfile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        sample_id        => $sample_id,
                        program_name     => $program_name,
                        infile           => $infile_prefix,
                        outdirectory     => $outsample_directory,
                        outfile          => $qc_bwa_log,
                    }
                );
            }

            slurm_submit_job_sample_id_dependency_step_in_parallel(
                {
                    job_id_href             => $job_id_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    family_id               => $family_id,
                    sample_id               => $sample_id,
                    path                    => $job_id_chain,
                    log                     => $log,
                    sbatch_file_name        => $file_name,
                    sbatch_script_tracker   => $infile_index
                }
            );
        }
    }
    return;
}

sub _select_bwamem_binary {

## Function : Detect version and source of the human_genome_reference: Source (hg19 or GRCh) and return the correct bwa_mem binary
## Returns  :
## Arguments: $human_genome_reference_source  => Human genome reference source
##          : $human_genome_reference_version => Human genome reference version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $human_genome_reference_source;
    my $human_genome_reference_version;

    my $tmpl = {
        human_genome_reference_source => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$human_genome_reference_source
        },
        human_genome_reference_version => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$human_genome_reference_version
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    Readonly my $GENOME_BUILD_VERSION_GRCH_PRIOR_ALTS => 37;
    Readonly my $GENOME_BUILD_VERSION_HG_PRIOR_ALTS   => 19;

    if ( $human_genome_reference_source eq q{GRCh} ) {

        if ( $human_genome_reference_version >
            $GENOME_BUILD_VERSION_GRCH_PRIOR_ALTS )
        {

            return q{run-bwamem};
        }
        else {
            # Human genome version less than GrCh37

            return q{bwa mem};
        }
    }
    else {
        # hgXX build

        if ( $human_genome_reference_version >
            $GENOME_BUILD_VERSION_HG_PRIOR_ALTS )
        {

            return q{run-bwamem};
        }
        else {
            # Human genome version less than hg19

            return q{bwa mem};
        }
    }
    return;
}

sub _add_percentage_mapped_reads_from_samtools {

## Function : Collect raw total sequences and reads mapped from samtools stats and calculate the percentage. Write it to stdout.
## Returns  :
## Arguments: $outfile_path => Outfile path
##          : $FILEHANDLE   => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $outfile_path;
    my $FILEHANDLE;

    my $tmpl = {
        outfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_path
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Add percentage mapped reads to samtools stats output
    # Execute perl
    print {$FILEHANDLE} q?perl -ne '?;

    # Initiate variables
    print {$FILEHANDLE} q?$raw; $map; ?;

    # Remove newline
    print {$FILEHANDLE} q?chomp $_; ?;

    # Always relay incoming line to stdout
    print {$FILEHANDLE} q?print $_, qq{\n}; ?;

    # Find raw total sequences
    print {$FILEHANDLE} q?if ($_=~/raw total sequences:\s+(\d+)/) { ?;

    # Assign raw total sequence
    print {$FILEHANDLE} q?$raw = $1; ?;

    # End if
    print {$FILEHANDLE} q?} ?;

    # Find reads mapped
    print {$FILEHANDLE} q?elsif ($_=~/reads mapped:\s+(\d+)/) { ?;

    # Assign reads mapped
    print {$FILEHANDLE} q?$map = $1; ?;

    # Calculate percentage
    print {$FILEHANDLE} q?my $percentage = ($map / $raw ) * 100; ?;

    # Write calculation to stdout
    print {$FILEHANDLE}
      q?print qq{percentage mapped reads:\t} . $percentage . qq{\n}?;

    # End elsif
    print {$FILEHANDLE} q?} ?;

    # End oneliner
    print {$FILEHANDLE} q?' ?;
    say   {$FILEHANDLE} q{>} . $SPACE . $outfile_path;

    return;
}

1;
