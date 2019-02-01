package MIP::Recipes::Analysis::Bwa_mem;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname fileparse };
use File::Spec::Functions qw{ catdir catfile devnull splitpath };
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
    our $VERSION = 1.11;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_bwa_mem };

}

## Constants
Readonly my $ASTERISK     => q{*};
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
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
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
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
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
    use MIP::Get::Parameter
      qw{ get_recipe_parameters get_recipe_attributes get_read_group };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Alignment::Bwa qw{ bwa_mem run_bwamem };
    use MIP::Program::Alignment::Samtools qw{ samtools_stats samtools_view };
    use MIP::Program::Alignment::Sambamba qw{ sambamba_sort };
    use MIP::QC::Sample_info
      qw{ set_processing_metafile_in_sample_info set_recipe_metafile_in_sample_info set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP recipe name
    my $recipe_mode = $active_parameter_href->{$recipe_name};

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my @infile_paths         = @{ $io{in}{file_paths} };
    my @infile_names         = @{ $io{in}{file_names} };
    my @infile_name_prefixes = @{ $io{in}{file_name_prefixes} };
    my @temp_infile_paths    = @{ $io{temp}{file_paths} };

    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my $job_id_chain            = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
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

    ## Assign file tags
    my $outfile_tag =
      $file_info_href->{$sample_id}{$recipe_name}{file_tag};

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
          $sample_info_href->{sample}{$sample_id}{file}{$infile_prefix}{interleaved};

        ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
        my ( $recipe_file_path, $recipe_info_path ) = setup_script(
            {
                active_parameter_href           => $active_parameter_href,
                core_number                     => $core_number,
                directory_id                    => $sample_id,
                FILEHANDLE                      => $FILEHANDLE,
                job_id_href                     => $job_id_href,
                log                             => $log,
                recipe_directory                => $recipe_name,
                recipe_name                     => $recipe_name,
                process_time                    => $time,
                sleep                           => 1,
                source_environment_commands_ref => \@source_environment_cmds,
                temp_directory                  => $temp_directory,
            }
        );

        # Split to enable submission to %sample_info_qc later
        my ( $volume, $directory, $stderr_file ) =
          splitpath( $recipe_info_path . $DOT . q{stderr.txt} );

        ### SHELL:

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

        ### BWA MEM
        say {$FILEHANDLE} q{## Aligning reads with }
          . $recipe_name
          . q{ and sorting via Sambamba};

        ## Detect version and source of the human_genome_reference: Source (hg19 or GRCh) and return the correct bwa_mem binary
        my $bwa_binary = _get_bwamem_binary(
            {
                human_genome_reference_source =>
                  $file_info_href->{human_genome_reference_source},
                human_genome_reference_version =>
                  $file_info_href->{human_genome_reference_version},
            }
        );

        ### Get parameters

        ## Infile(s)
        my $fastq_file_path = $temp_infile_paths[$paired_end_tracker];
        my $second_fastq_file_path;

        # If second read direction is present
        if ( $sequence_run_mode eq q{paired-end} ) {

            # Increment to collect correct read 2
            $paired_end_tracker     = $paired_end_tracker + 1;
            $second_fastq_file_path = $temp_infile_paths[$paired_end_tracker];
        }

        ## Read group header line
        my %read_group = get_read_group(
            {
                infile_prefix    => $infile_prefix,
                platform         => $active_parameter_href->{platform},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        my @read_group_headers = (
            $DOUBLE_QUOTE . q{@RG} . q{\t},
            q{ID:} . $read_group{id} . q{\t},
            q{SM:} . $read_group{sm} . q{\t},
            q{PL:} . $read_group{pl} . q{\t},
            q{PU:} . $read_group{pu} . q{\t},
            q{LB:} . $read_group{lb} . $DOUBLE_QUOTE,
        );

        ## Prepare for downstream processing
        # Can be either infile or instream
        my $sambamba_sort_infile;

        # Prior to ALTs in reference genome
        if ( $bwa_binary eq q{bwa mem} ) {

            bwa_mem(
                {
                    FILEHANDLE              => $FILEHANDLE,
                    idxbase                 => $referencefile_path,
                    infile_path             => $fastq_file_path,
                    interleaved_fastq_file  => $interleaved_fastq_file,
                    mark_split_as_secondary => 1,
                    read_group_header       => join( $EMPTY_STR, @read_group_headers ),
                    soft_clip_sup_align =>
                      $active_parameter_href->{bwa_soft_clip_sup_align},
                    second_infile_path => $second_fastq_file_path,
                    thread_number      => $core_number,
                }
            );

            # Pipe SAM to BAM conversion of aligned reads
            print {$FILEHANDLE} $PIPE . $SPACE;

            samtools_view(
                {
                    auto_detect_input_format => 1,
                    FILEHANDLE               => $FILEHANDLE,
                    infile_path              => q{-},
                    thread_number            => $core_number,
                    uncompressed_bam_output  => $uncompressed_bam_output,
                    with_header              => 1,
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
                    FILEHANDLE           => $FILEHANDLE,
                    hla_typing           => $active_parameter_href->{bwa_mem_hla},
                    infile_path          => $fastq_file_path,
                    idxbase              => $referencefile_path,
                    outfiles_prefix_path => $file_path_prefix,
                    read_group_header    => join( $EMPTY_STR, @read_group_headers ),
                    second_infile_path   => $second_fastq_file_path,
                    thread_number        => $core_number,
                }
            );
            print {$FILEHANDLE} $PIPE . $SPACE;
            print {$FILEHANDLE} q{sh} . $SPACE;
            say   {$FILEHANDLE} $NEWLINE;

            ## Set sambamba sort input; Sort directly from run-bwakit
            $sambamba_sort_infile = $file_path_prefix . $DOT . q{aln} . $outfile_suffix;
        }
        ## Increment paired end tracker
        $paired_end_tracker++;

        ## Sort the output from bwa mem|run-bwamem
        sambamba_sort(
            {
                FILEHANDLE    => $FILEHANDLE,
                infile_path   => $sambamba_sort_infile,
                memory_limit  => $active_parameter_href->{bwa_sambamba_sort_memory_limit},
                outfile_path  => $file_path,
                show_progress => 1,
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
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => $file_path_prefix . $DOT . $ASTERISK,
                    ,
                    outfile_path => $outdir_path,
                }
            );
            say {$FILEHANDLE} q{wait}, $NEWLINE;
        }
        if ( $bwa_binary eq q{run-bwamem} ) {

            my @run_bwa_files = (
                $file_path . substr( $outfile_suffix, 0, 2 ) . $ASTERISK,
                $file_path_prefix . $DOT . q{log} . $ASTERISK,
                $file_path_prefix . $DOT . q{hla} . $ASTERISK,
            );

          OUTFILE:
            foreach my $run_bwa_file (@run_bwa_files) {

                migrate_file(
                    {
                        FILEHANDLE   => $FILEHANDLE,
                        infile_path  => $run_bwa_file,
                        outfile_path => $outdir_path,
                    }
                );
            }
            say {$FILEHANDLE} q{wait}, $NEWLINE;
        }

        if ( $active_parameter_href->{bwa_mem_bamstats} ) {

            samtools_stats(
                {
                    auto_detect_input_format => 1,
                    FILEHANDLE               => $FILEHANDLE,
                    infile_path              => $file_path,
                    remove_overlap           => 1,
                }
            );
            print {$FILEHANDLE} $PIPE . $SPACE;

            ## Collect raw total sequences and reads mapped from samtools stats and calculate the percentage. Write it to stdout
            _add_percentage_mapped_reads_from_samtools(
                {
                    FILEHANDLE   => $FILEHANDLE,
                    outfile_path => $file_path_prefix . $DOT . q{stats},
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            ## Copies file from temporary directory.
            say {$FILEHANDLE} q{## Copy file from temporary directory};
            migrate_file(
                {
                    FILEHANDLE   => $FILEHANDLE,
                    infile_path  => $file_path_prefix . $DOT . q{stats},
                    outfile_path => $outdir_path,
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
                    FILEHANDLE         => $FILEHANDLE,
                    infile_path        => $file_path,
                    outfile_path       => $file_path_prefix . $DOT . q{cram},
                    output_format      => q{cram},
                    referencefile_path => $referencefile_path,
                    with_header        => 1,
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            ## Copies file from temporary directory.
            say {$FILEHANDLE} q{## Copy file from temporary directory};
            migrate_file(
                {
                    FILEHANDLE   => $FILEHANDLE,
                    infile_path  => $file_path_prefix . $DOT . q{cram},
                    outfile_path => $outdir_path,
                }
            );
            say {$FILEHANDLE} q{wait}, $NEWLINE;
        }

        close $FILEHANDLE;

        if ( $recipe_mode == 1 ) {

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

            if (   ( $active_parameter_href->{bwa_mem_cram} )
                && ( $outfile_suffix ne q{.cram} ) )
            {

                # Required for analysisRunStatus check downstream
                my $qc_cram_path = $outfile_path_prefix . $DOT . q{cram};
                set_recipe_metafile_in_sample_info(
                    {
                        infile           => $outfile_name_prefix,
                        metafile_tag     => q{cram},
                        path             => $qc_cram_path,
                        recipe_name      => $recipe_name,
                        sample_id        => $sample_id,
                        sample_info_href => $sample_info_href,
                    }
                );

            }
            if ( $active_parameter_href->{bwa_mem_bamstats} ) {

                ## Collect QC metadata info for later use
                my $qc_stats_outfile = $outfile_path_prefix . $DOT . q{stats};
                set_recipe_outfile_in_sample_info(
                    {
                        infile           => $outfile_name_prefix,
                        path             => $qc_stats_outfile,
                        recipe_name      => q{bamstats},
                        sample_id        => $sample_id,
                        sample_info_href => $sample_info_href,
                    }
                );
            }

            if ( $bwa_binary eq q{bwa mem} ) {

                set_recipe_outfile_in_sample_info(
                    {
                        infile           => $outfile_name_prefix,
                        path             => catfile( $directory, $stderr_file ),
                        recipe_name      => $recipe_name,
                        sample_id        => $sample_id,
                        sample_info_href => $sample_info_href,
                    }
                );
            }
            if ( $bwa_binary eq q{run-bwamem} ) {

                my $qc_bwa_log = $outfile_path_prefix . $DOT . q{log.bwamem};
                set_recipe_outfile_in_sample_info(
                    {
                        infile           => $outfile_name_prefix,
                        path             => $qc_bwa_log,
                        recipe_name      => $recipe_name,
                        sample_id        => $sample_id,
                        sample_info_href => $sample_info_href,
                    }
                );
            }

            submit_recipe(
                {
                    base_command            => $profile_base_command,
                    dependency_method       => q{sample_to_sample_parallel},
                    case_id                 => $case_id,
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

sub _get_bwamem_binary {

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
            defined     => 1,
            required    => 1,
            store       => \$human_genome_reference_source,
            strict_type => 1,
        },
        human_genome_reference_version => {
            defined     => 1,
            required    => 1,
            store       => \$human_genome_reference_version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    Readonly my $GENOME_BUILD_VERSION_GRCH_PRIOR_ALTS => 37;
    Readonly my $GENOME_BUILD_VERSION_HG_PRIOR_ALTS   => 19;

    if ( $human_genome_reference_source eq q{GRCh} ) {

        # Human genome version > GRCh37
        return q{run-bwamem}
          if ( $human_genome_reference_version > $GENOME_BUILD_VERSION_GRCH_PRIOR_ALTS );

        # Human genome version <= GRCh37
        return q{bwa mem};
    }
    else {
        # hgXX build

        # Human genome version > hg19
        return q{run-bwamem}
          if ( $human_genome_reference_version > $GENOME_BUILD_VERSION_HG_PRIOR_ALTS );

        # Human genome version <= hg19
        return q{bwa mem};
    }
    return;
}

sub _add_percentage_mapped_reads_from_samtools {

## Function : Collect raw total sequences and reads mapped from samtools stats and calculate the percentage. Write it to stdout.
## Returns  :
## Arguments: $FILEHANDLE   => Filehandle to write to
##          : $outfile_path => Outfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $outfile_path;

    my $tmpl = {
        FILEHANDLE   => { defined => 1, required => 1, store => \$FILEHANDLE, },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
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
    print {$FILEHANDLE} q?print qq{percentage mapped reads:\t} . $percentage . qq{\n}?;

    # End elsif
    print {$FILEHANDLE} q?} ?;

    # End oneliner
    print {$FILEHANDLE} q?' ?;
    say   {$FILEHANDLE} q{>} . $SPACE . $outfile_path;

    return;
}

1;
