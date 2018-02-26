package MIP::Recipes::Analysis::Samtools_subsample_MT;

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
    our @EXPORT_OK = qw{ analysis_samtools_subsample_MT };

}

## Constants
Readonly my $BACKTICK                    => q{`};
Readonly my $DOT                         => q{.};
Readonly my $NEWLINE                     => qq{\n};
Readonly my $MAX_LIMIT_SEED              => 100;
Readonly my $PIPE                        => q{|};
Readonly my $SAMTOOLS_UNMAPPED_READ_FLAG => 4;
Readonly my $SPACE                       => q{ };
Readonly my $UNDERSCORE                  => q{_};

sub analysis_samtools_subsample_MT {

## Function : Creates a BAM file containing a subset of the MT alignments
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $insample_directory      => In sample directory
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outsample_directory     => Out sample directory
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Program::Alignment::Samtools
      qw{ samtools_depth samtools_index samtools_view };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_dead_end };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$sample_id}{pgatk_baserecalibration}{file_tag};
    my $outfile_tag =
      $file_info_href->{$sample_id}{$mip_program_name}{file_tag};

    my $infile_prefix  = $merged_infile_prefix . $infile_tag;
    my $outfile_prefix = $merged_infile_prefix . $outfile_tag;

    ## Get infile_suffix from baserecalibration jobid chain
    my $infile_suffix = get_file_suffix(
        {
            jobid_chain    => $parameter_href->{pgatk_baserecalibration}{chain},
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
        }
    );
    my $outfile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            program_name   => $mip_program_name,
            suffix_key     => q{outfile_suffix},
        }
    );

    ## Files
    my $infile_name  = $infile_prefix . q{_MT} . $infile_suffix;
    my $outfile_name = $outfile_prefix . $outfile_suffix;

    ## Paths
    my $infile_path  = catfile( $insample_directory,  $infile_name );
    my $outfile_path = catfile( $outsample_directory, $outfile_name );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_name, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $sample_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            process_time                    => $time,
            program_directory               => catfile($outaligner_dir),
            program_name                    => $program_name,
            source_environment_commands_ref => [$source_environment_cmd],
        }
    );

    ## Set up seed and fraction combination
    say $FILEHANDLE q{## Creating subsample filter for samtools view};

    ## Get average coverage over MT bases
    print $FILEHANDLE q{MT_COVERAGE=} . $BACKTICK;

    # Get depth per base
    samtools_depth(
        {
            FILEHANDLE         => $FILEHANDLE,
            infile_path        => $infile_path,
            max_depth_treshold => q{500000},
        }
    );

    # Pipe to AWK
    print $FILEHANDLE $PIPE . $SPACE;

    # Add AWK statment for calculation of avgerage coverage
    print $FILEHANDLE _awk_calculate_average_coverage();

    # Close statment
    say $FILEHANDLE $BACKTICK;

    ## Get subsample depth
    my $mt_subsample_depth =
      $active_parameter_href->{samtools_subsample_mt_depth};

    ## Get random seed
    my $seed = int rand $MAX_LIMIT_SEED;

    ## Add seed to fraction for ~100x
    # Create bash variable
    say $FILEHANDLE q{SEED_FRACTION=}

      # Open statment
      . $BACKTICK

      # Lauch perl and print
      . q?perl -e "print ?

      # Add the random seed number to..
      . $seed . q{ + }

     # ...the subsample fraction, consisting of the desired subsample coverag...
      . $mt_subsample_depth

      # ...divided by the starting coverage
      . q? / $MT_COVERAGE"?

      # Close statment
      . $BACKTICK . $NEWLINE;

    ## Filter the bam file to only include a subset of reads that maps to the MT
    say $FILEHANDLE q{## Filter the BAM file};
    samtools_view(
        {
            exclude_reads_with_these_flags => $SAMTOOLS_UNMAPPED_READ_FLAG,
            FILEHANDLE                     => $FILEHANDLE,
            fraction                       => q{"$SEED_FRACTION"},
            infile_path                    => $infile_path,
            outfile_path                   => $outfile_path,
            with_header                    => 1,
        }
    );
    say $FILEHANDLE $NEWLINE;

    ## Index new bam file
    say $FILEHANDLE q{## Index the subsampled BAM file};
    samtools_index(
        {
            bai_format  => 1,
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $outfile_path,
        }
    );
    say $FILEHANDLE $NEWLINE;

    ## Close FILEHANDLES
    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $mip_program_mode == 1 ) {

        my $program_outfile_path = catfile( $outsample_directory,
            $outfile_prefix . $DOT . q{bam} );
        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                infile           => $merged_infile_prefix,
                path             => $program_outfile_path,
                program_name     => q{samtools_subsample_mt},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        slurm_submit_job_sample_id_dependency_dead_end(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_id               => $sample_id,
                sbatch_file_name        => $file_name,
            }
        );
    }
    return;
}

sub _awk_calculate_average_coverage {

## Function : Writes an awk expression to an open filehandle. The awk expression calculates the average coverage based on input from samtools depth and prints it.
## Returns  : $awk_statment

    my $awk_statment =

      # Start awk
      q?awk '?

      # Sum the coverage data for each base ()
      . q?{cov += $3}?

      # Add end rule
      . q?END?

      # Divide the total coverage sum with the number of covered
      # bases (rows of output from samtools depth),
      # stored in the awk built in "NR"
      . q?{ if (NR > 0) print cov / NR }'?;

    return $awk_statment;
}

1;
