package MIP::Recipes::Analysis::Samtools_subsample_mt;

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
use List::MoreUtils qw{ first_value };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $BACKTICK $DOT $LOG_NAME $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.16;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_samtools_subsample_mt };

}

## Constants
Readonly my $MAX_LIMIT_SEED              => 100;
Readonly my $SAMTOOLS_UNMAPPED_READ_FLAG => 4;

sub analysis_samtools_subsample_mt {

## Function : Creates a BAM file containing a subset of the MT alignments
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $profile_base_command;
    my $recipe_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $case_id;

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
    use MIP::Language::Awk qw{ awk };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Program::Samtools qw{ samtools_index samtools_view };
    use MIP::Program::Bedtools qw{ bedtools_genomecov };
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
            recipe_name    => $recipe_name,
            stream         => q{in},
        }
    );
    my $infile_name_prefix   = $io{in}{file_name_prefix};
    my @infile_name_prefixes = @{ $io{in}{file_name_prefixes} };
    my @infile_paths         = @{ $io{in}{file_paths} };

    ## Find Mitochondrial contig infile_path
    my $infile_path = first_value { / $infile_name_prefix [.]M|chrM /sxm } @infile_paths;

    if ( not $infile_path ) {

        $log->warn(
qq{Mitochondrial contig is not part of analysis contig set - skipping $recipe_name}
        );
        return 1;
    }

    my $job_id_chain = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $mt_subsample_depth = $active_parameter_href->{samtools_subsample_mt_depth};
    my $recipe_mode        = $active_parameter_href->{$recipe_name};
    my %recipe_resource    = get_recipe_resources(
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
                file_name_prefixes_ref => \@infile_name_prefixes,
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );

    my $outfile_name_prefix = $io{out}{file_name_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my $outfile_path        = catfile( $outfile_path_prefix . $outfile_suffix );

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
            memory_allocation               => $recipe_resource{memory},
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
        }
    );

    ### SHELL:

    ## Set up seed and fraction combination
    say {$filehandle} q{## Creating subsample filter for samtools view};

    ## Get average coverage over MT bases
    print {$filehandle} q{MT_COVERAGE=} . $BACKTICK;

    # Get depth per base
    bedtools_genomecov(
        {
            depth_each_position => 1,
            filehandle          => $filehandle,
            bam_infile_path     => $infile_path,
        }
    );

    # Pipe to AWK
    print {$filehandle} $PIPE . $SPACE;

    # Add AWK statment for calculation of avgerage coverage
    my $awk_statment = _awk_calculate_average_coverage();
    awk(
        {
            filehandle => $filehandle,
            statement  => $awk_statment,
        }
    );

    # Close statment
    say {$filehandle} $BACKTICK;

    ## Get random seed
    my $seed = int rand $MAX_LIMIT_SEED;

    ## Add seed to fraction for ~100x
    # Create bash variable
    say {$filehandle} q{SEED_FRACTION=}

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
    say {$filehandle} q{## Filter the BAM file};
    samtools_view(
        {
            exclude_reads_with_these_flags => $SAMTOOLS_UNMAPPED_READ_FLAG,
            filehandle                     => $filehandle,
            fraction                       => q{"$SEED_FRACTION"},
            infile_path                    => $infile_path,
            outfile_path                   => $outfile_path,
            with_header                    => 1,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Index new bam file
    say {$filehandle} q{## Index the subsampled BAM file};
    samtools_index(
        {
            bai_format  => 1,
            filehandle  => $filehandle,
            infile_path => $outfile_path,
        }
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

        set_file_path_to_store(
            {
                format           => q{bam},
                id               => $sample_id,
                path             => $outfile_path,
                path_index       => $outfile_path . $DOT . q{bai},
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command         => $profile_base_command,
                case_id              => $case_id,
                dependency_method    => q{sample_to_island},
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

sub _awk_calculate_average_coverage {

## Function : Writes an awk expression to an open filehandle. The awk expression calculates the average coverage based on input from samtools depth and prints it.
## Returns  : $awk_statment

    my $awk_statment =

      # Sum the coverage data for each base ()
      q?{cov += $3}?

      # Add end rule
      . q?END?

      # Divide the total coverage sum with the number of covered
      # bases (rows of output from samtools depth),
      # stored in the awk built in "NR"
      . q?{ if (NR > 0) print cov / NR }?;

    return $awk_statment;
}

1;
