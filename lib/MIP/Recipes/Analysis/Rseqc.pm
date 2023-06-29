package MIP::Recipes::Analysis::Rseqc;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename fileparse };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use POSIX qw{ floor };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_rseqc };

}

sub analysis_rseqc {

## Function : Rseqc analysis for RNA-seq
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
    my $temp_directory;
    my $xargs_file_counter;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id_ref => {
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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr{ ^\d+$ }xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Rseqc
      qw{ rseqc_bam_stat rseqc_infer_experiment rseqc_inner_distance rseqc_junction_annotation rseqc_junction_saturation rseqc_read_distribution rseqc_read_duplication };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Reference qw{ get_transcript_annotation_file_path };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
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
    my $infile_path_prefix        = $io{in}{file_path_prefix};
    my @infile_name_prefixes      = @{ $io{in}{file_name_prefixes} };
    my $infile_suffix             = $io{in}{file_suffix};
    my $infile_path               = $infile_path_prefix . $infile_suffix;
    my $transcript_annotation_bed = get_transcript_annotation_file_path(
        {
            active_parameter_href => $active_parameter_href,
            file_format           => q{bed},
            file_info_href        => $file_info_href,
        }
    );
    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $recipe{job_id_chain},
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                outdata_dir            => $active_parameter_href->{outdata_dir},
                file_name_prefixes_ref => \@infile_name_prefixes,
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
                temp_directory         => $temp_directory,
            }
        )
    );
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe{core_number},
            directory_id          => $sample_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe{memory},
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            temp_directory        => $temp_directory,
        }
    );

    ## SHELL

    ## Rseq
    say {$filehandle} q{## Rseq infer_experiment.py};
    rseqc_infer_experiment(
        {
            bed_file_path   => $transcript_annotation_bed,
            filehandle      => $filehandle,
            infile_path     => $infile_path,
            stdoutfile_path => $outfile_path_prefix
              . $UNDERSCORE
              . q{infer_experiment}
              . $outfile_suffix,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Rseq junction_annotation.py};
    rseqc_junction_annotation(
        {
            bed_file_path        => $transcript_annotation_bed,
            filehandle           => $filehandle,
            infile_path          => $infile_path,
            outfiles_path_prefix => $outfile_path_prefix . $UNDERSCORE . q{junction_annotation},
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Rseq junction_saturation.py};
    rseqc_junction_saturation(
        {
            bed_file_path        => $transcript_annotation_bed,
            filehandle           => $filehandle,
            infile_path          => $infile_path,
            outfiles_path_prefix => $outfile_path_prefix . $UNDERSCORE . q{junction_saturation},
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Rseq inner_distance.py};
    rseqc_inner_distance(
        {
            bed_file_path        => $transcript_annotation_bed,
            filehandle           => $filehandle,
            infile_path          => $infile_path,
            outfiles_path_prefix => $outfile_path_prefix . $UNDERSCORE . q{inner_distance},
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Rseq bam_stat.py};
    rseqc_bam_stat(
        {
            filehandle      => $filehandle,
            infile_path     => $infile_path,
            stdoutfile_path => $outfile_path_prefix . $UNDERSCORE . q{bam_stat} . $outfile_suffix,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Rseq read_distribution.py};
    rseqc_read_distribution(
        {
            bed_file_path   => $transcript_annotation_bed,
            filehandle      => $filehandle,
            infile_path     => $infile_path,
            stdoutfile_path => $outfile_path_prefix
              . $UNDERSCORE
              . q{read_distribution}
              . $outfile_suffix,
        }
    );
    say {$filehandle} $NEWLINE;

    close $filehandle;

    if ( $recipe{mode} == 1 ) {

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                dependency_method                 => q{sample_to_island},
                case_id                           => $case_id,
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
                log                               => $log,
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

1;
