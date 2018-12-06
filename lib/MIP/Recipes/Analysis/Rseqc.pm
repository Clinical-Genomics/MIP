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
    our @EXPORT_OK = qw{ analysis_rseqc };

}

## Constants
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_rseqc {

## Function : Rseqc analysis for RNA-seq
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id               => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name            => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

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

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_parameters get_recipe_attributes };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Qc::Rseqc
      qw{ rseqc_bam_stat rseqc_infer_experiment rseqc_inner_distance rseqc_junction_annotation rseqc_junction_saturation rseqc_read_distribution rseqc_read_duplication };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

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
    my $infile_path_prefix   = $io{in}{file_path_prefix};
    my @infile_name_prefixes = @{ $io{in}{file_name_prefixes} };
    my $infile_suffix        = $io{in}{file_suffix};
    my $infile_path          = $infile_path_prefix . $infile_suffix;
    my $recipe_mode          = $active_parameter_href->{$recipe_name};
    my $bed_file_path        = $active_parameter_href->{rseqc_transcripts_file};
    my $job_id_chain         = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
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
    my $FILEHANDLE = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $sample_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ## SHELL

    ## Rseq
    say {$FILEHANDLE} q{## Rseq infer_experiment.py};
    rseqc_infer_experiment(
        {
            bed_file_path   => $bed_file_path,
            FILEHANDLE      => $FILEHANDLE,
            infile_path     => $infile_path,
            stdoutfile_path => $outfile_path_prefix
              . $UNDERSCORE
              . q{infer_experiment}
              . $outfile_suffix,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Rseq junction_annotation.py};
    rseqc_junction_annotation(
        {
            bed_file_path        => $bed_file_path,
            FILEHANDLE           => $FILEHANDLE,
            infile_path          => $infile_path,
            outfiles_path_prefix => $outfile_path_prefix
              . $UNDERSCORE
              . q{junction_annotation},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Rseq junction_saturation.py};
    rseqc_junction_saturation(
        {
            bed_file_path        => $bed_file_path,
            FILEHANDLE           => $FILEHANDLE,
            infile_path          => $infile_path,
            outfiles_path_prefix => $outfile_path_prefix
              . $UNDERSCORE
              . q{junction_saturation},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Rseq inner_distance.py};
    rseqc_inner_distance(
        {
            bed_file_path        => $bed_file_path,
            FILEHANDLE           => $FILEHANDLE,
            infile_path          => $infile_path,
            outfiles_path_prefix => $outfile_path_prefix
              . $UNDERSCORE
              . q{inner_distance},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Rseq bam_stat.py};
    rseqc_bam_stat(
        {
            FILEHANDLE      => $FILEHANDLE,
            infile_path     => $infile_path,
            stdoutfile_path => $outfile_path_prefix
              . $UNDERSCORE
              . q{bam_stat}
              . $outfile_suffix,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Rseq read_distribution.py};
    rseqc_read_distribution(
        {
            bed_file_path   => $bed_file_path,
            FILEHANDLE      => $FILEHANDLE,
            infile_path     => $infile_path,
            stdoutfile_path => $outfile_path_prefix
              . $UNDERSCORE
              . q{read_distribution}
              . $outfile_suffix,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Rseqc read_duplication};
    rseqc_read_duplication(
        {
            FILEHANDLE           => $FILEHANDLE,
            infile_path          => $infile_path,
            outfiles_path_prefix => $outfile_path_prefix
              . $UNDERSCORE
              . q{read_duplication},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    close $FILEHANDLE;

    if ( $recipe_mode == 1 ) {

        submit_recipe(
            {
                dependency_method       => q{sample_to_sample},
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                job_id_chain            => $job_id_chain,
                recipe_file_path        => $recipe_file_path,
                sample_id               => $sample_id,
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return;
}

1;
