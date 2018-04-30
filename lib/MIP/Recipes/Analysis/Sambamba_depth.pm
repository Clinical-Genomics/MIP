package MIP::Recipes::Analysis::Sambamba_depth;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile devnull };
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
    our @EXPORT_OK = qw{ analysis_sambamba_depth };

}

## Constants
Readonly my $ASTERIX      => q{*};
Readonly my $NEWLINE      => qq{\n};
Readonly my $SINGLE_QUOTE => q{'};
Readonly my $SPACE        => q{ };

sub analysis_sambamba_depth {

## Function : Generate coverage bed outfile for each individual.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => The file_info hash {REF}
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
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        sample_id => {
            default     => 1,
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

    use List::MoreUtils qw { any };
    use MIP::Script::Setup_script qw{ setup_script};
    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file};
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_dead_end };
    use MIP::Program::Alignment::Sambamba qw{ sambamba_depth };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my ( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Get merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

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

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$sample_id}{pgatk_baserecalibration}{file_tag};
    my $outfile_tag =
      $file_info_href->{$sample_id}{$mip_program_name}{file_tag};
    my $infile_prefix  = $merged_infile_prefix . $infile_tag;
    my $outfile_prefix = $merged_infile_prefix . $outfile_tag;

    ## Files
    my $infile_name  = $infile_prefix . $infile_suffix;
    my $outfile_name = $outfile_prefix . $outfile_suffix;

    ## Paths
    # q{.bam} -> ".b*" for getting index as well)
    my $infile_path = catfile( $insample_directory,
        $infile_prefix . substr( $infile_suffix, 0, 2 ) . $ASTERIX );
    my $file_path_prefix = catfile( $temp_directory, $infile_prefix );
    my $outfile_path     = catfile( $temp_directory, $outfile_name );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $sample_id,
            FILEHANDLE            => $FILEHANDLE,
            job_id_href           => $job_id_href,
            process_time          => $time,
            program_directory => catfile( $outaligner_dir, q{coveragereport} ),
            program_name      => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $infile_path,
            outfile_path => $temp_directory
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## sambamba_depth
    say {$FILEHANDLE} q{## Annotating bed from alignment};

    ## Get parameters
    my $sambamba_filter =
        $SINGLE_QUOTE
      . q{mapping_quality >=}
      . $SPACE
      . $active_parameter_href->{sambamba_depth_mapping_quality}
      . $SPACE;

    #Do not include duplicates in coverage calculation
    if ( $active_parameter_href->{sambamba_depth_noduplicates} ) {

        $sambamba_filter .= q{and not duplicate} . $SPACE;
    }

    #Do not include failed quality control reads in coverage calculation
    if ( $active_parameter_href->{sambamba_depth_quality_control} ) {

        $sambamba_filter .= q{and not failed_quality_control};
    }
    $sambamba_filter .= $SINGLE_QUOTE;

    sambamba_depth(
        {
            depth_cutoffs_ref =>
              \@{ $active_parameter_href->{sambamba_depth_cutoffs} },
            FILEHANDLE       => $FILEHANDLE,
            filter           => $sambamba_filter,
            fix_mate_overlap => 1,
            infile_path      => $file_path_prefix . $infile_suffix,
            min_base_quality =>
              $active_parameter_href->{sambamba_depth_base_quality},
            mode         => $active_parameter_href->{sambamba_depth_mode},
            outfile_path => $outfile_path,
            region       => $active_parameter_href->{sambamba_depth_bed},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $outfile_path,
            outfile_path => $outsample_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;
    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {

        my $qc_sambamba_path = catfile( $outsample_directory, $outfile_name );

        add_program_outfile_to_sample_info(
            {
                infile           => $merged_infile_prefix,
                path             => $qc_sambamba_path,
                program_name     => $program_name,
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
                sbatch_file_name        => $file_path,
            }
        );
    }
    return;
}

1;
