package MIP::Recipes::Sambamba_depth;

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
use File::Spec::Functions qw{ catdir catfile devnull };

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_sambamba_depth };

}

##Constants
Readonly my $NEWLINE => qq{\n};

sub analysis_sambamba_depth {

## Function : Generate coverage bed outfile for each individual.
## Returns  : ""
## Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $sample_id, $program_name, family_id, $temp_directory, $outaligner_dir
##          : $parameter_href             => Parameter hash {REF}
##          : $active_parameter_href      => Active parameters for this analysis hash {REF}
##          : $sample_info_href           => Info on samples and family hash {REF}
##          : $file_info_href             => The file_info hash {REF}
##          : $infile_lane_prefix_href    => Infile(s) without the ".ending" {REF}
##          : $job_id_href                => Job id hash {REF}
##          : $sample_id                  => Sample id
##          : $outaligner_dir             => Outaligner_dir used in the analysis
##          : $program_name               => Program name
##          : $family_id                  => Family id
##          : $temp_directory             => Temporary directory
##          : $outaligner_dir             => Outaligner_dir used in the analysis

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id;
    my $temp_directory;
    my $outaligner_dir;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
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
        sample_id => {
            required    => 1,
            defined     => 1,
            default     => 1,
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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory
        },
        outaligner_dir => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use List::MoreUtils qw { any };
    use MIP::Script::Setup_script qw{ setup_script};
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::IO::Files qw{ migrate_file};
    use MIP::Program::Alignment::Sambamba qw{ sambamba_depth };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_dead_end };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Assign directories
    my $insample_directory = catdir( $active_parameter_href->{outdata_dir},
        $sample_id, $outaligner_dir );
    my $outsample_directory = catdir( $active_parameter_href->{outdata_dir},
        $sample_id, $outaligner_dir, q{coveragereport} );

    ## Add merged infile name after merging all BAM files per sample_id
    my $infile = $file_info_href->{$sample_id}{merge_infile};

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$sample_id}{pgatk_baserecalibration}{file_tag};
    my $outfile_tag =
      $file_info_href->{$sample_id}{$mip_program_name}{file_tag};

    my $infile_prefix       = $infile . $infile_tag;
    my $outfile_prefix      = $infile . $outfile_tag;
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = $file_path_prefix . $outfile_tag;

    ## Get infile_suffix from baserecalibration jobid chain
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
            jobid_chain    => $parameter_href->{pgatk_baserecalibration}{chain},
        }
    );
    my $outfile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{outfile_suffix},
            program_name   => $mip_program_name,
        }
    );

    # q{.bam} -> ".b*" for getting index as well)
    my $infile_path = catfile( $insample_directory,
        $infile_prefix . substr( $infile_suffix, 0, 2 ) . q{*} );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $sample_id,
            program_name          => $program_name,
            program_directory =>
              catfile( lc($outaligner_dir), q{coveragereport} ),
            core_number    => $core_number,
            process_time   => $time,
            temp_directory => $temp_directory,
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
    my $sambamba_filter = q?'mapping_quality >= ?
      . $active_parameter_href->{sambamba_depth_mapping_quality} . q? ?;

    #Do not include duplicates in coverage calculation
    if ( $active_parameter_href->{sambamba_depth_noduplicates} ) {

        $sambamba_filter .= q{and not duplicate };
    }

    #Do not include failed quality control reads in coverage calculation
    if ( $active_parameter_href->{sambamba_depth_quality_control} ) {

        $sambamba_filter .= q{and not failed_quality_control};
    }
    $sambamba_filter .= q?'?;

    sambamba_depth(
        {
            depth_cutoffs_ref =>
              \@{ $active_parameter_href->{sambamba_depth_cutoffs} },
            infile_path      => $file_path_prefix . $infile_suffix,
            outfile_path     => $outfile_path_prefix . $outfile_suffix,
            mode             => $active_parameter_href->{sambamba_depth_mode},
            fix_mate_overlap => 1,
            min_base_quality =>
              $active_parameter_href->{sambamba_depth_base_quality},
            filter     => $sambamba_filter,
            region     => $active_parameter_href->{sambamba_depth_bed},
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            infile_path  => $outfile_path_prefix . $outfile_suffix,
            outfile_path => $outsample_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;
    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {

        my $qc_sambamba_path =
          catfile( $outsample_directory, $outfile_prefix . $outfile_suffix );
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                sample_id        => $sample_id,
                program_name     => $program_name,
                infile           => $infile,
                path             => $qc_sambamba_path,
            }
        );

        slurm_submit_job_sample_id_dependency_dead_end(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                family_id               => $family_id,
                sample_id               => $sample_id,
                path                    => $job_id_chain,
                sbatch_file_name        => $file_path,
                log                     => $log,
            }
        );
    }

    return;
}

1;
