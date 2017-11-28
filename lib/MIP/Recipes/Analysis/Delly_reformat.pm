package MIP::Recipes::Analysis::Delly_reformat;

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
    our @EXPORT_OK = qw{ analysis_delly_reformat };

}

## Constants
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_delly_reformat {

## Function : Merge, regenotype, and filter using Delly
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => Variant call type
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $insample_directory      => In sample directory
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outsample_directory     => Out sample directory
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $reference_dir           => MIP reference directory
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter


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
    my $call_type;
    my $family_id;
    my $outaligner_dir;
    my $reference_dir;
    my $temp_directory;
    my $xargs_file_counter;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        call_type =>
          { default => q{SV}, strict_type => 1, store => \$call_type, },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        insample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$insample_directory,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
        },
        outsample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outsample_directory,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        reference_dir_ref => {
            default     => $arg_href->{active_parameter_href}{reference_dir},
            strict_type => 1,
            store       => \$reference_dir_ref,
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        temp_directory_ref => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref,
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ### SET IMPORT IN ALPHABETIC ORDER
    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::PATH::TO::PROGRAMS qw{ COMMANDS_SUB };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::QC::Record
      qw{ add_program_metafile_to_sample_info add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};

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
      $file_info_href->{$sample_id}{UPPSTREAM_DEPENDENCY_PROGRAM}{file_tag};
    my $outfile_tag =
      $file_info_href->{$sample_id}{$mip_program_name}{file_tag};

    my $infile_prefix  = $merged_infile_prefix . $infile_tag;
    my $outfile_prefix = $merged_infile_prefix . $outfile_tag;

    ## Get infile_suffix from baserecalibration jobid chain
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
            jobid_chain =>
              $parameter_href->{UPPSTREAM_DEPENDENCY_PROGRAM}{chain},
        }
    );
    my $outfile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{outfile_suffix},
            program_name   => $mip_program_name,
        }
    );

    ## Files
    my $infile_name  = $infile_prefix . $infile_suffix;
    my $outfile_name = $outfile_prefix . $outfile_suffix;

    ## Paths
    my $infile_path  = catfile( $insample_directory,  $infile_name );
    my $outfile_path = catfile( $outsample_directory, $outfile_name );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_name, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $sample_id,
            program_name          => $program_name,
            program_directory => catfile( $outaligner_dir, q{PROGRAM_NAME} ),
            core_number       => $core_number,
            process_time      => $time,
        }
    );

###############################
###RECIPE TOOL COMMANDS HERE###
###############################

    ## Close FILEHANDLES
    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $mip_program_mode == 1 ) {

        my $program_outfile_path = catfile( $outsample_directory,
            $outfile_prefix . $UNDERSCORE . q{ENDING} );
        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                sample_id        => $sample_id,
                program_name     => q{PROGRAM_NAME},
                infile           => $merged_infile_prefix,
                path             => $program_outfile_path,
            }
        );

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

        ## MODIY THE CHOICE OF SUB ACCORDING TO HOW YOU WANT SLURM TO PROCESSES IT AND DOWNSTREAM DEPENDENCIES
        slurm_submit_job_sample_id_dependency_add_to_sample(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                family_id               => $family_id,
                sample_id               => $sample_id,
                path                    => $job_id_chain,
                log                     => $log,
                sbatch_file_name        => $file_path
            }
        );
    }
    return;
}

1;
