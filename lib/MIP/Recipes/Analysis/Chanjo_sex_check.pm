package MIP::Recipes::Analysis::Chanjo_sex_check;

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
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_chanjo_sex_check };

}

##Constants
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_chanjo_sex_check {

## analysis_chanjo_sex_check

## Function : Predicts gender from BAM files.
## Returns  : ""
## Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $sample_id, family_id, $insample_directory, $outsample_directory, $outaligner_dir, $program_name
##          : $parameter_href             => Parameter hash {REF}
##          : $active_parameter_href      => Active parameters for this analysis hash {REF}
##          : $sample_info_href           => Info on samples and family hash {REF}
##          : $file_info_href             => The file_info hash {REF}
##          : $infile_lane_prefix_href    => Infile(s) without the ".ending" {REF}
##          : $job_id_href                => Job id hash {REF}
##          : $sample_id                  => Sample id
##          : $insample_directory         => In sample directory
##          : $outsample_directory        => Out sample directory
##          : $program_name               => Program name
##          : $family_id                  => Family id
##          : $outaligner_dir             => Outaligner_dir used in the analysis

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $sample_id;
    my $insample_directory;
    my $outsample_directory;
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
            strict_type => 1,
            store       => \$sample_id
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use List::MoreUtils qw { any };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Program::Alignment::Chanjo qw{ chanjo_sex };
    use MIP::QC::Record
      qw{ add_program_outfile_to_sample_info add_program_metafile_to_sample_info };
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

    ## Files
    my $infile_name  = $infile_prefix . $infile_suffix;
    my $outfile_name = $outfile_prefix . $outfile_suffix;

    ## Paths
    my $infile_path  = catfile( $insample_directory,  $infile_name );
    my $outfile_path = catfile( $outsample_directory, $outfile_name );
    my $log_file_path = catfile( $outsample_directory,
        $infile_prefix . $UNDERSCORE . q{chanjo_sexcheck.log} );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_name, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $sample_id,
            program_name          => $program_name,
            program_directory =>
              catfile( lc($outaligner_dir), q{coveragereport} ),
            core_number  => $core_number,
            process_time => $time,
        }
    );

    ## chanjo_sexcheck
    say {$FILEHANDLE} q{## Predicting sex from alignment};

    ## Get parameters

    my $chr_prefix;

    # If element is part of array
    if ( any { $_ eq q{chrX} } @{ $file_info_href->{contigs_size_ordered} } ) {
        $chr_prefix = q{chr};
    }
    chanjo_sex(
        {
            infile_path  => $infile_path,
            outfile_path => $outfile_path,
            log_level    => $active_parameter_href->{chanjo_sexcheck_log_level},
            log_file_path => $log_file_path,
            chr_prefix    => $chr_prefix,
            FILEHANDLE    => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;
    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                sample_id        => $sample_id,
                program_name     => q{chanjo_sexcheck},
                infile           => $merged_infile_prefix,
                outdirectory     => $outsample_directory,
                outfile          => $outfile_name,
            }
        );
        add_program_metafile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                sample_id        => $sample_id,
                program_name     => q{chanjo_sexcheck},
                infile           => $merged_infile_prefix,
                metafile_tag     => q{log},
                directory        => $outsample_directory,
                file => $infile_prefix . $UNDERSCORE . q{chanjo_sexcheck.log},
            }
        );
        slurm_submit_job_sample_id_dependency_dead_end(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                family_id               => $family_id,
                sample_id               => $sample_id,
                path                    => $job_id_chain,
                sbatch_file_name        => $file_name,
                log                     => $log,
            }
        );
    }

    return;
}

1;
