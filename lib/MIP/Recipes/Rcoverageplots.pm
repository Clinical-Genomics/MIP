package MIP::Recipes::Rcoverageplots;

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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_rcoverageplots };

}

##Constants
Readonly my $AMPERSAND  => q{&};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_rcoverageplots {

## analysis_rcoverageplots

## Function : Generates sbatch scripts for R scripts: 1. covplots_genome.R 2. covplots_exome.R; on files generated from bedtools genomecov.
## Returns  :
## Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $lane_href, $infile_lane_prefix_href, $job_id_href, $sample_id, $insample_directory, $outsample_directory, $program_name, $family_id, $temp_directory, $outaligner_dir,
##          : $parameter_href          => Parameter hash {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $lane_href               => The lane info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $sample_id               => Sample id
##          : $insample_directory      => In sample directory
##          : $outsample_directory     => Out sample directory
##          : $program_name            => Program name
##          : $family_id               => Family id
##          : $temp_directory          => Temporary directory
##          : $outaligner_dir          => Outaligner_dir used in the analysis

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
    my $lane_href;
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
        lane_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$lane_href
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
            store       => \$insample_directory,
        },
        outsample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outsample_directory,
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
            store       => \$family_id,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_merged_infile_prefix };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_dead_end };
    use MIP::Program::Qc::Rcoverageplots qw{ covplots_genome };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};

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
      $file_info_href->{$sample_id}{pbedtools_genomecov}{file_tag};
    my $outfile_tag =
      $file_info_href->{$sample_id}{pgatk_baserecalibration}{file_tag};

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $sample_id,
            program_name          => $program_name,
            program_directory => catfile( $outaligner_dir, q{coveragereport} ),
            temp_directory    => $temp_directory,
        }
    );

    my $infile_path =
      catfile( $insample_directory, $merged_infile_prefix . $infile_tag );
    covplots_genome(
        {
            infile_path => $infile_path,
            sample_id   => $sample_id,
            max_coverage_depth =>
              $active_parameter_href->{bedtools_genomecov_max_coverage},
            outdirectory_path => $outsample_directory,
            FILEHANDLE        => $FILEHANDLE,
        }
    );

    say $FILEHANDLE $AMPERSAND, $NEWLINE;
    say $FILEHANDLE q{wait}, $NEWLINE;

    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {

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
