package MIP::Recipes::Analysis::Vardict;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use List::MoreUtils qw { any };
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
    our @EXPORT_OK = qw{ analysis_vardict };

}

## Constants
Readonly my $STDOUT     => q{>};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $PIPE       => q{|};
Readonly my $UNDERSCORE => q{_};

sub analysis_vardict {

## Function : Analysis recipe for varDict variant calling
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
    use MIP::Program::Variantcalling::Vardict
      qw{ vardict vardict_var2vcf_single vardict_var2vcf_paired vardict_testsomatic vardict_teststrandbias };
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
      $file_info_href->{$sample_id}{UPPSTREAM_DEPENDENCY_PROGRAM}{file_tag};
    my $outfile_tag =
      $file_info_href->{$sample_id}{$mip_program_name}{file_tag};

    my $infile_prefix  = $merged_infile_prefix . $infile_tag;
    my $outfile_prefix = $merged_infile_prefix . $outfile_tag;

    ## Get infile_suffix from baserecalibration jobid chain
    my $infile_suffix = get_file_suffix(
        {
            jobid_chain =>
              $parameter_href->{UPPSTREAM_DEPENDENCY_PROGRAM}{chain},
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
    my $infile_name  = $infile_prefix . $infile_suffix;
    my $outfile_name = $outfile_prefix . $outfile_suffix;

    ## Paths
    my $infile_path  = catfile( $insample_directory,  $infile_name );
    my $outfile_path = catfile( $outsample_directory, $outfile_name );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_name, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $sample_id,
            FILEHANDLE            => $FILEHANDLE,
            job_id_href           => $job_id_href,
            process_time          => $time,
            program_directory => catfile( $outaligner_dir, q{PROGRAM_NAME} ),
            program_name      => $program_name,
            source_environment_commands_ref => [$source_environment_cmd],
        }
    );

    ## Create input filename bam file for $infile_paths_ref in vardict module. This has to be an array where first
    ## is tumor and the second one is normal tissue.
    ##TODO: sort alphabetically
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $af_threshold       = $active_parameter_href->{af_threshold};
    my $sample_name        = $active_parameter_href->{sample_name};
    my $chrom_start        = $active_parameter_href->{chrom_start};
    my $region_start       = $active_parameter_href->{region_start};
    my $region_end         = $active_parameter_href->{region_end};
    my $segment_annotn     = $active_parameter_href->{segment_annotn};
    my $input_bed_file     = $active_parameter_href->{input_bed_file};

    say {$FILEHANDLE} q{## varDict variant calling};

    if ( $active_parameter_href->{variantcall_varDict_PE} ) {
        vardict(
            {
                af_threshold           => $af_threshold,
                FILEHANDLE             => $FILEHANDLE,
                infile_paths_ref       => $infile_path,
                out_chrom_start        => $chrom_start,
                out_region_start       => $region_start,
                out_region_end         => $region_end,
                out_segment_annotn     => $segment_annotn,
                referencefile_path     => $referencefile_path,
                sample_name            => $sample_name,
                infile_bed_region_info => $input_bed_file,
            }
        );

        print {$FILEHANDLE} $PIPE . $SPACE;

        vardict_teststrandbias(
            {
                FILEHANDLE => $FILEHANDLE,

            }
        );

        print {$FILEHANDLE} $PIPE . $SPACE;

        vardict_var2vcf_paired(
            {
                af_threshold => $af_threshold,
                FILEHANDLE   => $FILEHANDLE,
                sample_name  => $sample_name,
            }
        );

        print {$FILEHANDLE} $SPACE . $STDOUT . print{$outfile_path};

    }

    if ( $active_parameter_href->{variantcall_varDict_SE} ) {
        vardict(
            {
                af_threshold           => $af_threshold,
                FILEHANDLE             => $FILEHANDLE,
                infile_paths_ref       => $infile_path,
                out_chrom_start        => $chrom_start,
                out_region_start       => $region_start,
                out_region_end         => $region_end,
                out_segment_annotn     => $segment_annotn,
                referencefile_path     => $referencefile_path,
                sample_name            => $sample_name,
                infile_bed_region_info => $input_bed_file,
            }
        );

        print {$FILEHANDLE} $PIPE . $SPACE;

        vardict_testsomatic(
            {
                FILEHANDLE => $FILEHANDLE,

            }
        );

        print {$FILEHANDLE} $PIPE . $SPACE;

        vardict_var2vcf_single(
            {
                af_threshold => $af_threshold,
                FILEHANDLE   => $FILEHANDLE,
                sample_name  => $sample_name,
            }
        );

        print {$FILEHANDLE} $SPACE . $STDOUT . print{$outfile_path};

    }

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
                infile           => $merged_infile_prefix,
                path             => $program_outfile_path,
                program_name     => q{PROGRAM_NAME},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        my $most_complete_format_key =
          q{most_complete} . $UNDERSCORE . substr $outfile_suffix, 1;
        my $qc_metafile_path =
          catfile( $outsample_directory, $infile_prefix . $outfile_suffix );
        add_processing_metafile_to_sample_info(
            {
                metafile_tag     => $most_complete_format_key,
                path             => $qc_metafile_path,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        ## MODIY THE CHOICE OF SUB ACCORDING TO HOW YOU WANT SLURM TO PROCESSES IT AND DOWNSTREAM DEPENDENCIES
        slurm_submit_job_sample_id_dependency_add_to_sample(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_id               => $sample_id,
                sbatch_file_name        => $file_name
            }
        );
    }
    return;
}

1;
