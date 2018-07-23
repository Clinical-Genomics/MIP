package MIP::Recipes::Analysis::Gatk_asereadcounter;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
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

## MIPs lib/
use MIP::Check::Cluster qw{ check_max_core_number };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_gatk_asereadcounter };

}

## Constants
Readonly my $ASTERIX => q{*};
Readonly my $NEWLINE => qq{\n};

sub analysis_gatk_asereadcounter {

## Function : Gatk asereadcounter analysis for rna recipe
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
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
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        family_id_ref => {
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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Pedigree qw{ create_fam_file gatk_pedigree_flag };
    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::Program::Variantcalling::Gatk qw{ gatk_asereadcounter };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };

    ## Constants
    Readonly my $JAVA_MEMORY_ALLOCATION => 8;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $program_mode = $active_parameter_href->{$program_name};

    ## Alias
    my $analysis_type = \$active_parameter_href->{analysis_type}{$sample_id};
    my $job_id_chain  = $parameter_href->{$program_name}{chain};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Assign directories
    my $insample_directory = catdir( $active_parameter_href->{outdata_dir},
        $sample_id, $active_parameter_href->{outaligner_dir} );
    my $outsample_directory = catdir( $active_parameter_href->{outdata_dir},
        $sample_id, $active_parameter_href->{outaligner_dir} );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $sample_id,
            FILEHANDLE            => $FILEHANDLE,
            job_id_href           => $job_id_href,
            log                   => $log,
            process_time          => $time,
            program_directory =>
              catfile( $active_parameter_href->{outaligner_dir}, q{gatk} ),
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    # Division by X according to the java heap
    $core_number = floor(
        $active_parameter_href->{node_ram_memory} / $JAVA_MEMORY_ALLOCATION );

    ## Limit number of cores requested to the maximum number of cores available per node
    $core_number = check_max_core_number(
        {
            core_number_requested => $core_number,
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    ## Assign directories
    # For ".fam" file
    my $outfamily_file_directory =
      catdir( $active_parameter_href->{outdata_dir}, $family_id );

    ## Used downstream
    $parameter_href->{$program_name}{$sample_id}{indirectory} =
      $outsample_directory;

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$sample_id}{gatk_baserecalibration}{file_tag};
    my $outfile_tag =
      $file_info_href->{$sample_id}{$program_name}{file_tag};
    my $sitesfile_tag =
      $file_info_href->{$sample_id}{gatk_variantfiltration}{file_tag};

    ## Files
    my $infile_prefix    = $merged_infile_prefix . $infile_tag;
    my $outfile_prefix   = $merged_infile_prefix . $outfile_tag;
    my $sitesfile_prefix = $merged_infile_prefix . $sitesfile_tag;

    ## Paths
    my $file_path_prefix      = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix   = catfile( $temp_directory, $outfile_prefix );
    my $sitesfile_path_prefix = catfile( $temp_directory, $sitesfile_prefix );

    ## Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            jobid_chain    => $parameter_href->{gatk_baserecalibration}{chain},
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
        }
    );

    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            file_suffix    => $parameter_href->{$program_name}{outfile_suffix},
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{ase_file_suffix},
        }
    );

    my $sitesfile_suffix = get_file_suffix(
        {
            jobid_chain    => $parameter_href->{gatk_variantfiltration}{chain},
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => catfile(
                $insample_directory,
                $infile_prefix . substr( $infile_suffix, 0, 2 ) . $ASTERIX
            ),
            outfile_path => $temp_directory,
        }
    );
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => catfile(
                $insample_directory,
                $sitesfile_prefix
                  . substr( $sitesfile_suffix, 0, 2 )
                  . $ASTERIX
            ),
            outfile_path => $temp_directory,
        }
    );
    say {$FILEHANDLE} q{wait} . $NEWLINE;

    ## GATK ASEReadCounter
    say {$FILEHANDLE} q{## GATK ASEReadCounter};

    ## Set file paths
    my $infile_path  = $file_path_prefix . $infile_suffix;
    my $outfile_path = $outfile_path_prefix . $outfile_suffix;

    gatk_asereadcounter(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $infile_path,
            java_jar    => catfile(
                $active_parameter_href->{gatk_path},
                q{GenomeAnalysisTK.jar},
            ),
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            logging_level     => $active_parameter_href->{gatk_logging_level},
            memory_allocation => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
            outfile_path      => $outfile_path,
            referencefile_path =>
              $active_parameter_href->{human_genome_reference},
            gatk_sites_vcffile => $sitesfile_path_prefix . $sitesfile_suffix,
            temp_directory     => $temp_directory,
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
    say {$FILEHANDLE} q{wait};
    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        my $program_outfile_path =
          catfile( $outsample_directory, $outfile_prefix . $outfile_suffix );

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                infile           => $infile_prefix,
                path             => $program_outfile_path,
                program_name     => q{gatk_asereadcounter},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        slurm_submit_job_sample_id_dependency_add_to_sample(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_id               => $sample_id,
                sbatch_file_name        => $file_path
            }
        );
    }
    return;
}
1;
