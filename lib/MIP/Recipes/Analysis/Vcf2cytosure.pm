package MIP::Recipes::Analysis::Vcf2cytosure;

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
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_vcf2cytosure };

}

## Constants
Readonly my $AMPERSAND    => q{&};
Readonly my $ASTERISK     => q{*};
Readonly my $DOT          => q{.};
Readonly my $NEWLINE      => qq{\n};
Readonly my $SINGLE_QUOTE => q{'};
Readonly my $SPACE        => q{ };
Readonly my $SV_LENGTH    => 3000;
Readonly my $UNDERSCORE   => q{_};

sub analysis_vcf2cytosure {

## Function : Convert VCF with structural variations to the “.CGH” format used by the CytoSure Interpret Software
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $bin_size                => Bin size
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outfamily_directory     => Out family directory
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $outfamily_directory;
    my $parameter_href;
    my $program_name;
    my $sample_info_href;

    ## Default(s)
    my $bin_size;
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
        bin_size => {
            default     => $arg_href->{active_parameter_href}{tiddit_bin_size},
            strict_type => 1,
            store       => \$bin_size
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
        outfamily_directory => {
            defined     => 1,
            required    => 1,
            store       => \$outfamily_directory,
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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Get::Parameter qw{ get_module_parameters get_program_parameters };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Program::Variantcalling::Vcf2cytosure qw{ vcf2cytosure_convert };
    use MIP::Processmanagement::Processes qw{ print_wait };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_family_dead_end };
    use MIP::Program::Variantcalling::Bcftools qw{ bcftools_view };
    use MIP::Program::Variantcalling::Tiddit qw{ tiddit_coverage };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $program_mode = $active_parameter_href->{$program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$program_name}{chain};
    my $program_outdirectory_name =
      $parameter_href->{$program_name}{outdir_name};
    my ( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name      => $program_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $family_id,
            FILEHANDLE            => $FILEHANDLE,
            job_id_href           => $job_id_href,
            process_time          => $time,
            program_directory =>
              catfile( $outaligner_dir, $program_outdirectory_name ),
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ## Tags
    my $infile_tag;
    my $outfile_tag;

    ## Files
    my $infile_prefix;
    my $sample_outfile_prefix;
    my $merged_sv_vcf;

    ## Paths
    my %file_path_prefix;
    my $merged_sv_vcf_path;

    # Copy family-merged SV VCF file in temporary directory:
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir );
    $infile_tag =
      $file_info_href->{$family_id}{sv_combinevariantcallsets}{file_tag};

    $merged_sv_vcf = $family_id . $infile_tag . q{SV} . $DOT . q{vcf};
    $merged_sv_vcf_path = catfile( $infamily_directory, $merged_sv_vcf );

    say {$FILEHANDLE} q{## Log vcf2cytosure version - use dummy parameters}
      . $NEWLINE;
    my $stderrfile_path = $program_info_path . $DOT . q{stderr.txt};
    vcf2cytosure_convert(
        {
            coverage_file   => q{Na},
            FILEHANDLE      => $FILEHANDLE,
            stderrfile_path => $stderrfile_path,
            vcf_infile_path => q{Na},
            version         => 1,
        }
    );

    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE}
      q{## Copy family-level merged SV VCF file to temporary directory}
      . $NEWLINE;

    migrate_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $merged_sv_vcf_path,
            outfile_path => $temp_directory,
        }
    );

    my $process_batches_count = 1;

    # Loop over all samples
    while ( my ( $sample_id_index, $sample_id ) =
        each @{ $active_parameter_href->{sample_ids} } )
    {

        say {$FILEHANDLE} q{## Processing sample} . $SPACE . $sample_id;

     # Using tiddit coverage, create coverage file from .bam file of this sample
        my $insample_directory = catdir( $active_parameter_href->{outdata_dir},
            $sample_id, $outaligner_dir );

        ## Add merged infile name prefix after merging all BAM files per sample_id
        my $merged_infile_prefix = get_merged_infile_prefix(
            {
                file_info_href => $file_info_href,
                sample_id      => $sample_id,
            }
        );

        ## Assign file_tags
        $infile_tag =
          $file_info_href->{$sample_id}{gatk_baserecalibration}{file_tag};
        $infile_prefix = $merged_infile_prefix . $infile_tag;
        $outfile_tag =
          $file_info_href->{$family_id}{sv_combinevariantcallsets}{file_tag};
        $sample_outfile_prefix = $merged_infile_prefix . $outfile_tag;

        ## Assign suffix
        my $infile_suffix = get_file_suffix(
            {
                jobid_chain =>
                  $parameter_href->{gatk_baserecalibration}{chain},
                parameter_href => $parameter_href,
                suffix_key     => q{alignment_file_suffix},
            }
        );

        ## Set file suffix for next module within jobid chain
        my $cov_outfile_suffix = get_file_suffix(
            {
                parameter_href => $parameter_href,
                program_name   => q{tiddit},
                suffix_key     => q{coverage_file_suffix},
            }
        );

        my $outfile_suffix = get_file_suffix(
            {
                parameter_href => $parameter_href,
                program_name   => $program_name,
                suffix_key     => q{outfile_suffix},
            }
        );

        $file_path_prefix{$sample_id}{in} =
          catfile( $temp_directory, $infile_prefix );
        $file_path_prefix{$sample_id}{out} =
          catfile( $temp_directory, $sample_outfile_prefix );

        # q{.bam} -> ".b*" for getting index as well
        my $infile_path = catfile( $insample_directory,
            $infile_prefix . substr( $infile_suffix, 0, 2 ) . $ASTERISK );

        $process_batches_count = print_wait(
            {
                FILEHANDLE            => $FILEHANDLE,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                process_counter       => $sample_id_index,
            }
        );

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy bam file to temporary directory};
        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $infile_path,
                outfile_path => $temp_directory,
            }
        );
        say {$FILEHANDLE} q{wait}, $NEWLINE;
        say {$FILEHANDLE}
          q{## Creating coverage file with tiddit -cov for sample}
          . $SPACE
          . $sample_id;

        ## Tiddit coverage
        tiddit_coverage(
            {
                bin_size    => $bin_size,
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $file_path_prefix{$sample_id}{in}
                  . $infile_suffix,
                outfile_path_prefix => $file_path_prefix{$sample_id}{out},
            }
        );
        say {$FILEHANDLE} $AMPERSAND . $SPACE . $NEWLINE;
        say {$FILEHANDLE} q{wait}, $NEWLINE;

        # Extract SV from this sample from merged SV VCF file
        say {$FILEHANDLE} q{## Using bcftools_view to extract SVs for sample}
          . $SPACE
          . $sample_id
          . $NEWLINE;

        $infile_tag =
          $file_info_href->{$sample_id}{sv_combinevariantcallsets}{file_tag};
        my $sample_vcf_file = $sample_id . $infile_tag . q{SV} . $DOT . q{vcf};

        # Bcftools view
        bcftools_view(
            {
                exclude =>
                  $active_parameter_href->{vcf2cytosure_exclude_filter},
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => catfile( $temp_directory, $merged_sv_vcf ),
                samples_ref  => [$sample_id],
                outfile_path => catfile( $temp_directory, $sample_vcf_file ),
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        say {$FILEHANDLE}
q{## Converting sample's SV VCF file into cytosure, using Vcf2cytosure}
          . $NEWLINE;

        my $cgh_outfile_path = catfile( $temp_directory,
            $sample_id . $infile_tag . q{SV} . $DOT . q{cgh} );

        vcf2cytosure_convert(
            {
                coverage_file => $file_path_prefix{$sample_id}{out}
                  . $cov_outfile_suffix,
                FILEHANDLE      => $FILEHANDLE,
                outfile_path    => $cgh_outfile_path,
                vcf_infile_path => catfile( $temp_directory, $sample_vcf_file ),
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $cgh_outfile_path,
                outfile_path => $outfamily_directory,
            }
        );
        say {$FILEHANDLE} q{wait}, $NEWLINE;

        if ( $program_mode == 1 ) {

            add_program_outfile_to_sample_info(
                {
                    infile    => $merged_infile_prefix,
                    sample_id => $sample_id,
                    path      => catfile(
                        $outfamily_directory,
                        $sample_id . $infile_tag . q{SV} . $DOT . q{cgh}
                    ),
                    program_name     => q{vcf2cytosure},
                    sample_info_href => $sample_info_href,
                }
            );

            ## For logging version - until present in cgh file
            add_program_outfile_to_sample_info(
                {
                    infile           => $merged_infile_prefix,
                    sample_id        => $sample_id,
                    path             => $stderrfile_path,
                    program_name     => q{vcf2cytosure_version},
                    sample_info_href => $sample_info_href,
                }
            );

        }
    }

    if ( $program_mode == 1 ) {

        slurm_submit_job_sample_id_dependency_family_dead_end(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                sbatch_file_name => $file_path,
            }
        );
    }

    return;

}

1;
