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
    our @EXPORT_OK = qw{ analysis_vardict_tn };

}

## Constants
Readonly my $STDOUT     => q{>};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $PIPE       => q{|};
Readonly my $UNDERSCORE => q{_};
Readonly my $ASTERISK   => q{*};

sub analysis_vardict_tn {

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
##          : $sample_origin           => Info on sample origin (tumor vs normal)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $outfamily_directory;
    my $outsample_directory;
    my $parameter_href;
    my $program_name;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory,
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Cluster qw{ check_max_core_number };
    use MIP::Cluster qw{ get_core_number };
    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files
      qw{ migrate_file migrate_files xargs_migrate_contig_files };
    use MIP::Program::Variantcalling::Vardict
      qw{ vardict vardict_var2vcf_single vardict_var2vcf_paired vardict_testsomatic vardict_teststrandbias };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::QC::Record
      qw{ add_program_metafile_to_sample_info add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw { set_file_suffix };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    ## Add get core number not more than max
    $core_number = get_core_number(
        {
            modifier_core_number =>
              scalar( @{ $active_parameter_href->{sample_ids} } ),
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $family_id,
            FILEHANDLE            => $FILEHANDLE,
            job_id_href           => $job_id_href,
            process_time          => $time,
            program_directory     => catfile( $outaligner_dir, q{vardict} ),
            program_name          => $program_name,
            source_environment_commands_ref => [$source_environment_cmd],
            temp_directory                  => $temp_directory,
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

    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            file_suffix => $parameter_href->{$mip_program_name}{outfile_suffix},
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{outfile_suffix},
        }
    );

    ## Assign file_tags
    my $outfile_tag =
      $file_info_href->{$family_id}{$mip_program_name}{file_tag};

    my $outfile_prefix = $family_id . $outfile_tag;

    ## Files
    my $outfile_name = $outfile_prefix . $outfile_suffix;

    ## Paths
    my $outfile_path = catfile( $outsample_directory, $outfile_name );

    my %file_path_prefix;

  SAMPLE:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Assign directories
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
        my $infile_tag =
          $file_info_href->{$sample_id}{pgatk_baserecalibration}{file_tag};

        ## Files
        my $infile_prefix = $merged_infile_prefix . $infile_tag;

        my $infile_name = $infile_prefix . $infile_suffix;

        my $infile_path = catfile( $insample_directory, $infile_name );

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $infile_path,
                outfile_path => $temp_directory,
            }
        );

        ## Paths
        $file_path_prefix{$sample_id} =
          catfile( $temp_directory, $infile_prefix );

    }

    ## Prepare an array of inputs for vardict
    my @infile_path =
      map { $file_path_prefix{$_} . $UNDERSCORE . $infile_suffix }
      @{ $active_parameter_href->{sample_ids} };

    ## Check if the array size is less than two
    ## TODO: sub in vardict module to check for size of array and prepare Tumor vs Normal set.
    ## TODO: sub in vardict module to create pairs of analysis: TvsN and NvsN.
    ## TODO: add merge bam files for more than one normal sample id. This should be added to picardtools_mergesamfiles as a sub.

    ## Vardict
    ## Note: Vardict doesn't runs the analysis based on the input bed file, so if the input bed file points towards the
    ## whole genome, Vardict will run it on whole genome, but mind you, this might take a long time (up to 3 weeks). So
    ## in order to reduce this time two measures need to be taken: 1) use vardict-java to support multiple threads, 2)
    ## run bedtools genomecov and run it on callable regions the same way GATK is doing. After having the callable
    ## regions, by having multiple bed files, it can create a more proper form of parallelization. For now vardict's
    ## recipe is designed to accommodate small bed files with small regions (it might be good to have the regions
    ## overlap with each other).

    ## Create input filename bam file for $infile_paths_ref in vardict module. This has to be an array where first
    ## is tumor and the second one is normal tissue.
    ## Get parameters
    # Unpack
    my $af_threshold       = $active_parameter_href->{vrd_af_threshold};
    my $chrom_start        = $active_parameter_href->{vrd_chrom_start};
    my $input_bed_file     = $active_parameter_href->{vrd_input_bed_file};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $region_end         = $active_parameter_href->{vrd_region_end};
    my $region_start       = $active_parameter_href->{vrd_region_start};

    ## Assign family_id to sample_name
    my $sample_name    = $active_parameter_href->{family_id};
    my $segment_annotn = $active_parameter_href->{vrd_segment_annotn};

    say {$FILEHANDLE} q{## varDict variant calling};

    vardict(
        {
            af_threshold           => $af_threshold,
            FILEHANDLE             => $FILEHANDLE,
            infile_paths_ref       => \@infile_path,
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
            af_threshold    => $af_threshold,
            FILEHANDLE      => $FILEHANDLE,
            sample_name     => $sample_name,
            stdoutfile_path => $outfile_path,
        }
    );

    print {$FILEHANDLE} $STDOUT . print {$outfile_path} . $NEWLINE

      ## Copies file from temporary directory.
      say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $outfile_prefix . $outfile_suffix . $ASTERISK,
            outfile_path => $outfamily_directory,
        }
    );
    say {$FILEHANDLE} q{wait} . $NEWLINE;

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});
      or $log->logcroak(q{Could not close XARGSFILEHANDLE});

    if ( $mip_program_mode == 1 ) {

        add_program_outfile_to_sample_info(
            {
                path => catfile(
                    $outfamily_directory, $outfile_prefix . $outfile_suffix
                ),
                program_name     => q{vardict},
                sample_info_href => $sample_info_href,
            }
        );

        slurm_submit_job_sample_id_dependency_add_to_family(
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
