package MIP::Recipes::Analysis::Gatk_haplotypecaller;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename };
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
    our $VERSION = 1.07;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ analysis_gatk_haplotypecaller analysis_gatk_haplotypecaller_rna };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_gatk_haplotypecaller {

## Function : Gatk haplotypecaller analysis recipe
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
##          : $xargs_file_counter      => The xargs file counter

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
    my $xargs_file_counter;

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
        xargs_file_counter => {
            allow       => qr{ ^\d+$ }xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Pedigree qw{ create_fam_file };
    use MIP::File::Interval qw{ generate_contig_interval_file };
    use MIP::Get::File qw{ get_exom_target_bed_file get_io_files };
    use MIP::Get::Parameter qw{ get_module_parameters get_program_attributes };
    use MIP::IO::Files qw{ xargs_migrate_contig_files };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::Program::Alignment::Gatk qw{ gatk_haplotypecaller };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $JAVA_MEMORY_ALLOCATION        => 8;
    Readonly my $STANDARD_MIN_CONFIDENCE_THRSD => 10;

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            program_name   => $program_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my $indir_path_prefix       = $io{in}{dir_path_prefix};
    my $infile_suffix           = $io{in}{file_suffix};
    my $infile_name_prefix      = $io{in}{file_name_prefix};
    my $temp_infile_name_prefix = $io{temp}{file_name_prefix};
    my %temp_infile_path        = %{ $io{temp}{file_path_href} };

    my $job_id_chain = get_program_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            program_name   => $program_name,
        }
    );
    my $program_mode       = $active_parameter_href->{$program_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $analysis_type = $active_parameter_href->{analysis_type}{$sample_id};
    my $xargs_file_path_prefix;

    ## Get module parameters
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
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

    ## Outpaths
    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $sample_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => $file_info_href->{contigs_size_ordered},
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                program_name     => $program_name,
                temp_directory   => $temp_directory,
            }
        )
    );
    my @outfile_paths       = @{ $io{out}{file_paths} };
    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_name_prefix = $io{out}{file_name_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my %temp_outfile_path   = %{ $io{temp}{file_path_href} };

    ## For wes
    my $file_ending_pad_interval_list = $file_info_href->{exome_target_bed}[2];
    my $exome_target_bed_file         = get_exom_target_bed_file(
        {
            exome_target_bed_href => $active_parameter_href->{exome_target_bed},
            file_ending           => $file_ending_pad_interval_list,
            log                   => $log,
            sample_id             => $sample_id,
        }
    );

    ### SHELL:

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $sample_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            program_directory               => $program_name,
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    # For ".fam" file
    my $outfamily_file_directory =
      catdir( $active_parameter_href->{outdata_dir}, $family_id );

    ## Create .fam file to be used in variant calling analyses
    my $fam_file_path =
      catfile( $outfamily_file_directory, $family_id . $DOT . q{fam} );
    create_fam_file(
        {
            active_parameter_href => $active_parameter_href,
            fam_file_path         => $fam_file_path,
            FILEHANDLE            => $FILEHANDLE,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Get exome_target_bed file for specfic sample_id and add file_ending from file_info hash if supplied

    ## Exome analysis
    if ( $analysis_type eq q{wes} ) {

        ## Generate contig specific interval_list
        generate_contig_interval_file(
            {
                contigs_ref => \@{ $file_info_href->{contigs_size_ordered} },
                exome_target_bed_file => $exome_target_bed_file,
                FILEHANDLE            => $FILEHANDLE,
                max_cores_per_node    => $core_number,
                outdirectory          => $temp_directory,
                reference_dir => $active_parameter_href->{reference_dir},
            }
        );

        ## Reroute to only filename
        $exome_target_bed_file = basename($exome_target_bed_file);
    }

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    ($xargs_file_counter) = xargs_migrate_contig_files(
        {
            contigs_ref        => \@{ $file_info_href->{contigs_size_ordered} },
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_ending        => substr( $infile_suffix, 0, 2 ) . $ASTERISK,
            file_path          => $file_path,
            indirectory        => $indir_path_prefix,
            infile             => $infile_name_prefix,
            program_info_path  => $program_info_path,
            temp_directory     => $temp_directory,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## GATK HaplotypeCaller
    say {$FILEHANDLE} q{## GATK HaplotypeCaller};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## Get parameters
        my $emit_ref_confidence = q{GVCF};
        my @intervals;
        my $pcr_indel_model;
        ## Exome analysis
        if ( $analysis_type eq q{wes} ) {

            ## Limit to targets kit target file
            @intervals = (
                catfile(
                    $temp_directory,
                    $contig . $UNDERSCORE . $exome_target_bed_file
                )
            );
        }
        else {
            ## wgs

            ## Per contig
            @intervals = ($contig);

            if (
                $active_parameter_href->{gatk_haplotypecaller_pcr_indel_model} )
            {

                ## Assume that we run pcr-free sequencing (true for Rapid WGS and X-ten)
                $pcr_indel_model =
                  $active_parameter_href
                  ->{gatk_haplotypecaller_pcr_indel_model};
            }
        }
        if ( $analysis_type eq q{wts} ) {
            $emit_ref_confidence = q{NONE};
        }

        ## GATK Haplotypecaller
        my $stderrfile_path => $xargs_file_path_prefix . $DOT . $contig . $DOT
          . q{stderr.txt};
        gatk_haplotypecaller(
            {
                annotations_ref =>
                  \@{ $active_parameter_href->{gatk_haplotypecaller_annotation}
                  },
                dbsnp_path =>
                  $active_parameter_href->{gatk_haplotypecaller_snp_known_set},
                dont_use_soft_clipped_bases => $active_parameter_href
                  ->{gatk_haplotypecaller_no_soft_clipped_bases},
                emit_ref_confidence => $emit_ref_confidence,
                FILEHANDLE          => $XARGSFILEHANDLE,
                infile_path         => $temp_infile_path{$contig},
                intervals_ref       => \@intervals,
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                memory_allocation  => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
                outfile_path       => $temp_outfile_path{$contig},
                pcr_indel_model    => $pcr_indel_model,
                pedigree           => $fam_file_path,
                referencefile_path => $referencefile_path,
                standard_min_confidence_threshold_for_calling =>
                  $STANDARD_MIN_CONFIDENCE_THRSD,
                stderrfile_path => $stderrfile_path,
                temp_directory  => $temp_directory,
                verbosity       => $active_parameter_href->{gatk_logging_level},
                xargs_mode      => 1,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## Copies file from temporary directory. Per contig
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    ($xargs_file_counter) = xargs_migrate_contig_files(
        {
            core_number        => $core_number,
            contigs_ref        => \@{ $file_info_href->{contigs_size_ordered} },
            FILEHANDLE         => $FILEHANDLE,
            file_ending        => $outfile_suffix . $ASTERISK,
            file_path          => $file_path,
            outdirectory       => $outdir_path_prefix,
            outfile            => $outfile_name_prefix,
            program_info_path  => $program_info_path,
            temp_directory     => $temp_directory,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );
    close $FILEHANDLE;
    close $XARGSFILEHANDLE;

    if ( $program_mode == 1 ) {

        my $first_outfile_path = $outfile_paths[0];
        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                infile           => $outfile_name_prefix,
                path             => $first_outfile_path,
                program_name     => $program_name,
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
