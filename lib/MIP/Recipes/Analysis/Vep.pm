package MIP::Recipes::Analysis::Vep;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ fileparse };
use File::Spec::Functions qw{ catdir catfile splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use POSIX;
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
    our $VERSION = 1.07;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_vep analysis_vep_rio analysis_vep_sv };

}

## Constants
Readonly my $ASTERIX                => q{*};
Readonly my $DOT                    => q{.};
Readonly my $EMPTY_STR              => q{};
Readonly my $NEWLINE                => qq{\n};
Readonly my $SPACE                  => q{ };
Readonly my $UNDERSCORE             => q{_};
Readonly my $ANNOTATION_DISTANCE    => 5000;
Readonly my $ANNOTATION_DISTANCE_MT => 10;

sub analysis_vep {

## Function : Varianteffectpredictor performs effect predictions and annotation of variants.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => The variant call type
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $program_info_path       => Program info path
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $file_path;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $program_info_path;
    my $sample_info_href;

    ## Default(s)
    my $call_type;
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
        call_type =>
          { default => q{BOTH}, store => \$call_type, strict_type => 1, },
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
        file_path               => { store => \$file_path, strict_type => 1, },
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
        program_info_path =>
          { store => \$program_info_path, strict_type => 1, },
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
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw{ get_core_number };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Vep qw{ variant_effect_predictor };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };
    use MIP::QC::Record
      qw{ add_program_metafile_to_sample_info add_program_outfile_to_sample_info };

    ## Constants
    Readonly my $VEP_FORK_NUMBER => 4;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set program mode
    my $program_mode = $active_parameter_href->{$program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$program_name}{chain};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );
    my $xargs_file_path_prefix;

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
            modifier_core_number => scalar @{ $file_info_href->{contigs} },
            module_core_number =>
              $active_parameter_href->{module_core_number}{$program_name},
        }
    );

    # Adjust for the number of forks vep forks
    $core_number = floor( $core_number / $VEP_FORK_NUMBER );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            call_type                       => $call_type,
            core_number                     => $core_number,
            directory_id                    => $family_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            program_directory               => $outaligner_dir,
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );
    my $stderr_path = $program_info_path . $DOT . q{stderr.txt};

    # Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) = splitpath($stderr_path);

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir );
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir );

    ## Used downstream
    $parameter_href->{$program_name}{indirectory} = $outfamily_directory;

    ## Assign file_tags
    my $infile_tag = $file_info_href->{$family_id}{frequency_filter}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$program_name}{file_tag};
    my $infile_prefix       = $family_id . $infile_tag . $call_type;
    my $outfile_prefix      = $family_id . $outfile_tag . $call_type;
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain
    my $infile_suffix = get_file_suffix(
        {
            jobid_chain    => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );
    my $outfile_suffix = set_file_suffix(
        {
            file_suffix    => $parameter_href->{$program_name}{outfile_suffix},
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    ($xargs_file_counter) = xargs_migrate_contig_files(
        {
            contigs_ref        => \@{ $file_info_href->{contigs_size_ordered} },
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $file_path,
            indirectory        => $infamily_directory,
            infile             => $infile_prefix,
            program_info_path  => $program_info_path,
            temp_directory     => $temp_directory,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Varianteffectpredictor
    say {$FILEHANDLE} q{## Varianteffectpredictor};

    my $assembly_version = $file_info_href->{human_genome_reference_source}
      . $file_info_href->{human_genome_reference_version};

    ## Get genome source and version to be compatible with VEP
    $assembly_version =
      _get_assembly_name( { assembly_version => $assembly_version } );

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
        my $distance = $ANNOTATION_DISTANCE;

        # Special case for mitochondrial contig annotation
        if ( $contig =~ / MT|M /xsm ) {

            $distance = $ANNOTATION_DISTANCE_MT;
        }

        # VEP plugins
        my @plugins;
      PLUGIN:
        foreach my $plugin ( @{ $active_parameter_href->{vep_plugins} } ) {

            if ( $plugin eq q{LoF} ) {

                my $lof_parameter = q{,human_ancestor_fa:}
                  . catfile( $active_parameter_href->{vep_directory_cache},
                    q{Plugins}, q{human_ancestor.fa,filter_position:0.05} );
                push @plugins, $plugin . $lof_parameter;
            }
            elsif ( $plugin eq q{MaxEntScan} ) {

                my $max_ent_scan_parameter = q{,}
                  . catfile( $active_parameter_href->{vep_directory_cache},
                    qw{ Plugins fordownload } );
                push @plugins, $plugin . $max_ent_scan_parameter;
            }
            else {

                push @plugins, $plugin;
            }
        }

        ## VEP features
        my @vep_features_ref;
      FEATURE:
        foreach my $vep_feature ( @{ $active_parameter_href->{vep_features} } )
        {

            # Add VEP features to the output.
            push @vep_features_ref, $vep_feature;

            # Special case for mitochondrial contig annotation
            if ( $contig =~ / MT|M /xsm && $vep_feature eq q{refseq} ) {

                push @vep_features_ref, q{all_refseq};
            }
        }

        my $infile_path =
          $file_path_prefix . $UNDERSCORE . $contig . $infile_suffix;
        my $outfile_path =
          $outfile_path_prefix . $UNDERSCORE . $contig . $infile_suffix;
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        my $stdoutfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stdout.txt};
        variant_effect_predictor(
            {
                assembly    => $assembly_version,
                buffer_size => 20_000,
                cache_directory =>
                  $active_parameter_href->{vep_directory_cache},
                distance       => $distance,
                FILEHANDLE     => $XARGSFILEHANDLE,
                fork           => $VEP_FORK_NUMBER,
                infile_format  => substr( $outfile_suffix, 1 ),
                infile_path    => $infile_path,
                outfile_format => substr( $outfile_suffix, 1 ),
                outfile_path   => $outfile_path,
                plugins_dir_path =>
                  $active_parameter_href->{vep_plugins_dir_path},
                plugins_ref => \@plugins,
                reference_path =>
                  $active_parameter_href->{human_genome_reference},
                regions_ref      => [$contig],
                stderrfile_path  => $stderrfile_path,
                stdoutfile_path  => $stdoutfile_path,
                vep_features_ref => \@vep_features_ref,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## QC Data File(s)
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $outfile_path_prefix
              . $UNDERSCORE
              . $ASTERIX
              . $infile_suffix
              . $UNDERSCORE . q{s}
              . $ASTERIX,
            outfile_path => $outfamily_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $XARGSFILEHANDLE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $outfile_path_prefix
              . $UNDERSCORE
              . $ASTERIX
              . $infile_suffix
              . $ASTERIX,
            outfile_path => $outfamily_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        ## Collect QC metadata info for later use
        my $qc_vep_summary_outfile_path = catfile( $outfamily_directory,
                $outfile_prefix
              . $UNDERSCORE
              . $file_info_href->{contigs_size_ordered}[0]
              . $infile_suffix
              . $UNDERSCORE
              . q{summary.html} );
        add_program_metafile_to_sample_info(
            {
                metafile_tag     => q{summary},
                path             => $qc_vep_summary_outfile_path,
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );
        add_program_metafile_to_sample_info(
            {
                metafile_tag     => q{stderrfile},
                path             => catfile( $directory, $stderr_file ),
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );

        ## Collect QC metadata info for later use
        my $qc_vep_outfile_path = catfile( $outfamily_directory,
                $outfile_prefix
              . $UNDERSCORE
              . $file_info_href->{contigs_size_ordered}[0]
              . $infile_suffix );
        add_program_outfile_to_sample_info(
            {
                path             => $qc_vep_outfile_path,
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                log                     => $log,
                job_id_href             => $job_id_href,
                path                    => $job_id_chain,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                sbatch_file_name => $file_path,
            }
        );
    }
    return;
}

sub analysis_vep_rio {

## Function : Varianteffectpredictor performs effect predictions and annotation of variants.
## Returns  : undef | $xargs_file_counter
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => The variant call type
##          : $family_id               => Family id
##          : $FILEHANDLE              => Filehandle to write to
##          : $file_info_href          => The file_info hash {REF}
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $program_info_path       => The program info path
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $stderr_path             => Stderr path of the block script
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $FILEHANDLE;
    my $file_info_href;
    my $file_path;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $program_info_path;
    my $sample_info_href;
    my $stderr_path;

    ## Default(s)
    my $call_type;
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
        call_type =>
          { default => q{BOTH}, store => \$call_type, strict_type => 1, },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
            strict_type => 1,
        },
        FILEHANDLE     => { required => 1, store => \$FILEHANDLE, },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_path               => { store => \$file_path, strict_type => 1, },
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
        program_info_path =>
          { store => \$program_info_path, strict_type => 1, },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        stderr_path => {
            defined     => 1,
            required    => 1,
            store       => \$stderr_path,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
        xargs_file_counter => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw{ get_core_number };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Vep qw{ variant_effect_predictor };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Set::File qw{ set_file_suffix };
    use MIP::Script::Setup_script
      qw{ write_return_to_conda_environment write_source_environment_command };
    use MIP::QC::Record
      qw{ add_program_metafile_to_sample_info add_program_outfile_to_sample_info };

    ## Constants
    Readonly my $VEP_FORK_NUMBER => 4;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set program mode
    my $program_mode = $active_parameter_href->{$program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$program_name}{chain};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    my $reduce_io = $active_parameter_href->{reduce_io};
    my $xargs_file_path_prefix;

    ## Filehandles
    # Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
            modifier_core_number => scalar @{ $file_info_href->{contigs} },
            module_core_number =>
              $active_parameter_href->{module_core_number}{$program_name},
        }
    );

    # Adjust for the number of forks vep forks
    $core_number = floor( $core_number / $VEP_FORK_NUMBER );

    ## If program needs special environment variables set
    if (@source_environment_cmds) {

        write_source_environment_command(
            {
                FILEHANDLE                      => $FILEHANDLE,
                source_environment_commands_ref => \@source_environment_cmds,
            }
        );
    }

    # Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) = splitpath($stderr_path);

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir );
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir );

    ## Used downstream
    $parameter_href->{$program_name}{indirectory} = $outfamily_directory;

    ## Assign file_tags
    my $infile_tag = $file_info_href->{$family_id}{frequency_filter}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$program_name}{file_tag};
    my $infile_prefix       = $family_id . $infile_tag . $call_type;
    my $outfile_prefix      = $family_id . $outfile_tag . $call_type;
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain
    my $infile_suffix = get_file_suffix(
        {
            jobid_chain    => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );
    my $outfile_suffix = set_file_suffix(
        {
            file_suffix    => $parameter_href->{$program_name}{outfile_suffix},
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

    ## Varianteffectpredictor
    say {$FILEHANDLE} q{## Varianteffectpredictor};

    my $assembly_version = $file_info_href->{human_genome_reference_source}
      . $file_info_href->{human_genome_reference_version};

    ## Get genome source and version to be compatible with VEP
    $assembly_version =
      _get_assembly_name( { assembly_version => $assembly_version } );

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
        my $distance = $ANNOTATION_DISTANCE;

        # Special case for mitochondrial contig annotation
        if ( $contig =~ / MT|M /xsm ) {

            $distance = $ANNOTATION_DISTANCE_MT;
        }

        # VEP plugins
        my @plugins;
      PLUGIN:
        foreach my $plugin ( @{ $active_parameter_href->{vep_plugins} } ) {

            if ( $plugin eq q{LoF} ) {

                my $lof_parameter = q{,human_ancestor_fa:}
                  . catfile( $active_parameter_href->{vep_directory_cache},
                    q{Plugins}, q{human_ancestor.fa,filter_position:0.05} );
                push @plugins, $plugin . $lof_parameter;
            }
            elsif ( $plugin eq q{MaxEntScan} ) {

                my $lof_parameter = q{,}
                  . catfile( $active_parameter_href->{vep_directory_cache},
                    qw{ Plugins fordownload } );
                push @plugins, $plugin . $lof_parameter;
            }
            else {

                push @plugins, $plugin;
            }
        }

        ## VEP features
        my @vep_features_ref;
      FEATURE:
        foreach my $vep_feature ( @{ $active_parameter_href->{vep_features} } )
        {

            # Add VEP features to the output.
            push @vep_features_ref, $vep_feature;

            # Special case for mitochondrial contig annotation
            if ( ( $contig =~ /MT|M/xsm ) && ( $vep_feature eq q{refseq} ) ) {

                push @vep_features_ref, q{all_refseq};
            }
        }

        my $infile_path =
          $file_path_prefix . $UNDERSCORE . $contig . $infile_suffix;
        my $outfile_path =
          $outfile_path_prefix . $UNDERSCORE . $contig . $infile_suffix;
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        my $stdoutfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stdout.txt};
        variant_effect_predictor(
            {
                assembly => $assembly_version,
                cache_directory =>
                  $active_parameter_href->{vep_directory_cache},
                buffer_size    => 20_000,
                distance       => $distance,
                FILEHANDLE     => $XARGSFILEHANDLE,
                fork           => $VEP_FORK_NUMBER,
                infile_format  => substr( $outfile_suffix, 1 ),
                infile_path    => $infile_path,
                outfile_format => substr( $outfile_suffix, 1 ),
                outfile_path   => $outfile_path,
                plugins_dir_path =>
                  $active_parameter_href->{vep_plugins_dir_path},
                plugins_ref => \@plugins,
                reference_path =>
                  $active_parameter_href->{human_genome_reference},
                regions_ref      => [$contig],
                stderrfile_path  => $stderrfile_path,
                stdoutfile_path  => $stdoutfile_path,
                vep_features_ref => \@vep_features_ref,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    close $XARGSFILEHANDLE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};

    # QC Data File(s)
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $outfile_path_prefix
              . $UNDERSCORE
              . $ASTERIX
              . $infile_suffix
              . $UNDERSCORE . q{s}
              . $ASTERIX,
            outfile_path => $outfamily_directory,
        }
    );

    # Move file for downstream collection of VEP version

    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $outfile_path_prefix
              . $UNDERSCORE
              . $file_info_href->{contigs_size_ordered}[0]
              . $infile_suffix,
            outfile_path => $outfamily_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Return to main or default environment using conda
    write_return_to_conda_environment(
        {
            FILEHANDLE => $FILEHANDLE,
            source_main_environment_commands_ref =>
              \@{ $active_parameter_href->{source_main_environment_commands} },
        }
    );

    if ( $program_mode == 1 ) {

        ## Collect QC metadata info for later use
        my $qc_vep_summary_outfile_path = catfile( $outfamily_directory,
                $outfile_prefix
              . $UNDERSCORE
              . $file_info_href->{contigs_size_ordered}[0]
              . $infile_suffix
              . $UNDERSCORE
              . q{summary.html} );
        add_program_metafile_to_sample_info(
            {
                metafile_tag     => q{summary},
                path             => $qc_vep_summary_outfile_path,
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );

        add_program_metafile_to_sample_info(
            {
                metafile_tag     => q{stderrfile},
                path             => catfile( $directory, $stderr_file ),
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );
        my $qc_vep_outfile_path = catfile( $outfamily_directory,
                $outfile_prefix
              . $UNDERSCORE
              . $file_info_href->{contigs_size_ordered}[0]
              . $infile_suffix );
        add_program_outfile_to_sample_info(
            {
                path             => $qc_vep_outfile_path,
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );
    }

    # Track the number of created xargs scripts per module for Block algorithm
    return $xargs_file_counter;
}

sub analysis_vep_sv {

## Function : Varianteffectpredictor performs annotation of SV variants.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => The variant call type
##          : $family_id               => Family id
##          : $FILEHANDLE              => Filehandle to write to
##          : $file_info_href          => The file_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $program_info_path       => The program info path
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
    my $sample_info_href;

    ## Default(s)
    my $call_type;
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
        call_type => {
            allow       => [qw{ SV }],
            default     => q{SV},
            store       => \$call_type,
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
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir,
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
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw{ get_core_number };
    use MIP::Delete::List qw{ delete_contig_elements delete_male_contig };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };

    use MIP::Program::Variantcalling::Vep qw{ variant_effect_predictor };

    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::QC::Record
      qw{ add_program_metafile_to_sample_info add_program_outfile_to_sample_info };

    ## Constants
    Readonly my $VEP_FORK_NUMBER => 4;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program mode
    my $program_mode = $active_parameter_href->{$program_name};

    ## Alias
    my $contigs_ref = \@{ $file_info_href->{contigs} };
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $job_id_chain = $parameter_href->{$program_name}{chain};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );
    my $xargs_file_path_prefix;

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
            modifier_core_number => scalar @{$contigs_ref},
            module_core_number =>
              $active_parameter_href->{module_core_number}{$program_name},
        }
    );

    # Adjust for the number of forks vep forks
    $core_number = floor( $core_number / $VEP_FORK_NUMBER );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            call_type                       => $call_type,
            core_number                     => $core_number,
            directory_id                    => $family_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            program_directory               => catfile( lc $outaligner_dir ),
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );
    my $stderr_path = $program_info_path . $DOT . q{stderr.txt};

    # Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) = splitpath($stderr_path);

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir );
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir );

    # Used downstream
    $parameter_href->{$program_name}{indirectory} = $outfamily_directory;

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$family_id}{sv_annotate}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$program_name}{file_tag};
    my $infile_prefix       = $family_id . $infile_tag . $call_type;
    my $outfile_prefix      = $family_id . $outfile_tag . $call_type;
    my $infile_path_prefix  = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ## Assign suffix
    my $file_suffix = get_file_suffix(
        {
            jobid_chain    => $job_id_chain,
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
                $infamily_directory, $infile_prefix . $file_suffix . $ASTERIX
            ),
            outfile_path => $temp_directory
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Reformat SV with no length as these will fail in the annotation with VEP
    _reformat_sv_with_no_length(
        {
            FILEHANDLE         => $FILEHANDLE,
            file_suffix        => $file_suffix,
            infile_path_prefix => $infile_path_prefix,
        }
    );

    ## varianteffectpredictor
    say {$FILEHANDLE} q{## Varianteffectpredictor};

    my $assembly_version = $file_info_href->{human_genome_reference_source}
      . $file_info_href->{human_genome_reference_version};

    ## Get genome source and version to be compatible with VEP
    $assembly_version =
      _get_assembly_name( { assembly_version => $assembly_version } );

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
    foreach my $contig ( @{$contigs_ref} ) {

        ## Get parameters
        my $distance = $ANNOTATION_DISTANCE;

        # Special case for mitochondrial contig annotation
        if ( $contig =~ / MT|M /xsm ) {

            $distance = $ANNOTATION_DISTANCE_MT;
        }

        my $vep_outfile_prefix         = $outfile_prefix;
        my $vep_xargs_file_path_prefix = $xargs_file_path_prefix;
        my @regions;
        my $mt_name;

        ## Contig specific
        # Update endings with contig info
        if (   $consensus_analysis_type eq q{wgs}
            || $consensus_analysis_type eq q{mixed}
            || $consensus_analysis_type eq q{vrn} )
        {

            $vep_outfile_prefix = $outfile_prefix . $UNDERSCORE . $contig;
            $vep_xargs_file_path_prefix =
              $xargs_file_path_prefix . $DOT . $contig;
            push @regions, $contig;

            ## The MT needs to be analyzed together with another contig to guarantee that SV:s are found.
            ## Otherwise VEP will complain.
            if ( $contig =~ /chrM/xms ) {
                unshift @regions, q{chr21};
                $mt_name = $contig;
            }
            if ( $contig =~ /MT/xms ) {
                unshift @regions, q{21};
                $mt_name = $contig;
            }
        }

        ## VEP plugins
        my @plugins;

      PLUGIN:
        foreach my $plugin ( @{ $active_parameter_href->{sv_vep_plugins} } ) {

            if ( $plugin eq q{LoF} ) {

                my $lof_parameter = q{,human_ancestor_fa:}
                  . catfile(
                    $active_parameter_href->{vep_directory_cache},
                    q{human_ancestor.fa,filter_position:0.05}
                  );
                push @plugins, $plugin . $lof_parameter;
            }
            else {

                push @plugins, $plugin;
            }
        }

        ## VEPFeatures
        my @vep_features_ref;

      FEATURE:
        foreach
          my $vep_feature ( @{ $active_parameter_href->{sv_vep_features} } )
        {

            # Add VEP features to the output.
            push @vep_features_ref, $vep_feature;

            # Special case for mitochondrial contig annotation
            if ( ( $contig =~ /MT|M/sxm ) && ( $vep_feature eq q{refseq} ) ) {

                push @vep_features_ref, q{all_refseq};
            }
        }

        my $infile_path =
          $infile_path_prefix . $UNDERSCORE . q{fixedsvlength} . $file_suffix;
        my $outfile_path =
          catfile( $temp_directory, $vep_outfile_prefix . $file_suffix );
        my $stderrfile_path =
          $vep_xargs_file_path_prefix . $DOT . q{stderr.txt};
        my $stdoutfile_path =
          $vep_xargs_file_path_prefix . $DOT . q{stdout.txt};
        variant_effect_predictor(
            {
                assembly    => $assembly_version,
                buffer_size => 100,
                cache_directory =>
                  $active_parameter_href->{vep_directory_cache},
                distance       => $distance,
                FILEHANDLE     => $XARGSFILEHANDLE,
                fork           => $VEP_FORK_NUMBER,
                infile_format  => substr( $file_suffix, 1 ),
                infile_path    => $infile_path,
                outfile_format => substr( $file_suffix, 1 ),
                outfile_path   => $outfile_path,
                plugins_dir_path =>
                  $active_parameter_href->{vep_plugins_dir_path},
                plugins_ref => \@plugins,
                regions_ref => \@regions,
                reference_path =>
                  $active_parameter_href->{human_genome_reference},
                stderrfile_path  => $stderrfile_path,
                stdoutfile_path  => $stdoutfile_path,
                vep_features_ref => \@vep_features_ref,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;

# Only perform once for exome samples to avoid risking contigs lacking variants throwing errors
        if ( $consensus_analysis_type eq q{wes} ) {

            last CONTIG;
        }

        ## Filter out the MT annotations from the combined chr21 chrM call
        if ($mt_name) {

            say {$FILEHANDLE} q{## Filter out MT annotations};
            _subset_vcf(
                {
                    FILEHANDLE   => $FILEHANDLE,
                    infile_path  => $outfile_path,
                    outfile_path => $outfile_path,
                    regions_ref  => [qw{ chrM MT }],
                }
            );
        }
    }

    close $XARGSFILEHANDLE;

    ## QC Data File(s)
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $outfile_path_prefix
              . $ASTERIX
              . $file_suffix
              . $UNDERSCORE . q{s}
              . $ASTERIX,
            outfile_path => $outfamily_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $outfile_path_prefix
              . $ASTERIX
              . $file_suffix
              . $ASTERIX,
            outfile_path => $outfamily_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        my $outfile_sample_info_prefix = $outfile_prefix;

        # Update endings with contig info
        if (   $consensus_analysis_type eq q{wgs}
            || $consensus_analysis_type eq q{mixed}
            || $consensus_analysis_type eq q{vrn} )
        {

            $outfile_sample_info_prefix .= $UNDERSCORE . $contigs_ref->[0];
        }

        ## Collect QC metadata info for later use
        my $qc_vep_summary_outfile =
          $outfile_sample_info_prefix . $DOT . q{vcf_summary.html};
        add_program_metafile_to_sample_info(
            {
                directory        => $outfamily_directory,
                file             => $qc_vep_summary_outfile,
                metafile_tag     => q{summary},
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                path => catfile(
                    $outfamily_directory,
                    $outfile_sample_info_prefix . $file_suffix
                ),
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                log                     => $log,
                job_id_href             => $job_id_href,
                path                    => $job_id_chain,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                sbatch_file_name => $file_path,
            }
        );
    }
    return;
}

sub _get_assembly_name {

## Function : Get genome source and version to be compatible with VEP
## Returns  : $assembly_version
## Arguments: $assembly_version => The genome source and version to be checked

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $assembly_version;

    my $tmpl = {
        assembly_version => {
            defined     => 1,
            required    => 1,
            store       => \$assembly_version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    Readonly my $VEP_CACHE_GENOME_BUILD_NAME_SWITCH => 20;

    if ( $assembly_version =~ /hg(\d+)/xsm ) {

        my $version_number = $1;

        if ( $version_number > $VEP_CACHE_GENOME_BUILD_NAME_SWITCH ) {

            $assembly_version = q{GRCh} . $version_number;
        }
    }
    return $assembly_version;
}

sub _reformat_sv_with_no_length {

## Function : Reformat SV with no length as these will fail in the annotation with VEP
## Returns  :
## Arguments: $infile_path_prefix => Infile path prefix
##          : $file_suffix        => File suffix
##          : $FILEHANDLE         => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path_prefix;
    my $file_suffix;
    my $FILEHANDLE;

    my $tmpl = {
        infile_path_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path_prefix,
            strict_type => 1,
        },
        file_suffix => {
            defined     => 1,
            required    => 1,
            store       => \$file_suffix,
            strict_type => 1,
        },
        FILEHANDLE => { required => 1, store => \$FILEHANDLE, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Fix SV with no length as these will fail in the annotation with VEP
    my $perl_fix_sv_nolengths;

    # Set up perl
    $perl_fix_sv_nolengths .= q?perl -nae '?;

    # Initate variables
    $perl_fix_sv_nolengths .=
      q?my %info; my $start; my $end; my $alt; my @data; ?;

    # Split line
    $perl_fix_sv_nolengths .= q?@data=split("\t", $_); ?;

    # Add $start position
    $perl_fix_sv_nolengths .= q?$start = $data[1]; $start++; ?;

    # Add $alt allele
    $perl_fix_sv_nolengths .= q?$alt=$data[4]; ?;

    # Add INFO field to %data using $key->$value
    $perl_fix_sv_nolengths .=
q?foreach my $bit (split /\;/, $data[7]) { my ($key, $value) = split /\=/, $bit; $info{$key} = $value; } ?;

    # Add $end position
    $perl_fix_sv_nolengths .=
      q?if(defined($info{END})) { $end = $info{END}; } ?;

# If SV, strip SV type entry and check if no length, then do not print variant else print
    $perl_fix_sv_nolengths .=
q?if($alt=~ /\<|\[|\]|\>/) { $alt=~ s/\<|\>//g; $alt=~ s/\:.+//g; if($start >= $end && $alt=~ /del/i) {} else {print $_} } ?;

    # All other lines - print
    $perl_fix_sv_nolengths .= q?else {print $_}' ?;

    print {$FILEHANDLE} $perl_fix_sv_nolengths . $SPACE;
    print {$FILEHANDLE} $infile_path_prefix . $file_suffix . $SPACE;
    say   {$FILEHANDLE} q{>}
      . $SPACE
      . $infile_path_prefix
      . $UNDERSCORE
      . q{fixedsvlength}
      . $file_suffix
      . $SPACE, $NEWLINE;
    return;
}

sub _subset_vcf {

## Function : Subsets the input vcf to only contain a subset, defined in the regions variable
## Returns  :
## Arguments: $FILEHANDLE   => Filehandle
##          : $infile_path  => Path to infile
##          : $outfile_path => Path to outfile
##          : $regions_ref  => Array with regions to be included in the subset {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $regions_ref;

    my $tmpl = {
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        regions_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$regions_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Utility::Htslib qw{ htslib_bgzip };
    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_index bcftools_view };

    ## Prepare for bcftools_view
    htslib_bgzip(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $infile_path,
        }
    );
    print {$FILEHANDLE} $NEWLINE;
    bcftools_index(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $infile_path . q{.gz},
        }
    );
    print {$FILEHANDLE} $NEWLINE;
    bcftools_view(
        {
            regions_ref  => $regions_ref,
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $infile_path . q{.gz},
            outfile_path => $outfile_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    return;
}

1;
