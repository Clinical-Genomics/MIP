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
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_vep analysis_vep_rio analysis_vep_sv };

}

## Constants
Readonly my $ASTERIX    => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

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
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href
        },
        file_path               => { strict_type => 1, store => \$file_path },
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
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        sample_info_href  => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter
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
      qw{ add_program_outfile_to_sample_info add_program_metafile_to_sample_info};

    ## Constants
    Readonly my $VEP_FORK_NUMBER => 4;

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
    my $xargs_file_path_prefix;

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            module_core_number =>
              $active_parameter_href->{module_core_number}{$mip_program_name},
            modifier_core_number => scalar @{ $file_info_href->{contigs} },
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    # Adjust for the number of forks vep forks
    $core_number = floor( $core_number / $VEP_FORK_NUMBER );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            job_id_href                     => $job_id_href,
            FILEHANDLE                      => $FILEHANDLE,
            directory_id                    => $family_id,
            program_name                    => $program_name,
            program_directory               => $outaligner_dir,
            call_type                       => $call_type,
            core_number                     => $core_number,
            process_time                    => $time,
            source_environment_commands_ref => [$source_environment_cmd],
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
    $parameter_href->{$mip_program_name}{indirectory} = $outfamily_directory;

    ## Assign file_tags
    my $infile_tag = $file_info_href->{$family_id}{pvt}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$mip_program_name}{file_tag};
    my $infile_prefix       = $family_id . $infile_tag . $call_type;
    my $outfile_prefix      = $family_id . $outfile_tag . $call_type;
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
            jobid_chain    => $job_id_chain,
        }
    );
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
            job_id_chain   => $job_id_chain,
            file_suffix => $parameter_href->{$mip_program_name}{outfile_suffix},
        }
    );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    ($xargs_file_counter) = xargs_migrate_contig_files(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            contigs_ref        => \@{ $file_info_href->{contigs_size_ordered} },
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
            infile             => $infile_prefix,
            indirectory        => $infamily_directory,
            temp_directory     => $temp_directory,
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
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## Get parameters
        # VEP plugins
        my @plugins;
      PLUGIN:
        foreach my $plugin ( @{ $active_parameter_href->{vep_plugins} } ) {

            if ( $plugin eq q{LoF} ) {

                my $lof_parameter = q{,human_ancestor_fa:}
                  . catfile(
                    $active_parameter_href->{vep_directory_cache},
                    q{human_ancestor.fa,filter_position:0.05}
                  );
                push @plugins, $plugin . $lof_parameter;
            }
            elsif ( $plugin eq q{MaxEntScan} ) {

                my $lof_parameter = q{,}
                  . catfile( $active_parameter_href->{vep_directory_cache},
                    q{fordownload} );
                push @plugins, $plugin . $lof_parameter;
            }
            elsif ( $plugin eq q{UpDownDistance} ) {

                # Special case for mitochondrial contig annotation

                if ( $contig =~ / MT|M /xsm ) {

                    push @plugins, q{UpDownDistance,10,10};
                }
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
            if ( ( $contig =~ / MT|M /xsm ) && ( $vep_feature eq q{refseq} ) ) {

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
                regions_ref      => [$contig],
                plugins_ref      => \@plugins,
                vep_features_ref => \@vep_features_ref,
                assembly         => $assembly_version,
                cache_directory =>
                  $active_parameter_href->{vep_directory_cache},
                reference_path =>
                  $active_parameter_href->{human_genome_reference},
                infile_format   => substr( $outfile_suffix, 1 ),
                outfile_format  => substr( $outfile_suffix, 1 ),
                fork            => $VEP_FORK_NUMBER,
                buffer_size     => 20_000,
                infile_path     => $infile_path,
                outfile_path    => $outfile_path,
                stderrfile_path => $stderrfile_path,
                stdoutfile_path => $stdoutfile_path,
                FILEHANDLE      => $XARGSFILEHANDLE,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## QC Data File(s)
    migrate_file(
        {
            infile_path => $outfile_path_prefix
              . $UNDERSCORE
              . $ASTERIX
              . $infile_suffix
              . $UNDERSCORE . q{s}
              . $ASTERIX,
            outfile_path => $outfamily_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $XARGSFILEHANDLE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            infile_path => $outfile_path_prefix
              . $UNDERSCORE
              . $ASTERIX
              . $infile_suffix
              . $ASTERIX,
            outfile_path => $outfamily_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {

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
                sample_info_href => $sample_info_href,
                program_name     => $program_name,
                metafile_tag     => q{summary},
                path             => $qc_vep_summary_outfile_path,
            }
        );
        add_program_metafile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => $program_name,
                metafile_tag     => q{stderrfile},
                path             => catfile( $directory, $stderr_file ),
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
                sample_info_href => $sample_info_href,
                program_name     => $program_name,
                path             => $qc_vep_outfile_path,
            }
        );

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $family_id,
                path             => $job_id_chain,
                log              => $log,
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
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id
        },
        FILEHANDLE     => { required => 1, store => \$FILEHANDLE },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href
        },
        file_path               => { strict_type => 1, store => \$file_path },
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
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        sample_info_href  => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href
        },
        stderr_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$stderr_path
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter
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
      qw{ add_program_outfile_to_sample_info add_program_metafile_to_sample_info};

    ## Constants
    Readonly my $VEP_FORK_NUMBER => 4;

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

    my $reduce_io = $active_parameter_href->{reduce_io};
    my $xargs_file_path_prefix;

    ## Filehandles
    # Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            module_core_number =>
              $active_parameter_href->{module_core_number}{$mip_program_name},
            modifier_core_number => scalar @{ $file_info_href->{contigs} },
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    # Adjust for the number of forks vep forks
    $core_number = floor( $core_number / $VEP_FORK_NUMBER );

    ## If program needs special environment variables set
    if ($source_environment_cmd) {

        write_source_environment_command(
            {
                FILEHANDLE                      => $FILEHANDLE,
                source_environment_commands_ref => [$source_environment_cmd],
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
    $parameter_href->{$mip_program_name}{indirectory} = $outfamily_directory;

    ## Assign file_tags
    my $infile_tag = $file_info_href->{$family_id}{pvt}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$mip_program_name}{file_tag};
    my $infile_prefix       = $family_id . $infile_tag . $call_type;
    my $outfile_prefix      = $family_id . $outfile_tag . $call_type;
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
            jobid_chain    => $job_id_chain,
        }
    );
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
            job_id_chain   => $job_id_chain,
            file_suffix => $parameter_href->{$mip_program_name}{outfile_suffix},
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
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## Get parameters
        # VEP plugins
        my @plugins;
      PLUGIN:
        foreach my $plugin ( @{ $active_parameter_href->{vep_plugins} } ) {

            if ( $plugin eq q{LoF} ) {

                my $lof_parameter = q{,human_ancestor_fa:}
                  . catfile(
                    $active_parameter_href->{vep_directory_cache},
                    q{human_ancestor.fa,filter_position:0.05}
                  );
                push @plugins, $plugin . $lof_parameter;
            }
            elsif ( $plugin eq q{MaxEntScan} ) {

                my $lof_parameter = q{,}
                  . catfile( $active_parameter_href->{vep_directory_cache},
                    q{fordownload} );
                push @plugins, $plugin . $lof_parameter;
            }
            elsif ( $plugin eq q{UpDownDistance} ) {

                # Special case for mitochondrial contig annotation

                if ( $contig =~ /MT|M/xsm ) {

                    push @plugins, q{UpDownDistance,10,10};
                }
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
                regions_ref      => [$contig],
                plugins_ref      => \@plugins,
                vep_features_ref => \@vep_features_ref,
                assembly         => $assembly_version,
                cache_directory =>
                  $active_parameter_href->{vep_directory_cache},
                reference_path =>
                  $active_parameter_href->{human_genome_reference},
                infile_format   => substr( $outfile_suffix, 1 ),
                outfile_format  => substr( $outfile_suffix, 1 ),
                fork            => $VEP_FORK_NUMBER,
                buffer_size     => 20_000,
                infile_path     => $infile_path,
                outfile_path    => $outfile_path,
                stderrfile_path => $stderrfile_path,
                stdoutfile_path => $stdoutfile_path,
                FILEHANDLE      => $XARGSFILEHANDLE,
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
            infile_path => $outfile_path_prefix
              . $UNDERSCORE
              . $ASTERIX
              . $infile_suffix
              . $UNDERSCORE . q{s}
              . $ASTERIX,
            outfile_path => $outfamily_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    # Move file for downstream collection of VEP version

    migrate_file(
        {
            infile_path => $outfile_path_prefix
              . $UNDERSCORE
              . $file_info_href->{contigs_size_ordered}[0]
              . $infile_suffix,
            outfile_path => $outfamily_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Return to main or default environment using conda
    write_return_to_conda_environment(
        {
            source_main_environment_commands_ref =>
              \@{ $active_parameter_href->{source_main_environment_commands} },
            FILEHANDLE => $FILEHANDLE,
        }
    );

    if ( $mip_program_mode == 1 ) {

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
                sample_info_href => $sample_info_href,
                program_name     => $program_name,
                metafile_tag     => q{summary},
                path             => $qc_vep_summary_outfile_path,
            }
        );

        add_program_metafile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => $program_name,
                metafile_tag     => q{stderrfile},
                path             => catfile( $directory, $stderr_file ),
            }
        );
        my $qc_vep_outfile_path = catfile( $outfamily_directory,
                $outfile_prefix
              . $UNDERSCORE
              . $file_info_href->{contigs_size_ordered}[0]
              . $infile_suffix );
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => $program_name,
                path             => $qc_vep_outfile_path,
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
##          : $contigs_ref             => Contigs to analyse
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
    my $contigs_ref;
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
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href
        },
        call_type => {
            default     => q{SV},
            allow       => [qw{ SV }],
            strict_type => 1,
            store       => \$call_type
        },
        contigs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$contigs_ref
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id
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
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter
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
      qw{ add_program_outfile_to_sample_info add_program_metafile_to_sample_info };

    ## Constants
    Readonly my $VEP_FORK_NUMBER => 4;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
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
            module_core_number =>
              $active_parameter_href->{module_core_number}{$mip_program_name},
            modifier_core_number => scalar @{$contigs_ref},
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    # Adjust for the number of forks vep forks
    $core_number = floor( $core_number / $VEP_FORK_NUMBER );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            job_id_href                     => $job_id_href,
            FILEHANDLE                      => $FILEHANDLE,
            directory_id                    => $family_id,
            program_name                    => $program_name,
            program_directory               => catfile( lc $outaligner_dir ),
            call_type                       => $call_type,
            core_number                     => $core_number,
            process_time                    => $time,
            source_environment_commands_ref => [$source_environment_cmd],
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
    $parameter_href->{$mip_program_name}{indirectory} = $outfamily_directory;

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$family_id}{psv_combinevariantcallsets}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$mip_program_name}{file_tag};
    my $infile_prefix       = $family_id . $infile_tag . $call_type;
    my $outfile_prefix      = $family_id . $outfile_tag . $call_type;
    my $infile_path_prefix  = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ## Assign suffix
    my $file_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
            jobid_chain    => $job_id_chain,
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
            infile_path_prefix => $infile_path_prefix,
            file_suffix        => $file_suffix,
            FILEHANDLE         => $FILEHANDLE,
        }
    );

    ## varianteffectpredictor
    say {$FILEHANDLE} q{## varianteffectpredictor};

    my $assembly_version = $file_info_href->{human_genome_reference_source}
      . $file_info_href->{human_genome_reference_version};

    ## Get genome source and version to be compatible with VEP
    $assembly_version =
      _get_assembly_name( { assembly_version => $assembly_version } );

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig ( @{$contigs_ref} ) {

        ## Get parameters
        my $vep_outfile_prefix         = $outfile_prefix;
        my $vep_xargs_file_path_prefix = $xargs_file_path_prefix;
        my @regions;

        ## Contig specific
        # Update endings with contig info
        if (   ( $consensus_analysis_type eq q{wgs} )
            || ( $consensus_analysis_type eq q{mixed} ) )
        {

            $vep_outfile_prefix = $outfile_prefix . $UNDERSCORE . $contig;
            $vep_xargs_file_path_prefix =
              $xargs_file_path_prefix . $DOT . $contig;
            push @regions, $contig;
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
            elsif ( $plugin eq q{UpDownDistance} ) {

                # Special case for mitochondrial contig annotation

                if ( $contig =~ /MT|M/xsm ) {

                    push @plugins, q{UpDownDistance,10,10};
                }
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
                regions_ref      => \@regions,
                plugins_ref      => \@plugins,
                vep_features_ref => \@vep_features_ref,
                assembly         => $assembly_version,
                cache_directory =>
                  $active_parameter_href->{vep_directory_cache},
                reference_path =>
                  $active_parameter_href->{human_genome_reference},
                infile_format   => substr( $file_suffix, 1 ),
                outfile_format  => substr( $file_suffix, 1 ),
                fork            => $VEP_FORK_NUMBER,
                buffer_size     => 100,
                infile_path     => $infile_path,
                outfile_path    => $outfile_path,
                stderrfile_path => $stderrfile_path,
                stdoutfile_path => $stdoutfile_path,
                FILEHANDLE      => $XARGSFILEHANDLE,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;

# Only perform once for exome samples to avoid risking contigs lacking variants throwing errors
        if ( $consensus_analysis_type eq q{wes} ) {

            last CONTIG;
        }
    }

    close $XARGSFILEHANDLE;

    ## QC Data File(s)
    migrate_file(
        {
            infile_path => $outfile_path_prefix
              . $ASTERIX
              . $file_suffix
              . $UNDERSCORE . q{s}
              . $ASTERIX,
            outfile_path => $outfamily_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            infile_path => $outfile_path_prefix
              . $ASTERIX
              . $file_suffix
              . $ASTERIX,
            outfile_path => $outfamily_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {

        my $outfile_sample_info_prefix = $outfile_prefix;

        # Update endings with contig info
        if (   ( $consensus_analysis_type eq q{wgs} )
            || ( $consensus_analysis_type eq q{mixed} ) )
        {

            $outfile_sample_info_prefix .= $UNDERSCORE . $contigs_ref->[0];
        }

        ## Collect QC metadata info for later use
        my $qc_vep_summary_outfile =
          $outfile_sample_info_prefix . $DOT . q{vcf_summary.html};
        add_program_metafile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => $program_name,
                metafile_tag     => q{summary},
                directory        => $outfamily_directory,
                file             => $qc_vep_summary_outfile,
            }
        );

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => $program_name,
                outdirectory     => $outfamily_directory,
                outfile          => $outfile_sample_info_prefix . $file_suffix,
            }
        );

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $family_id,
                path             => $job_id_chain,
                log              => $log,
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
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$assembly_version
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
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path_prefix
        },
        file_suffix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_suffix
        },
        FILEHANDLE => { store => \$FILEHANDLE },
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

1;
