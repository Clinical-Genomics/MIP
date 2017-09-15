package MIP::Recipes::Vep;

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
use File::Basename qw{ fileparse };
use File::Spec::Functions qw{ catdir catfile };
use POSIX;

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw( analysis_vep );

}

##Constants
Readonly my $ASTERIX    => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_vep {

##analysis_vep

##Function : Varianteffectpredictor performs effect predictions and annotation of variants.
##Returns  : "|$xargs_file_counter"
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, $program_info_path, $file_path, $stderr_path, $FILEHANDLE, family_id, $temp_directory, $outaligner_dir, $call_type, $xargs_file_counter
##         : $parameter_href          => Parameter hash {REF}
##         : $active_parameter_href   => Active parameters for this analysis hash {REF}
##         : $sample_info_href        => Info on samples and family hash {REF}
##         : $file_info_href          => The file_info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href             => Job id hash {REF}
##         : $program_name            => Program name
##         : $program_info_path       => The program info path
##         : $file_path               => File path
##         : $stderr_path             => Stderr path of the block script
##         : $FILEHANDLE              => Filehandle to write to
##         : $family_id               => Family id
##         : $temp_directory          => Temporary directory
##         : $outaligner_dir          => Outaligner_dir used in the analysis
##         : $call_type               => The variant call type
##         : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id;
    my $temp_directory;
    my $outaligner_dir;
    my $call_type;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;
    my $program_info_path;
    my $file_path;
    my $stderr_path;
    my $FILEHANDLE;

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
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        file_path         => { strict_type => 1, store => \$file_path },
        stderr_path       => { strict_type => 1, store => \$stderr_path },
        FILEHANDLE => { store => \$FILEHANDLE },
        family_id  => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id
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
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw{ get_core_number };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Set::File qw{ set_file_suffix };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Recipes::Xargs qw{ xargs_command };
    use MIP::Program::Variantcalling::Vep qw{ variant_effect_predictor };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info add_program_metafile_to_sample_info};
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };

    ##Constants
    Readonly my $VEP_FORK_NUMBER => 4;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $time         = $active_parameter_href->{module_time}{$mip_program_name};
    my $reduce_io    = $active_parameter_href->{reduce_io};
    my $xargs_file_path_prefix;

    ## Filehandles
    # Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();

    # Run as individual sbatch script
    if ( not defined $FILEHANDLE ) {

        # Create anonymous filehandle
        $FILEHANDLE = IO::Handle->new();
    }

    ## Get core number depending on user supplied input exists or not and max number of cores
    my $core_number = get_core_number(
        {
            module_core_number =>
              $active_parameter_href->{module_core_number}{$mip_program_name},
            modifier_core_number => scalar @{ $file_info_href->{contigs} },
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    # Adjust for the number of forks vep forks
    $core_number = floor( $core_number / $VEP_FORK_NUMBER );

    # Run as individual sbatch script
    if ( not $reduce_io ) {

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        ( $file_path, $program_info_path ) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                job_id_href           => $job_id_href,
                FILEHANDLE            => $FILEHANDLE,
                directory_id          => $family_id,
                program_name          => $program_name,
                program_directory     => $outaligner_dir,
                call_type             => $call_type,
                core_number           => $core_number,
                process_time          => $time,
                temp_directory        => $temp_directory
            }
        );
        $stderr_path = $program_info_path . $DOT . q{stderr.txt};
    }

    # Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) =
      File::Spec->splitpath($stderr_path);

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

    # Run as individual sbatch script
    if ( not $reduce_io ) {

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        ($xargs_file_counter) = xargs_migrate_contig_files(
            {
                FILEHANDLE      => $FILEHANDLE,
                XARGSFILEHANDLE => $XARGSFILEHANDLE,
                contigs_ref => \@{ $file_info_href->{contigs_size_ordered} },
                file_path   => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
                xargs_file_counter => $xargs_file_counter,
                infile             => $infile_prefix,
                indirectory        => $infamily_directory,
                temp_directory     => $temp_directory,
            }
        );
    }

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

  CONTIGS:
    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## Get parameters
        # VEP plugins
        my @plugins;
      PLUGINS:
        foreach my $plugin ( @{ $active_parameter_href->{vep_plugins} } ) {

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

        ## VEP features
        my @vep_features_ref;
      FEATURES:
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
                script_path      => catfile(
                    $active_parameter_href->{vep_directory_path},
                    q{variant_effect_predictor.pl}
                ),
                assembly => $assembly_version,
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

    # Run as individual sbatch script
    if ( not $reduce_io ) {

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
    }
    else {
        # Move file for downstream collection of VEP version

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
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
    }

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

        # Run as individual sbatch script
        if ( not $reduce_io ) {

            slurm_submit_job_sample_id_dependency_add_to_family(
                {
                    job_id_href             => $job_id_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    sample_ids_ref =>
                      \@{ $active_parameter_href->{sample_ids} },
                    family_id        => $family_id,
                    path             => $job_id_chain,
                    log              => $log,
                    sbatch_file_name => $file_path,
                }
            );
        }
        if ($reduce_io) {

    # Redirect qccollect search to Block File, since VEP will write stderr there
            $program_name = q{variantannotationblock};

            ## Collect QC metadata info for later use
            add_program_outfile_to_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    program_name     => $program_name,
                    path             => catfile( $directory, $stderr_file ),
                }
            );
        }
    }

    # Track the number of created xargs scripts per module for Block algorithm
    return $xargs_file_counter;
}

sub _get_assembly_name {

##_get_assembly_name

##Function : Get genome source and version to be compatible with VEP
##Returns  : $assembly_version
##Arguments: $assembly_version
##         : $assembly_version => The genome source and version to be checked

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

1;
