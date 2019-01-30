package MIP::Recipes::Analysis::Vep;

use 5.026;
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
    our $VERSION = 1.14;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ analysis_vep analysis_vep_rna analysis_vep_sv_wes analysis_vep_sv_wgs };

}

## Constants
Readonly my $ASTERISK               => q{*};
Readonly my $ANNOTATION_DISTANCE    => 5000;
Readonly my $ANNOTATION_DISTANCE_MT => 0;
Readonly my $COMMA                  => q{,};
Readonly my $DOT                    => q{.};
Readonly my $EMPTY_STR              => q{};
Readonly my $NEWLINE                => qq{\n};
Readonly my $SPACE                  => q{ };
Readonly my $UNDERSCORE             => q{_};

sub analysis_vep {

## Function : Varianteffectpredictor performs effect predictions and annotation of variants.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
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
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
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
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
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
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
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
            allow       => qr{ \A\d+\z }xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw{ get_core_number };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_parameters get_recipe_attributes };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Variantcalling::Vep qw{ variant_effect_predictor };
    use MIP::QC::Sample_info
      qw{ set_recipe_metafile_in_sample_info set_recipe_outfile_in_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constants
    Readonly my $VEP_FORK_NUMBER => 4;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $case_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my %infile_path        = %{ $io{in}{file_path_href} };
    my $infile_suffix      = $io{in}{file_suffix};

    my @contigs_size_ordered = @{ $file_info_href->{contigs_size_ordered} };
    my $job_id_chain         = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode = $active_parameter_href->{$recipe_name};
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $xargs_file_path_prefix;

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $case_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => \@contigs_size_ordered,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
                temp_directory   => $temp_directory,
            }
        )
    );

    my %outfile_path   = %{ $io{out}{file_path_href} };
    my @outfile_paths  = @{ $io{out}{file_paths} };
    my $outfile_suffix = $io{out}{file_suffix};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            max_cores_per_node   => $active_parameter_href->{max_cores_per_node},
            modifier_core_number => scalar @{ $file_info_href->{contigs} },
            recipe_core_number =>
              $active_parameter_href->{recipe_core_number}{$recipe_name},
        }
    );

    # Adjust for the number of forks vep forks
    $core_number = floor( $core_number / $VEP_FORK_NUMBER );

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $case_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );
    my $stderr_path = $recipe_info_path . $DOT . q{stderr.txt};

    ### SHELL:

    ## Varianteffectpredictor
    say {$FILEHANDLE} q{## Varianteffectpredictor};

    my $assembly_version = $file_info_href->{human_genome_reference_source}
      . $file_info_href->{human_genome_reference_version};

    ## Get genome source and version to be compatible with VEP
    $assembly_version = _get_assembly_name( { assembly_version => $assembly_version } );

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Get parameters
    # VEP custom annotations
    my @custom_annotations;
    if ( exists $active_parameter_href->{vep_custom_annotation} ) {

        @custom_annotations = _get_custom_annotation_cmds(
            {
                vep_custom_annotation_href =>
                  $active_parameter_href->{vep_custom_annotation},
            }
        );
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

            my $max_ent_scan_dir_path =
              catfile( $active_parameter_href->{vep_directory_cache},
                qw{ Plugins fordownload } );
            my @max_ent_scan_de_novos = qw{ SWA NCSS };

            push @plugins, join $COMMA,
              ( $plugin, $max_ent_scan_dir_path, @max_ent_scan_de_novos );
        }
        elsif ( $plugin eq q{ExACpLI}
            and exists $active_parameter_href->{vep_plugin_pli_value_file_path} )
        {

            my $pli_file_path =
              q{,} . $active_parameter_href->{vep_plugin_pli_value_file_path};
            push @plugins, $plugin . $pli_file_path;
        }
        else {

            push @plugins, $plugin;
        }
    }

  CONTIG:
    foreach my $contig (@contigs_size_ordered) {

        ## Get contig specific parameters
        my $distance = $ANNOTATION_DISTANCE;

        # Special case for mitochondrial contig annotation
        if ( $contig =~ / MT|M /xsm ) {

            $distance = $ANNOTATION_DISTANCE_MT;
        }

        ## VEP features
        my @vep_features_ref;
      FEATURE:
        foreach my $vep_feature ( @{ $active_parameter_href->{vep_features} } ) {

            # Add VEP features to the output.
            push @vep_features_ref, $vep_feature;

            # Special case for mitochondrial contig annotation
            if ( $contig =~ / MT|M /xsm && $vep_feature eq q{refseq} ) {

                push @vep_features_ref, q{all_refseq};
            }
        }

        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        my $stdoutfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stdout.txt};
        variant_effect_predictor(
            {
                assembly               => $assembly_version,
                buffer_size            => 20_000,
                cache_directory        => $active_parameter_href->{vep_directory_cache},
                custom_annotations_ref => \@custom_annotations,
                distance               => $distance,
                FILEHANDLE             => $XARGSFILEHANDLE,
                fork                   => $VEP_FORK_NUMBER,
                infile_format          => substr( $infile_suffix, 1 ),
                infile_path            => $infile_path{$contig},
                outfile_format         => substr( $outfile_suffix, 1 ),
                outfile_path           => $outfile_path{$contig},
                plugins_dir_path       => $active_parameter_href->{vep_plugins_dir_path},
                plugins_ref            => \@plugins,
                reference_path   => $active_parameter_href->{human_genome_reference},
                regions_ref      => [$contig],
                stderrfile_path  => $stderrfile_path,
                stdoutfile_path  => $stdoutfile_path,
                vep_features_ref => \@vep_features_ref,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});
    close $XARGSFILEHANDLE
      or $log->logcroak(q{Could not close XARGSFILEHANDLE});

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_metafile_in_sample_info(
            {
                metafile_tag     => q{stderrfile},
                path             => $stderr_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_paths[0],
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                dependency_method       => q{sample_to_case},
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                job_id_chain            => $job_id_chain,
                recipe_file_path        => $recipe_file_path,
                sample_ids_ref          => \@{ $active_parameter_href->{sample_ids} },
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return;
}

sub analysis_vep_sv_wes {

## Function : Varianteffectpredictor annotation of exome SV variants.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id               => Family id
##          : $FILEHANDLE              => Filehandle to write to
##          : $file_info_href          => The file_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name            => Program name
##          : $recipe_info_path       => Recipe info path
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
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
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
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
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
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
            allow       => qr{ \A\d+\z }xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw{ get_core_number };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_parameters get_recipe_attributes };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Variantcalling::Vep qw{ variant_effect_predictor };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::QC::Sample_info
      qw{ set_recipe_metafile_in_sample_info set_recipe_outfile_in_sample_info };

    ### PREPROCESSING:

    ## Constants
    Readonly my $VEP_FORK_NUMBER => 4;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $case_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );

    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path_prefix = $io{in}{file_path_prefix};
    my $infile_suffix      = $io{in}{file_suffix};
    my $infile_path        = $infile_path_prefix . $infile_suffix;

    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my $job_id_chain            = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode = $active_parameter_href->{$recipe_name};
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $xargs_file_path_prefix;
    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $case_id,
                file_info_href         => $file_info_href,
                outdata_dir            => $active_parameter_href->{outdata_dir},
                file_name_prefixes_ref => [$infile_name_prefix],
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
                temp_directory         => $temp_directory,
            }
        )
    );

    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my $outfile_path        = $outfile_path_prefix . $outfile_suffix;

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            max_cores_per_node   => $active_parameter_href->{max_cores_per_node},
            modifier_core_number => 1,
            recipe_core_number =>
              $active_parameter_href->{recipe_core_number}{$recipe_name},
        }
    );

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $case_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );
    my $stderr_path = $recipe_info_path . $DOT . q{stderr.txt};

    # Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) = splitpath($stderr_path);

    ### SHELL:

    ## Reformat SV with no length as these will fail in the annotation with VEP
    _reformat_sv_with_no_length(
        {
            FILEHANDLE         => $FILEHANDLE,
            file_suffix        => $infile_suffix,
            infile_path_prefix => $infile_path_prefix,
        }
    );

    ## varianteffectpredictor
    say {$FILEHANDLE} q{## Varianteffectpredictor};

    my $assembly_version = $file_info_href->{human_genome_reference_source}
      . $file_info_href->{human_genome_reference_version};

    ## Get genome source and version to be compatible with VEP
    $assembly_version = _get_assembly_name( { assembly_version => $assembly_version, } );

    # VEP custom annotations
    my @custom_annotations;
    if ( exists $active_parameter_href->{vep_custom_annotation} ) {

        @custom_annotations = _get_custom_annotation_cmds(
            {
                vep_custom_annotation_href =>
                  $active_parameter_href->{vep_custom_annotation},
            }
        );
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
        elsif ( $plugin eq q{ExACpLI}
            and exists $active_parameter_href->{vep_plugin_pli_value_file_path} )
        {

            my $pli_file_path =
              q{,} . $active_parameter_href->{vep_plugin_pli_value_file_path};
            push @plugins, $plugin . $pli_file_path;
        }
        else {

            push @plugins, $plugin;
        }
    }

    ## VEP features
    my @vep_features_ref;

  FEATURE:
    foreach my $vep_feature ( @{ $active_parameter_href->{sv_vep_features} } ) {

        # Add VEP features to the output.
        push @vep_features_ref, $vep_feature;
    }

    my $vep_infile_path =
      $infile_path_prefix . $UNDERSCORE . q{fixedsvlength} . $infile_suffix;
    my $stderrfile_path = $recipe_file_path . $DOT . q{stderr.txt};
    my $stdoutfile_path = $recipe_file_path . $DOT . q{stdout.txt};
    variant_effect_predictor(
        {
            assembly         => $assembly_version,
            buffer_size      => 100,
            cache_directory  => $active_parameter_href->{vep_directory_cache},
            FILEHANDLE       => $FILEHANDLE,
            fork             => $VEP_FORK_NUMBER,
            infile_format    => substr( $infile_suffix, 1 ),
            infile_path      => $vep_infile_path,
            outfile_format   => substr( $outfile_suffix, 1 ),
            outfile_path     => $outfile_path,
            plugins_dir_path => $active_parameter_href->{vep_plugins_dir_path},
            plugins_ref      => \@plugins,
            reference_path   => $active_parameter_href->{human_genome_reference},
            stderrfile_path  => $stderrfile_path,
            stdoutfile_path  => $stdoutfile_path,
            vep_features_ref => \@vep_features_ref,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    close $FILEHANDLE;

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );
        submit_recipe(
            {
                dependency_method       => q{sample_to_case},
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                job_id_chain            => $job_id_chain,
                recipe_file_path        => $recipe_file_path,
                sample_ids_ref          => \@{ $active_parameter_href->{sample_ids} },
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return;
}

sub analysis_vep_sv_wgs {

## Function : Varianteffectpredictor performs annotation of wgs SV variants.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id               => Family id
##          : $FILEHANDLE              => Filehandle to write to
##          : $file_info_href          => The file_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name            => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
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
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
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
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
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
            allow       => qr{ \A\d+\z }xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw{ get_core_number };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_parameters get_recipe_attributes };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Variantcalling::Vep qw{ variant_effect_predictor };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::QC::Sample_info
      qw{ set_recipe_metafile_in_sample_info set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constants
    Readonly my $VEP_FORK_NUMBER => 4;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $case_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );

    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path_prefix = $io{in}{file_path_prefix};
    my $infile_suffix      = $io{in}{file_suffix};
    my $infile_path        = $infile_path_prefix . $infile_suffix;

    my $contigs_ref             = \@{ $file_info_href->{contigs} };
    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my $job_id_chain            = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode = $active_parameter_href->{$recipe_name};
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $xargs_file_path_prefix;

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $case_id,
                file_info_href   => $file_info_href,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => $file_info_href->{contigs_size_ordered},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
                temp_directory   => $temp_directory,
            }
        )
    );

    my $outdir_path_prefix = $io{out}{dir_path_prefix};
    my @outfile_paths      = @{ $io{out}{file_paths} };
    my $outfile_suffix     = $io{out}{file_suffix};
    my %outfile_path       = %{ $io{out}{file_path_href} };

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            max_cores_per_node   => $active_parameter_href->{max_cores_per_node},
            modifier_core_number => scalar @{ $file_info_href->{contigs_size_ordered} },
            recipe_core_number =>
              $active_parameter_href->{recipe_core_number}{$recipe_name},
        }
    );

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $case_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );
    my $stderr_path = $recipe_info_path . $DOT . q{stderr.txt};

    # Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) = splitpath($stderr_path);

    ### SHELL:
    ## Reformat SV with no length as these will fail in the annotation with VEP
    _reformat_sv_with_no_length(
        {
            FILEHANDLE         => $FILEHANDLE,
            file_suffix        => $infile_suffix,
            infile_path_prefix => $infile_path_prefix,
        }
    );

    ## varianteffectpredictor
    say {$FILEHANDLE} q{## Varianteffectpredictor};

    my $assembly_version = $file_info_href->{human_genome_reference_source}
      . $file_info_href->{human_genome_reference_version};

    ## Get genome source and version to be compatible with VEP
    $assembly_version = _get_assembly_name( { assembly_version => $assembly_version, } );

    my $vep_infile_path =
      $infile_path_prefix . $UNDERSCORE . q{fixedsvlength} . $infile_suffix;

    # VEP custom annotations
    my @custom_annotations;
    if ( exists $active_parameter_href->{vep_custom_annotation} ) {

        @custom_annotations = _get_custom_annotation_cmds(
            {
                vep_custom_annotation_href =>
                  $active_parameter_href->{vep_custom_annotation},
            }
        );
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
        elsif ( $plugin eq q{ExACpLI}
            and exists $active_parameter_href->{vep_plugin_pli_value_file_path} )
        {

            my $pli_file_path =
              q{,} . $active_parameter_href->{vep_plugin_pli_value_file_path};
            push @plugins, $plugin . $pli_file_path;
        }
        else {

            push @plugins, $plugin;
        }
    }

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig ( @{$contigs_ref} ) {

        ## Get contig specific parameters
        my $distance = $ANNOTATION_DISTANCE;

        # Special case for mitochondrial contig annotation
        if ( $contig =~ / MT|M /xsm ) {

            $distance = $ANNOTATION_DISTANCE_MT;
        }

        my $mt_name;
        my @regions;
        my $vep_xargs_file_path_prefix = $xargs_file_path_prefix . $DOT . $contig;
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

        ## VEP features
        my @vep_features_ref;

      FEATURE:
        foreach my $vep_feature ( @{ $active_parameter_href->{sv_vep_features} } ) {

            # Add VEP features to the output.
            push @vep_features_ref, $vep_feature;

            # Special case for mitochondrial contig annotation
            if ( $contig =~ /MT|M/sxm && $vep_feature eq q{refseq} ) {

                push @vep_features_ref, q{all_refseq};
            }
        }

        my $stderrfile_path = $vep_xargs_file_path_prefix . $DOT . q{stderr.txt};
        my $stdoutfile_path = $vep_xargs_file_path_prefix . $DOT . q{stdout.txt};
        variant_effect_predictor(
            {
                assembly         => $assembly_version,
                buffer_size      => 100,
                cache_directory  => $active_parameter_href->{vep_directory_cache},
                distance         => $distance,
                FILEHANDLE       => $XARGSFILEHANDLE,
                fork             => $VEP_FORK_NUMBER,
                infile_format    => substr( $infile_suffix, 1 ),
                infile_path      => $vep_infile_path,
                outfile_format   => substr( $outfile_suffix, 1 ),
                outfile_path     => $outfile_path{$contig},
                plugins_dir_path => $active_parameter_href->{vep_plugins_dir_path},
                plugins_ref      => \@plugins,
                regions_ref      => \@regions,
                reference_path   => $active_parameter_href->{human_genome_reference},
                stderrfile_path  => $stderrfile_path,
                stdoutfile_path  => $stdoutfile_path,
                vep_features_ref => \@vep_features_ref,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;

        ## Filter out the MT annotations from the combined chr21 chrM call
        if ($mt_name) {

            say {$FILEHANDLE} q{## Filter out MT annotations};
            _subset_vcf(
                {
                    FILEHANDLE   => $FILEHANDLE,
                    infile_path  => $outfile_path{$contig},
                    outfile_path => $outfile_path{$contig},
                    regions_ref  => [qw{ chrM MT }],
                }
            );
        }
    }

    close $XARGSFILEHANDLE;
    close $FILEHANDLE;

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_paths[0],
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );
        submit_recipe(
            {
                dependency_method       => q{sample_to_case},
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                job_id_chain            => $job_id_chain,
                recipe_file_path        => $recipe_file_path,
                sample_ids_ref          => \@{ $active_parameter_href->{sample_ids} },
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return;
}

sub analysis_vep_rna {

## Function : Varianteffectpredictor performs effect predictions and annotation of variantsi from RNA-seq data.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
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
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
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
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
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

    use MIP::Cluster qw{ get_core_number };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_parameters get_recipe_attributes };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Variantcalling::Vep qw{ variant_effect_predictor };
    use MIP::QC::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constants
    Readonly my $VEP_FORK_NUMBER => 4;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $case_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path        = $io{in}{file_path};
    my $infile_suffix      = $io{in}{file_suffix};

    my $job_id_chain = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode = $active_parameter_href->{$recipe_name};
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set and get the io files per chain, id and stream
    %io = parse_io_outfiles(
        {
            chain_id               => $job_id_chain,
            id                     => $case_id,
            file_info_href         => $file_info_href,
            file_name_prefixes_ref => [$infile_name_prefix],
            outdata_dir            => $active_parameter_href->{outdata_dir},
            parameter_href         => $parameter_href,
            recipe_name            => $recipe_name,
        }
    );

    my $outfile_name   = $io{out}{file_names}->[0];
    my $outfile_path   = $io{out}{file_path};
    my $outfile_suffix = $io{out}{file_suffix};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            max_cores_per_node   => $active_parameter_href->{max_cores_per_node},
            modifier_core_number => scalar @{ $file_info_href->{contigs} },
            recipe_core_number =>
              $active_parameter_href->{recipe_core_number}{$recipe_name},
        }
    );

    # Adjust for the number of forks vep forks
    $core_number = floor( $core_number / $VEP_FORK_NUMBER );

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $case_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => \@source_environment_cmds,
        }
    );

    ### SHELL:

    ## Varianteffectpredictor
    say {$FILEHANDLE} q{## Varianteffectpredictor};

    my $assembly_version = $file_info_href->{human_genome_reference_source}
      . $file_info_href->{human_genome_reference_version};

    ## Get genome source and version to be compatible with VEP
    $assembly_version = _get_assembly_name( { assembly_version => $assembly_version } );

    ## Get parameters
    # VEP custom annotations
    my @custom_annotations;
    if ( exists $active_parameter_href->{vep_custom_annotation} ) {

        @custom_annotations = _get_custom_annotation_cmds(
            {
                vep_custom_annotation_href =>
                  $active_parameter_href->{vep_custom_annotation},
            }
        );
    }

    ## VEP features
    my @vep_features_ref;
  FEATURE:
    foreach my $vep_feature ( @{ $active_parameter_href->{vep_features} } ) {

        # Add VEP features to the output.
        push @vep_features_ref, $vep_feature;
    }

    variant_effect_predictor(
        {
            assembly               => $assembly_version,
            buffer_size            => 20_000,
            cache_directory        => $active_parameter_href->{vep_directory_cache},
            custom_annotations_ref => \@custom_annotations,
            distance               => $ANNOTATION_DISTANCE,
            FILEHANDLE             => $FILEHANDLE,
            fork                   => $VEP_FORK_NUMBER,
            infile_format          => substr( $infile_suffix, 1 ),
            infile_path            => $infile_path,
            outfile_format         => substr( $outfile_suffix, 1 ),
            outfile_path           => $outfile_path,
            reference_path         => $active_parameter_href->{human_genome_reference},
            vep_features_ref       => \@vep_features_ref,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        add_recipe_outfile_to_sample_info(
            {
                infile           => $outfile_name,
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                dependency_method       => q{sample_to_case},
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                job_id_chain            => $job_id_chain,
                recipe_file_path        => $recipe_file_path,
                sample_ids_ref          => \@{ $active_parameter_href->{sample_ids} },
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return;
}

sub _get_custom_annotation_cmds {

## Function : Build the custom annotation command per file for vep
## Returns  : @custom_annotations
## Arguments: $vep_custom_annotation_href => Custom annotation info {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vep_custom_annotation_href;

    my $tmpl = {
        vep_custom_annotation_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$vep_custom_annotation_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @custom_annotations;
    my @order_custom_options =
      qw{ path key file_type annotation_type force_report_coordinates };

  ANNOTATION:
    foreach my $annotation_href ( values %{$vep_custom_annotation_href} ) {

        my $cmd = join $COMMA, @{$annotation_href}{@order_custom_options};
        push @custom_annotations, $cmd;
    }
    return @custom_annotations;
}

sub _get_assembly_name {

## Function : Get genome source and version to be compatible with VEP
## Returns  : $assembly_version
## Arguments: $assembly_version => Genome source and version to be checked

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
## Arguments: $FILEHANDLE         => Filehandle to write to
##          : $file_suffix        => File suffix
##          : $infile_path_prefix => Infile path prefix

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $file_suffix;
    my $infile_path_prefix;

    my $tmpl = {
        FILEHANDLE  => { required => 1, store => \$FILEHANDLE, },
        file_suffix => {
            defined     => 1,
            required    => 1,
            store       => \$file_suffix,
            strict_type => 1,
        },
        infile_path_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path_prefix,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Fix SV with no length as these will fail in the annotation with VEP
    my $perl_fix_sv_nolengths;

    # Set up perl
    $perl_fix_sv_nolengths .= q?perl -nae '?;

    # Initate variables
    $perl_fix_sv_nolengths .= q?my %info; my $start; my $end; my $alt; my @data; ?;

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
    $perl_fix_sv_nolengths .= q?if(defined($info{END})) { $end = $info{END}; } ?;

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
    use MIP::Program::Variantcalling::Bcftools qw{ bcftools_index bcftools_view };

    ## Prepare for bcftools_view
    htslib_bgzip(
        {
            FILEHANDLE  => $FILEHANDLE,
            force       => 1,
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
