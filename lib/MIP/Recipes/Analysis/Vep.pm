package MIP::Recipes::Analysis::Vep;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ fileparse };
use File::Spec::Functions qw{ catfile splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use POSIX;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ %ANALYSIS $ASTERISK $COMMA $DOT $EMPTY_STR $LOG_NAME $NEWLINE %PRIMARY_CONTIG $SPACE $UNDERSCORE };
use MIP::File::Format::Vep qw{ create_vep_synonyms_file };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      analysis_vep
      analysis_vep_wgs
      analysis_vep_sv_wes
      analysis_vep_sv_wgs
    };

}

## Constants
Readonly my $ANNOTATION_DISTANCE    => $ANALYSIS{ANNOTATION_DISTANCE};
Readonly my $ANNOTATION_DISTANCE_MT => $ANALYSIS{ANNOTATION_DISTANCE_MT};
Readonly my $BUFFER_SIZE            => 20_000;
Readonly my $BUFFER_SIZE_SV         => 100;

sub analysis_vep_wgs {

## Function : Varianteffectpredictor performs effect predictions and annotation of variants.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $file_path               => File path
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $file_path;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
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
        file_path   => { store => \$file_path, strict_type => 1, },
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
        profile_base_command => {
            default     => q{sbatch},
            store       => \$profile_base_command,
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

    use MIP::Cluster qw{ get_core_number update_memory_allocation };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Vep qw{ variant_effect_predictor };
    use MIP::Sample_info qw{
      set_file_path_to_store
      set_recipe_metafile_in_sample_info
      set_recipe_outfile_in_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constants
    Readonly my $VEP_FORK_NUMBER => 4;

    my $log = Log::Log4perl->get_logger($LOG_NAME);

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
    my %infile_path        = %{ $io{in}{file_path_href} };
    my $infile_suffix      = $io{in}{file_suffix};

    my @contigs_size_ordered     = @{ $file_info_href->{contigs_size_ordered} };
    my $genome_reference_version = $file_info_href->{human_genome_reference_version};
    my $job_id_chain             = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number = $recipe_resource{core_number};
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
            }
        )
    );
    my $outdir_path_prefix = $io{out}{dir_path_prefix};
    my %outfile_path       = %{ $io{out}{file_path_href} };
    my @outfile_paths      = @{ $io{out}{file_paths} };
    my $outfile_suffix     = $io{out}{file_suffix};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle      = IO::Handle->new();
    my $xargsfilehandle = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            max_cores_per_node   => $active_parameter_href->{max_cores_per_node},
            modifier_core_number => scalar @{ $file_info_href->{contigs} },
            recipe_core_number   => $core_number,
        }
    );

    # Adjust for the number of forks vep forks
    my $parallel_processes = floor( $core_number / $VEP_FORK_NUMBER );

    ## Update memory depending on how many cores that are being used
    my $memory_allocation = update_memory_allocation(
        {
            node_ram_memory           => $active_parameter_href->{node_ram_memory},
            parallel_processes        => $parallel_processes,
            process_memory_allocation => $recipe_resource{memory},
        }
    );

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $case_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $memory_allocation,
            process_time          => $recipe_resource{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            temp_directory        => $temp_directory,
        }
    );
    my $stderr_path = $recipe_info_path . $DOT . q{stderr.txt};

    ### SHELL:

    ## Get the vep synonyms file path for if required (grch38)
    my $vep_synonyms_file_path = create_vep_synonyms_file(
        {
            outfile_path => catfile( $outdir_path_prefix, q{synonyms.tsv} ),
            version      => $genome_reference_version,
        }
    );

    ## Varianteffectpredictor
    say {$filehandle} q{## Varianteffectpredictor};

    my $assembly_version =
      $file_info_href->{human_genome_reference_source} . $genome_reference_version;

    ## Get genome source and version to be compatible with VEP
    $assembly_version = _get_assembly_name( { assembly_version => $assembly_version } );

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $parallel_processes,
            filehandle         => $filehandle,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            xargsfilehandle    => $xargsfilehandle,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Get parameters
    # VEP custom annotations
    my @custom_annotations = _get_custom_annotation_cmds(
        {
            vep_custom_annotation_href => $active_parameter_href->{vep_custom_annotation},
        }
    );

    # VEP plugins
    my @plugins =
      _get_plugin_cmds( { vep_plugin_href => $active_parameter_href->{vep_plugin}, } );

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

        my $stderrfile_path = $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        my $stdoutfile_path = $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stdout.txt};
        variant_effect_predictor(
            {
                assembly               => $assembly_version,
                buffer_size            => $BUFFER_SIZE,
                cache_directory        => $active_parameter_href->{vep_directory_cache},
                custom_annotations_ref => \@custom_annotations,
                distance               => $distance,
                filehandle             => $xargsfilehandle,
                fork                   => $VEP_FORK_NUMBER,
                infile_format          => substr( $infile_suffix, 1 ),
                infile_path            => $infile_path{$contig},
                outfile_format         => substr( $outfile_suffix, 1 ),
                outfile_path           => $outfile_path{$contig},
                plugins_dir_path       => $active_parameter_href->{vep_plugins_dir_path},
                plugins_ref            => \@plugins,
                reference_path         => $active_parameter_href->{human_genome_reference},
                regions_ref            => [$contig],
                stderrfile_path        => $stderrfile_path,
                stdoutfile_path        => $stdoutfile_path,
                synonyms_file_path     => $vep_synonyms_file_path,
                vep_features_ref       => \@vep_features_ref,
            }
        );
        say {$xargsfilehandle} $NEWLINE;
    }

    close $filehandle or $log->logcroak(q{Could not close filehandle});
    close $xargsfilehandle
      or $log->logcroak(q{Could not close xargsfilehandle});

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
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_case},
                job_id_chain                      => $job_id_chain,
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
                log                               => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_ids_ref     => \@{ $active_parameter_href->{sample_ids} },
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub analysis_vep_sv_wes {

## Function : Varianteffectpredictor annotation of exome SV variants.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $filehandle              => Filehandle to write to
##          : $file_info_href          => The file_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $recipe_info_path        => Recipe info path
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
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
        profile_base_command => {
            default     => q{sbatch},
            store       => \$profile_base_command,
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

    use MIP::Cluster qw{ get_core_number update_memory_allocation };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Vep qw{ variant_effect_predictor };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Sample_info qw{ set_recipe_metafile_in_sample_info set_recipe_outfile_in_sample_info };

    ### PREPROCESSING:

    ## Constants
    Readonly my $VEP_FORK_NUMBER => 4;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

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
    my $infile_path_prefix = $io{in}{file_path_prefix};
    my $infile_suffix      = $io{in}{file_suffix};
    my $infile_path        = $infile_path_prefix . $infile_suffix;

    my $consensus_analysis_type  = $parameter_href->{cache}{consensus_analysis_type};
    my $genome_reference_version = $file_info_href->{human_genome_reference_version};
    my $job_id_chain             = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number = $recipe_resource{core_number};

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
            }
        )
    );

    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my $outfile_path        = $outfile_path_prefix . $outfile_suffix;

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            max_cores_per_node   => $active_parameter_href->{max_cores_per_node},
            modifier_core_number => 1,
            recipe_core_number   => $core_number,
        }
    );

    ## Update memory depending on how many cores that are being used
    my $memory_allocation = update_memory_allocation(
        {
            node_ram_memory           => $active_parameter_href->{node_ram_memory},
            parallel_processes        => $core_number,
            process_memory_allocation => $recipe_resource{memory},
        }
    );

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $case_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $memory_allocation,
            process_time          => $recipe_resource{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            temp_directory        => $temp_directory,
        }
    );
    my $stderr_path = $recipe_info_path . $DOT . q{stderr.txt};

    # Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) = splitpath($stderr_path);

    ### SHELL:

    ## Reformat SV with no length as these will fail in the annotation with VEP
    _reformat_sv_with_no_length(
        {
            filehandle         => $filehandle,
            file_suffix        => $infile_suffix,
            infile_path_prefix => $infile_path_prefix,
        }
    );

    ## Get the vep synonyms file path for if required (grch38)
    my $vep_synonyms_file_path = create_vep_synonyms_file(
        {
            outfile_path => catfile( $outdir_path_prefix, q{synonyms.tsv} ),
            version      => $genome_reference_version,
        }
    );

    ## Varianteffectpredictor
    say {$filehandle} q{## Varianteffectpredictor};

    my $assembly_version =
      $file_info_href->{human_genome_reference_source} . $genome_reference_version;

    ## Get genome source and version to be compatible with VEP
    $assembly_version = _get_assembly_name( { assembly_version => $assembly_version, } );

    # VEP custom annotations
    my @custom_annotations = _get_custom_annotation_cmds(
        {
            vep_custom_annotation_href => $active_parameter_href->{vep_custom_annotation},
        }
    );

    ## VEP plugins
    my @plugins =
      _get_plugin_cmds( { vep_plugin_href => $active_parameter_href->{sv_vep_plugin}, } );

    ## VEP features
    my @vep_features_ref;

  FEATURE:
    foreach my $vep_feature ( @{ $active_parameter_href->{sv_vep_features} } ) {

        # Add VEP features to the output.
        push @vep_features_ref, $vep_feature;
    }

    my $vep_infile_path = $infile_path_prefix . $UNDERSCORE . q{fixedsvlength} . $infile_suffix;
    my $stderrfile_path = $recipe_file_path . $DOT . q{stderr.txt};
    my $stdoutfile_path = $recipe_file_path . $DOT . q{stdout.txt};
    variant_effect_predictor(
        {
            assembly           => $assembly_version,
            buffer_size        => $BUFFER_SIZE_SV,
            cache_directory    => $active_parameter_href->{vep_directory_cache},
            filehandle         => $filehandle,
            fork               => $VEP_FORK_NUMBER,
            infile_format      => substr( $infile_suffix, 1 ),
            infile_path        => $vep_infile_path,
            outfile_format     => substr( $outfile_suffix, 1 ),
            outfile_path       => $outfile_path,
            plugins_dir_path   => $active_parameter_href->{vep_plugins_dir_path},
            plugins_ref        => \@plugins,
            reference_path     => $active_parameter_href->{human_genome_reference},
            stderrfile_path    => $stderrfile_path,
            stdoutfile_path    => $stdoutfile_path,
            synonyms_file_path => $vep_synonyms_file_path,
            vep_features_ref   => \@vep_features_ref,
        }
    );
    say {$filehandle} $NEWLINE;

    close $filehandle;

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
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_case},
                job_id_chain                      => $job_id_chain,
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
                log                               => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_ids_ref     => \@{ $active_parameter_href->{sample_ids} },
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub analysis_vep_sv_wgs {

## Function : Varianteffectpredictor performs annotation of wgs SV variants.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $filehandle              => Filehandle to write to
##          : $file_info_href          => The file_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
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
        profile_base_command => {
            default     => q{sbatch},
            store       => \$profile_base_command,
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

    use MIP::Cluster qw{ get_core_number update_memory_allocation };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Vep qw{ variant_effect_predictor };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Sample_info qw{ set_recipe_metafile_in_sample_info set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constants
    Readonly my $VEP_FORK_NUMBER => 4;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

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
    my $infile_path_prefix = $io{in}{file_path_prefix};
    my $infile_suffix      = $io{in}{file_suffix};
    my $infile_path        = $infile_path_prefix . $infile_suffix;

    my $contigs_ref              = \@{ $file_info_href->{contigs} };
    my $consensus_analysis_type  = $parameter_href->{cache}{consensus_analysis_type};
    my $genome_reference_version = $file_info_href->{human_genome_reference_version};

    my $job_id_chain = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
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
            }
        )
    );

    my $outdir_path_prefix = $io{out}{dir_path_prefix};
    my @outfile_paths      = @{ $io{out}{file_paths} };
    my $outfile_suffix     = $io{out}{file_suffix};
    my %outfile_path       = %{ $io{out}{file_path_href} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle      = IO::Handle->new();
    my $xargsfilehandle = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    my $core_number = get_core_number(
        {
            max_cores_per_node   => $active_parameter_href->{max_cores_per_node},
            modifier_core_number => scalar @{ $file_info_href->{contigs_size_ordered} },
            recipe_core_number   => $recipe_resource{core_number},
        }
    );

    # Adjust for the number of forks vep forks
    my $parallel_processes = floor( $core_number / $VEP_FORK_NUMBER );

    ## Update memory depending on how many cores that are being used
    my $memory_allocation = update_memory_allocation(
        {
            node_ram_memory           => $active_parameter_href->{node_ram_memory},
            parallel_processes        => $parallel_processes,
            process_memory_allocation => $recipe_resource{memory},
        }
    );

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $case_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $memory_allocation,
            process_time          => $recipe_resource{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            temp_directory        => $temp_directory,
        }
    );
    my $stderr_path = $recipe_info_path . $DOT . q{stderr.txt};

    # Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) = splitpath($stderr_path);

    ### SHELL:
    ## Reformat SV with no length as these will fail in the annotation with VEP
    _reformat_sv_with_no_length(
        {
            filehandle         => $filehandle,
            file_suffix        => $infile_suffix,
            infile_path_prefix => $infile_path_prefix,
        }
    );

    ## Get the vep synonyms file path for if required (grch38)
    my $vep_synonyms_file_path = create_vep_synonyms_file(
        {
            outfile_path => catfile( $outdir_path_prefix, q{synonyms.tsv} ),
            version      => $genome_reference_version,
        }
    );

    ## Varianteffectpredictor
    say {$filehandle} q{## Varianteffectpredictor};

    my $assembly_version =
      $file_info_href->{human_genome_reference_source} . $genome_reference_version;

    ## Get genome source and version to be compatible with VEP
    $assembly_version = _get_assembly_name( { assembly_version => $assembly_version, } );

    my $vep_infile_path = $infile_path_prefix . $UNDERSCORE . q{fixedsvlength} . $infile_suffix;

    # VEP custom annotations
    my @custom_annotations = _get_custom_annotation_cmds(
        {
            vep_custom_annotation_href => $active_parameter_href->{vep_custom_annotation},
        }
    );

    ## VEP plugins
    my @plugins =
      _get_plugin_cmds( { vep_plugin_href => $active_parameter_href->{sv_vep_plugin}, } );

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $parallel_processes,
            filehandle         => $filehandle,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            xargsfilehandle    => $xargsfilehandle,
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
                assembly           => $assembly_version,
                buffer_size        => $BUFFER_SIZE_SV,
                cache_directory    => $active_parameter_href->{vep_directory_cache},
                distance           => $distance,
                filehandle         => $xargsfilehandle,
                fork               => $VEP_FORK_NUMBER,
                infile_format      => substr( $infile_suffix, 1 ),
                infile_path        => $vep_infile_path,
                outfile_format     => substr( $outfile_suffix, 1 ),
                outfile_path       => $outfile_path{$contig},
                plugins_dir_path   => $active_parameter_href->{vep_plugins_dir_path},
                plugins_ref        => \@plugins,
                regions_ref        => \@regions,
                reference_path     => $active_parameter_href->{human_genome_reference},
                stderrfile_path    => $stderrfile_path,
                stdoutfile_path    => $stdoutfile_path,
                synonyms_file_path => $vep_synonyms_file_path,
                vep_features_ref   => \@vep_features_ref,
            }
        );
        say {$xargsfilehandle} $NEWLINE;

        ## Filter out the MT annotations from the combined chr21 chrM call
        if ($mt_name) {

            say {$filehandle} q{## Filter out MT annotations};
            _subset_vcf(
                {
                    filehandle   => $filehandle,
                    infile_path  => $outfile_path{$contig},
                    outfile_path => $outfile_path{$contig},
                    regions_ref  => [qw{ chrM MT }],
                }
            );
        }
    }

    close $xargsfilehandle;
    close $filehandle;

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
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_case},
                job_id_chain                      => $job_id_chain,
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
                log                               => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_ids_ref     => \@{ $active_parameter_href->{sample_ids} },
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub analysis_vep {

## Function : Varianteffectpredictor performs effect predictions and annotation of variants.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $file_path               => File path
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $file_path;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;

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
        file_path   => { store => \$file_path, strict_type => 1, },
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
        profile_base_command => {
            default     => q{sbatch},
            store       => \$profile_base_command,
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

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::List qw{ get_splitted_lists };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Vep qw{ variant_effect_predictor };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constants
    Readonly my $VEP_FORK_NUMBER => 4;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

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

    my $consensus_analysis_type  = $parameter_href->{cache}{consensus_analysis_type};
    my $genome_reference_version = $file_info_href->{human_genome_reference_version};
    my $job_id_chain             = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
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

    my $outdir_path_prefix = $io{out}{dir_path_prefix};
    my $outfile_path       = $io{out}{file_path};
    my $outfile_suffix     = $io{out}{file_suffix};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe_resource{core_number},
            directory_id          => $case_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe_resource{memory},
            process_time          => $recipe_resource{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
        }
    );

    ### SHELL:

    ## Varianteffectpredictor
    say {$filehandle} q{## Varianteffectpredictor};

    my $assembly_version =
      $file_info_href->{human_genome_reference_source} . $genome_reference_version;

    ## Get genome source and version to be compatible with VEP
    $assembly_version = _get_assembly_name( { assembly_version => $assembly_version } );

    ## Get the vep synonyms file path for if required (grch38)
    my $vep_synonyms_file_path = create_vep_synonyms_file(
        {
            outfile_path => catfile( $outdir_path_prefix, q{synonyms.tsv} ),
            version      => $genome_reference_version,
        }
    );

    ## Get contigs
    my ( $mt_contig_ref, $contigs_ref ) = get_splitted_lists(
        {
            regexp   => qr/M/,
            list_ref => $file_info_href->{bam_contigs},
        }
    );

    ## Get plugins
    my @plugins =
      _get_plugin_cmds( { vep_plugin_href => $active_parameter_href->{vep_plugin}, } );

    # Get VEP custom annotations
    my @custom_annotations = _get_custom_annotation_cmds(
        {
            vep_custom_annotation_href => $active_parameter_href->{vep_custom_annotation},
        }
    );

    ## VEP features
    my ( @vep_features_ref, @vep_features_mt_ref );

  FEATURE:
    foreach my $vep_feature ( @{ $active_parameter_href->{vep_features} } ) {

        # Add VEP features to the output.
        push @vep_features_ref,    $vep_feature;
        push @vep_features_mt_ref, $vep_feature;

        # Special case for mitochondrial contig annotation
        if ( $vep_feature eq q{refseq} ) {

            push @vep_features_mt_ref, q{all_refseq};
        }
    }

    variant_effect_predictor(
        {
            assembly               => $assembly_version,
            buffer_size            => $BUFFER_SIZE,
            cache_directory        => $active_parameter_href->{vep_directory_cache},
            custom_annotations_ref => \@custom_annotations,
            distance               => $ANNOTATION_DISTANCE,
            filehandle             => $filehandle,
            fork                   => $VEP_FORK_NUMBER,
            infile_format          => substr( $infile_suffix, 1 ),
            infile_path            => $infile_path,
            outfile_format         => substr( $outfile_suffix, 1 ),
            outfile_path           => $outfile_path,
            plugins_dir_path       => $active_parameter_href->{vep_plugins_dir_path},
            plugins_ref            => \@plugins,
            reference_path         => $active_parameter_href->{human_genome_reference},
            regions_ref            => $contigs_ref,
            synonyms_file_path     => $vep_synonyms_file_path,
            vep_features_ref       => \@vep_features_ref,
        }
    );
    say {$filehandle} $NEWLINE;

    ## VEP for MT
    variant_effect_predictor(
        {
            assembly               => $assembly_version,
            buffer_size            => $BUFFER_SIZE,
            cache_directory        => $active_parameter_href->{vep_directory_cache},
            custom_annotations_ref => \@custom_annotations,
            distance               => $ANNOTATION_DISTANCE_MT,
            filehandle             => $filehandle,
            fork                   => $VEP_FORK_NUMBER,
            infile_format          => substr( $infile_suffix, 1 ),
            infile_path            => $infile_path,
            no_headers             => 1,
            outfile_format         => substr( $outfile_suffix, 1 ),
            outfile_path           => q{STDOUT},
            plugins_dir_path       => $active_parameter_href->{vep_plugins_dir_path},
            plugins_ref            => \@plugins,
            reference_path         => $active_parameter_href->{human_genome_reference},
            regions_ref            => $mt_contig_ref,
            synonyms_file_path     => $vep_synonyms_file_path,
            stdoutfile_path_append => $outfile_path,
            vep_features_ref       => \@vep_features_mt_ref,
        }
    );
    say {$filehandle} $NEWLINE;

    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        if ( $consensus_analysis_type eq q{wts} ) {

            set_file_path_to_store(
                {
                    format           => q{vcf},
                    id               => $case_id,
                    path             => $outfile_path,
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                }
            );
        }

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_case},
                job_id_chain                      => $job_id_chain,
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
                log                               => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_ids_ref     => \@{ $active_parameter_href->{sample_ids} },
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
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
            default     => $arg_href->{vep_custom_annotation_href} ||= undef,
            required    => 1,
            store       => \$vep_custom_annotation_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( not defined $vep_custom_annotation_href );

    my @custom_annotations;
    my @order_custom_options =
      qw{ path key file_type annotation_type force_report_coordinates vcf_fields};

  ANNOTATION:
    foreach my $annotation_href ( values %{$vep_custom_annotation_href} ) {

        ## Remove all undef elements and then join
        my $cmd = join $COMMA, grep { defined } @{$annotation_href}{@order_custom_options};
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
    if ( $assembly_version =~ /grch(\d+)/xsm ) {

        my $version_number = $1;

        if ( $version_number > $VEP_CACHE_GENOME_BUILD_NAME_SWITCH ) {

            $assembly_version = q{GRCh} . $version_number;
        }
    }
    return $assembly_version;
}

sub _get_plugin_cmds {

## Function : Build plugin command per plugin
## Returns  : @plugins
## Arguments: $vep_plugin_href => Plugin info {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vep_plugin_href;

    my $tmpl = {
        vep_plugin_href => {
            default     => $arg_href->{vep_plugin_href} ||= undef,
            required    => 1,
            store       => \$vep_plugin_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @plugins;

  PLUGIN:
    while ( my ( $plugin_name, $plugin_href ) = each %{$vep_plugin_href} ) {

        my $cmd = $plugin_name;

        if ( exists $plugin_href->{parameters} ) {

            $cmd .= $COMMA . join $COMMA, @{ $plugin_href->{parameters} };
        }
        push @plugins, $cmd;
    }
    return @plugins;
}

sub _reformat_sv_with_no_length {

## Function : Reformat SV with no length as these will fail in the annotation with VEP
## Returns  :
## Arguments: $filehandle         => Filehandle to write to
##          : $file_suffix        => File suffix
##          : $infile_path_prefix => Infile path prefix

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $file_suffix;
    my $infile_path_prefix;

    my $tmpl = {
        filehandle  => { required => 1, store => \$filehandle, },
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

    use MIP::Environment::Executable qw{ get_executable_base_command };

    ## Fix SV with no length as these will fail in the annotation with VEP
    my @commands = ( get_executable_base_command( { base_command => q{perl}, } ), );

    # Execute perl
    my $perl_fix_sv_nolengths = join $SPACE, @commands;

    # Set up perl
    $perl_fix_sv_nolengths .= q? -nae '?;

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

    print {$filehandle} $perl_fix_sv_nolengths . $SPACE;
    print {$filehandle} $infile_path_prefix . $file_suffix . $SPACE;
    say   {$filehandle} q{>}
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
## Arguments: $filehandle   => Filehandle
##          : $infile_path  => Path to infile
##          : $outfile_path => Path to outfile
##          : $regions_ref  => Array with regions to be included in the subset {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $regions_ref;

    my $tmpl = {
        filehandle => {
            required => 1,
            store    => \$filehandle,
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

    use MIP::Program::Htslib qw{ htslib_bgzip };
    use MIP::Program::Bcftools qw{ bcftools_index bcftools_view };

    ## Prepare for bcftools_view
    htslib_bgzip(
        {
            filehandle  => $filehandle,
            force       => 1,
            infile_path => $infile_path,
        }
    );
    print {$filehandle} $NEWLINE;
    bcftools_index(
        {
            filehandle  => $filehandle,
            infile_path => $infile_path . q{.gz},
        }
    );
    print {$filehandle} $NEWLINE;
    bcftools_view(
        {
            regions_ref  => $regions_ref,
            filehandle   => $filehandle,
            infile_path  => $infile_path . q{.gz},
            outfile_path => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    return;
}

1;
