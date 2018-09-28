package MIP::Recipes::Analysis::Tiddit;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile devnull };
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
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_tiddit analysis_tiddit_coverage };

}

## Constants
Readonly my $UNDERSCORE => q{_};
Readonly my $SPACE      => q{ };
Readonly my $ASTERISK   => q{*};
Readonly my $NEWLINE    => qq{\n};
Readonly my $AMPERSAND  => q{&};

sub analysis_tiddit {

## Function : Call structural variants using Tiddit 1.0.2
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => The file info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner directory used in the analysis
##          : $parameter_href          => Parameter hash {REF}
##          : $reference_dir           => MIP reference directory
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
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
    my $reference_dir;
    my $temp_directory;

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
        reference_dir => {
            default     => $arg_href->{active_parameter_href}{reference_dir},
            store       => \$reference_dir,
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

    use MIP::Cluster qw{ get_core_number };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_module_parameters get_program_attributes };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Processmanagement::Processes qw{ print_wait };
    use MIP::Program::Variantcalling::Svdb qw{ svdb_merge };
    use MIP::Program::Variantcalling::Tiddit qw{ tiddit_sv };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    my $job_id_chain = get_program_attributes(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
            attribute      => q{chain},
        }
    );
    my $max_cores_per_node = $active_parameter_href->{max_cores_per_node};
    my $modifier_core_number =
      scalar( @{ $active_parameter_href->{sample_ids} } );
    my $program_mode = $active_parameter_href->{$program_name};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Set and get the io files per chain, id and stream
    my %io = parse_io_outfiles(
        {
            chain_id               => $job_id_chain,
            id                     => $family_id,
            file_info_href         => $file_info_href,
            file_name_prefixes_ref => [$family_id],
            outdata_dir            => $active_parameter_href->{outdata_dir},
            parameter_href         => $parameter_href,
            program_name           => $program_name,
            temp_directory         => $temp_directory,
        }
    );

    my $outdir_path_prefix       = $io{out}{dir_path_prefix};
    my $outfile_path_prefix      = $io{out}{file_path_prefix};
    my $outfile_suffix           = $io{out}{file_suffix};
    my $outfile_path             = $outfile_path_prefix . $outfile_suffix;
    my $temp_outfile_path_prefix = $io{temp}{file_path_prefix};
    my $temp_outfile_suffix      = $io{temp}{file_suffix};
    my $temp_outfile_path = $temp_outfile_path_prefix . $temp_outfile_suffix;

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    $core_number = get_core_number(
        {
            max_cores_per_node   => $max_cores_per_node,
            modifier_core_number => $modifier_core_number,
            module_core_number   => $core_number,
        }
    );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $family_id,
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

    ### SHELL:

    my %tiddit_sample_file_info;
    my $process_batches_count = 1;

    ## Collect infiles for all sample_ids to enable migration to temporary directory
    while ( my ( $sample_id_index, $sample_id ) =
        each @{ $active_parameter_href->{sample_ids} } )
    {

        ## Get the io infiles per chain and id
        my %sample_io = get_io_files(
            {
                id             => $sample_id,
                file_info_href => $file_info_href,
                parameter_href => $parameter_href,
                program_name   => $program_name,
                stream         => q{in},
                temp_directory => $temp_directory,
            }
        );
        my $infile_path_prefix = $sample_io{in}{file_path_prefix};
        my $infile_suffix      = $sample_io{in}{file_suffix};
        my $infile_path =
          $infile_path_prefix . substr( $infile_suffix, 0, 2 ) . $ASTERISK;
        my $temp_infile_path_prefix = $sample_io{temp}{file_path_prefix};
        my $temp_infile_path        = $temp_infile_path_prefix . $infile_suffix;

        $tiddit_sample_file_info{$sample_id}{in} = $temp_infile_path;
        $tiddit_sample_file_info{$sample_id}{out} =
          $temp_infile_path_prefix . $UNDERSCORE . $sample_id;

        $process_batches_count = print_wait(
            {
                FILEHANDLE            => $FILEHANDLE,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                process_counter       => $sample_id_index,
            }
        );

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $infile_path,
                outfile_path => $temp_directory,
            }
        );
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    # Restart counter
    $process_batches_count = 1;

    ## Collect infiles for all sample_ids
    while ( my ( $sample_id_index, $sample_id ) =
        each @{ $active_parameter_href->{sample_ids} } )
    {

        $process_batches_count = print_wait(
            {
                FILEHANDLE            => $FILEHANDLE,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                process_counter       => $sample_id_index,
            }
        );

        ## Tiddit
        tiddit_sv(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $tiddit_sample_file_info{$sample_id}{in},
                minimum_number_supporting_pairs => $active_parameter_href
                  ->{tiddit_minimum_number_supporting_pairs},
                outfile_path_prefix =>
                  $tiddit_sample_file_info{$sample_id}{out},
                referencefile_path =>
                  $active_parameter_href->{human_genome_reference},
            }
        );
        say {$FILEHANDLE} $AMPERSAND . $SPACE . $NEWLINE;
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Get parameters
    ## Tiddit sample outfiles needs to be lexiographically sorted for svdb merge
    my @svdb_temp_infile_paths =
      map { $tiddit_sample_file_info{$_}{out} . $outfile_suffix } @{ $active_parameter_href->{sample_ids} };

    svdb_merge(
        {
            FILEHANDLE       => $FILEHANDLE,
            infile_paths_ref => \@svdb_temp_infile_paths,
            notag            => 1,
            stdoutfile_path  => $temp_outfile_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $temp_outfile_path . $ASTERISK,
            outfile_path => $outdir_path_prefix,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        add_program_outfile_to_sample_info(
            {
                path             => $outfile_path,
                program_name     => q{tiddit},
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

sub analysis_tiddit_coverage {

## Function : Calculate coverage using Tiddit 1.0.2
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $bin_size                => Bin size
##          : $family_id               => Family id
##          : $file_info_href          => The file info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner directory used in the analysis
##          : $outfamily_directory     => Out family directory
##          : $parameter_href          => Parameter hash {REF}
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
            store       => \$bin_size,
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

    use MIP::Cluster qw{ get_core_number };
    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Processmanagement::Processes qw{ print_wait };
    use MIP::Program::Variantcalling::Svdb qw{ svdb_merge };
    use MIP::Program::Variantcalling::Tiddit qw{ tiddit_coverage };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $program_mode = $active_parameter_href->{$program_name};

    ## Alias
    my $job_id_chain       = $parameter_href->{$program_name}{chain};
    my $max_cores_per_node = $active_parameter_href->{max_cores_per_node};
    my $modifier_core_number =
      scalar( @{ $active_parameter_href->{sample_ids} } );
    my $program_outdirectory_name =
      $parameter_href->{$program_name}{outdir_name};
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

    $core_number = get_core_number(
        {
            max_cores_per_node   => $max_cores_per_node,
            modifier_core_number => $modifier_core_number,
            module_core_number   => $core_number,
        }
    );

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

    # Used downstream
    $parameter_href->{$program_name}{indirectory} = $outfamily_directory;

    ## Assign file_tags
    my %file_path_prefix;
    my $outfile_tag =
      $file_info_href->{$family_id}{$program_name}{file_tag};
    my $outfile_prefix = $family_id . $outfile_tag;
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

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
            file_suffix =>
              $parameter_href->{$program_name}{coverage_file_suffix},
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

    my $process_batches_count = 1;

    ## Collect infiles for all sample_ids to enable migration to temporary directory
    while ( my ( $sample_id_index, $sample_id ) =
        each @{ $active_parameter_href->{sample_ids} } )
    {

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
          $file_info_href->{$sample_id}{gatk_baserecalibration}{file_tag};
        my $infile_prefix         = $merged_infile_prefix . $infile_tag;
        my $sample_outfile_prefix = $merged_infile_prefix . $outfile_tag;

        #q{.bam} -> ".b*" for getting index as well
        my $infile_path = catfile( $insample_directory,
            $infile_prefix . substr( $infile_suffix, 0, 2 ) . $ASTERISK );

        $file_path_prefix{$sample_id}{in} =
          catfile( $temp_directory, $infile_prefix );
        $file_path_prefix{$sample_id}{out} =
          catfile( $temp_directory, $sample_outfile_prefix );

        $process_batches_count = print_wait(
            {
                FILEHANDLE            => $FILEHANDLE,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                process_counter       => $sample_id_index,
            }
        );

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $infile_path,
                outfile_path => $temp_directory,
            }
        );
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    # Restart counter
    $process_batches_count = 1;

    ## Calculate coverage and create .cgh outfiles for all sample_ids
    while ( my ( $sample_id_index, $sample_id ) =
        each @{ $active_parameter_href->{sample_ids} } )
    {

        $process_batches_count = print_wait(
            {
                FILEHANDLE            => $FILEHANDLE,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                process_counter       => $sample_id_index,
            }
        );

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
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $outfile_path_prefix . $outfile_suffix . $ASTERISK,
            outfile_path => $outfamily_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        add_program_outfile_to_sample_info(
            {
                path => catfile(
                    $outfamily_directory, $outfile_prefix . $outfile_suffix
                ),
                program_name     => q{tiddit_coverage},
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
