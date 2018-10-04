package MIP::Recipes::Analysis::Expansionhunter;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
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
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_expansionhunter };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $AMPERSAND  => q{&};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_expansionhunter {

## Function : Call expansions of Short Tandem Repeats (STR) using Expansion Hunter
## Returns  :
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => The file_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
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
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ print_wait };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_rename_vcf_samples bcftools_view };
    use MIP::Program::Variantcalling::Expansionhunter qw{ expansionhunter };
    use MIP::Program::Variantcalling::Svdb qw{ svdb_merge };
    use MIP::Program::Variantcalling::Vt qw{ vt_decompose };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    my $max_cores_per_node = $active_parameter_href->{max_cores_per_node};
    my $modifier_core_number =
      scalar( @{ $active_parameter_href->{sample_ids} } );
    my $human_genome_reference =
      $arg_href->{active_parameter_href}{human_genome_reference};
    my $job_id_chain = get_program_attributes(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
            attribute      => q{chain},
        }
    );
    my $program_mode = $active_parameter_href->{$program_name};
    my $repeat_specs_dir_path =
      $active_parameter_href->{expansionhunter_repeat_specs_dir};
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
    my $outfile_suffix           = $io{out}{file_constant_suffix};
    my $outfile_path             = $outfile_path_prefix . $outfile_suffix;
    my $temp_outfile_path_prefix = $io{temp}{file_path_prefix};
    my $temp_outfile_suffix      = $io{temp}{file_suffix};
    my $temp_outfile_path = $temp_outfile_path_prefix . $temp_outfile_suffix;
    my $temp_file_suffix  = $DOT . q{vcf};

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

    my %exphun_sample_file_info;
    my $process_batches_count = 1;

    ## Collect infiles for all sample_ids to enable migration to temporary directory
  SAMPLE_ID:
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

        $exphun_sample_file_info{$sample_id}{in}  = $temp_infile_path;
        $exphun_sample_file_info{$sample_id}{out} = $temp_infile_path_prefix;

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

    ## Rename the bam file index file so that Expansion Hunter can find it
    say {$FILEHANDLE} q{## Rename index file};
  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        gnu_mv(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $exphun_sample_file_info{$sample_id}{out}
                  . q{.bai},
                outfile_path => $exphun_sample_file_info{$sample_id}{in}
                  . q{.bai},
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Run Expansion Hunter
    say {$FILEHANDLE} q{## Run ExpansionHunter};

    # Restart counter
    $process_batches_count = 1;

    ## Expansion hunter calling per sample id
  SAMPLE_ID:
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

        my $sample_sex = $sample_info_href->{sample}{$sample_id}{sex};
        expansionhunter(
            {
                FILEHANDLE        => $FILEHANDLE,
                infile_path       => $exphun_sample_file_info{$sample_id}{in},
                json_outfile_path => $exphun_sample_file_info{$sample_id}{out}
                  . $DOT . q{json},
                log_outfile_path => $exphun_sample_file_info{$sample_id}{out}
                  . $DOT . q{log},
                reference_genome_path => $human_genome_reference,
                repeat_specs_dir_path => $repeat_specs_dir_path,
                sex                   => $sample_sex,
                vcf_outfile_path => $exphun_sample_file_info{$sample_id}{out}
                  . $temp_file_suffix,
            }
        );
        say {$FILEHANDLE} $AMPERSAND, $NEWLINE;
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Get parameters
    ## Expansionhunter sample infiles needs to be lexiographically sorted for svdb merge
    my @svdb_temp_infile_paths =
      map { $exphun_sample_file_info{$_}{out} . $temp_file_suffix }
      @{ $active_parameter_href->{sample_ids} };
    my $svdb_temp_outfile_path =
        $temp_outfile_path_prefix
      . $UNDERSCORE
      . q{svdbmerge}
      . $temp_file_suffix;

    svdb_merge(
        {
            FILEHANDLE       => $FILEHANDLE,
            infile_paths_ref => \@svdb_temp_infile_paths,
            notag            => 1,
            stdoutfile_path  => $svdb_temp_outfile_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Split multiallelic variants
    say {$FILEHANDLE} q{## Split multiallelic variants};
    my $vt_temp_outfile_path =
        $temp_outfile_path_prefix
      . $UNDERSCORE
      . q{svdbmerg_vt}
      . $temp_file_suffix;
    vt_decompose(
        {
            FILEHANDLE          => $FILEHANDLE,
            infile_path         => $svdb_temp_outfile_path,
            outfile_path        => $vt_temp_outfile_path,
            smart_decomposition => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Adding sample id instead of file prefix};
    bcftools_rename_vcf_samples(
        {
            FILEHANDLE          => $FILEHANDLE,
            index               => 1,
            index_type          => q{csi},
            infile              => $vt_temp_outfile_path,
            outfile_path_prefix => $outfile_path_prefix,
            output_type         => q{z},
            temp_directory      => $temp_directory,
            sample_ids_ref      => \@{ $active_parameter_href->{sample_ids} },
        }
    );

    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => q{expansionhunter},
                path             => $outfile_path,
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
