package MIP::Recipes::Analysis::Manta;

use 5.026;
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
    our $VERSION = 1.08;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_manta };

}

## Constants
Readonly my $UNDERSCORE => q{_};
Readonly my $NEWLINE    => qq{\n};
Readonly my $ASTERISK   => q{*};

sub analysis_manta {

## Function : Joint analysis of structural variation
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => The file_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $referencefile_path      => Path to reference file
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

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
    my $profile_base_command;
    my $referencefile_path;
    my $temp_directory;

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
        referencefile_path => {
            default     => $arg_href->{active_parameter_href}{human_genome_reference},
            store       => \$referencefile_path,
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

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_parameters get_recipe_attributes };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ print_wait submit_recipe };
    use MIP::Program::Compression::Gzip qw{ gzip };
    use MIP::Program::Variantcalling::Manta qw{ manta_config manta_workflow };
    use MIP::QC::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
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

    ## Set and get the io files per chain, id and stream
    my %io = parse_io_outfiles(
        {
            chain_id               => $job_id_chain,
            id                     => $case_id,
            file_info_href         => $file_info_href,
            file_name_prefixes_ref => [$case_id],
            outdata_dir            => $active_parameter_href->{outdata_dir},
            parameter_href         => $parameter_href,
            recipe_name            => $recipe_name,
            temp_directory         => $temp_directory,
        }
    );

    my $outdir_path_prefix       = $io{out}{dir_path_prefix};
    my $outfile_path_prefix      = $io{out}{file_path_prefix};
    my $outfile_suffix           = $io{out}{file_suffix};
    my $outfile_path             = $outfile_path_prefix . $outfile_suffix;
    my $temp_outfile_path_prefix = $io{temp}{file_path_prefix};
    my $temp_outfile_suffix      = $io{temp}{file_suffix};
    my $temp_outfile_path        = $temp_outfile_path_prefix . $temp_outfile_suffix;

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ($recipe_file_path) = setup_script(
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

    ### SHELL:

    ## Collect infiles for all sample_ids to enable migration to temporary directory
    my @manta_temp_infile_paths;
    my $process_batches_count = 1;
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

        ## Get the io infiles per chain and id
        my %sample_io = get_io_files(
            {
                id             => $sample_id,
                file_info_href => $file_info_href,
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
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

        ## Store temp infile path for each sample_id
        push @manta_temp_infile_paths, $temp_infile_path;

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

    ## Manta
    say {$FILEHANDLE} q{## Manta};

    ## Get parameters
    my $exome_analysis;
    if ( $consensus_analysis_type ne q{wgs} ) {

        $exome_analysis = 1;
    }

    manta_config(
        {
            exome_analysis     => $exome_analysis,
            FILEHANDLE         => $FILEHANDLE,
            infile_paths_ref   => \@manta_temp_infile_paths,
            outdirectory_path  => $temp_directory,
            referencefile_path => $referencefile_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Manta workflow};
    manta_workflow(
        {
            FILEHANDLE        => $FILEHANDLE,
            mode              => q{local},
            outdirectory_path => $temp_directory,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    my $manta_temp_outfile_path =
      catfile( $temp_directory, qw{ results variants diploidSV.vcf.gz } );

    ## Perl wrapper for writing gzip recipe to $FILEHANDLE
    gzip(
        {
            decompress   => 1,
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $manta_temp_outfile_path,
            outfile_path => $temp_outfile_path,
            stdout       => 1,
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

    if ( $recipe_mode == 1 ) {

        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => q{manta},
                sample_info_href => $sample_info_href,
            }
        );
        submit_recipe(
            {
                base_command            => $profile_base_command,
                case_id                 => $case_id,
                dependency_method       => q{sample_to_case},
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
    return 1;
}

1;
