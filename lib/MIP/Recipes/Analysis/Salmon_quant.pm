package MIP::Recipes::Analysis::Salmon_quant;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
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
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_salmon_quant };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_salmon_quant {

## Function : Transcript quantification using salmon quant
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name            => Program name
##          : $sample_id               => Sample id
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
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
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
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_parameters get_recipe_attributes };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Program::Variantcalling::Salmon qw{ salmon_quant };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::QC::Sample_info qw{ set_recipe_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            file_info_href => $file_info_href,
            id             => $sample_id,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my @infile_paths      = @{ $io{in}{file_paths} };
    my @temp_infile_paths = @{ $io{temp}{file_paths} };
    my $recipe_mode       = $active_parameter_href->{$recipe_name};
    my $referencefile_dir_path =
        $active_parameter_href->{salmon_quant_reference_genome}
      . $file_info_href->{salmon_quant_reference_genome}[0];
    my $job_id_chain = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => \@{ $infile_lane_prefix_href->{$sample_id} },
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
                temp_directory         => $temp_directory,
            }
        )
    );
    my $outdir_path                = $io{out}{dir_path};
    my @outfile_name_prefixes      = @{ $io{out}{file_name_prefixes} };
    my @outfile_paths              = @{ $io{out}{file_paths} };
    my @temp_outfile_path_prefixes = @{ $io{temp}{file_path_prefixes} };
    my @temp_outfile_paths         = @{ $io{temp}{file_paths} };

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    # Too avoid adjusting infile_index in submitting to jobs
    my $paired_end_tracker = 0;

    ## Perform per single-end or read pair
  INFILE_PREFIX:
    while ( my ( $infile_index, $infile_prefix ) =
        each @{ $infile_lane_prefix_href->{$sample_id} } )
    {

        ## Assign file features
        my $file_path           = $temp_outfile_paths[$infile_index];
        my $file_path_prefix    = $temp_outfile_path_prefixes[$infile_index];
        my $outfile_name_prefix = $outfile_name_prefixes[$infile_index];
        my $outfile_path        = $outfile_paths[$infile_index];

        # Collect paired-end or single-end sequence run mode
        my $sequence_run_mode =
          $sample_info_href->{sample}{$sample_id}{file}{$infile_prefix}
          {sequence_run_type};

        ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
        my ( $recipe_file_path, $recipe_info_path ) = setup_script(
            {
                active_parameter_href           => $active_parameter_href,
                core_number                     => $core_number,
                directory_id                    => $sample_id,
                FILEHANDLE                      => $FILEHANDLE,
                job_id_href                     => $job_id_href,
                log                             => $log,
                recipe_directory                => $recipe_name,
                recipe_name                     => $recipe_name,
                process_time                    => $time,
                source_environment_commands_ref => \@source_environment_cmds,
                temp_directory                  => $temp_directory,
            }
        );

        ### SHELL

        ## Copies file to temporary directory.
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};

        # Read 1
        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $infile_paths[$paired_end_tracker],
                outfile_path => $temp_directory,
            }
        );

        # If second read direction is present
        if ( $sequence_run_mode eq q{paired-end} ) {

            # Read 2
            migrate_file(
                {
                    FILEHANDLE   => $FILEHANDLE,
                    infile_path  => $infile_paths[ $paired_end_tracker + 1 ],
                    outfile_path => $temp_directory,
                }
            );
        }
        say {$FILEHANDLE} q{wait}, $NEWLINE;

        ## Salmon quant
        say {$FILEHANDLE} q{## Quantifying transcripts using } . $recipe_name;

        ### Get parameters
        ## Infile(s)
        my @fastq_files = $temp_infile_paths[$paired_end_tracker];

        # If second read direction is present
        if ( $sequence_run_mode eq q{paired-end} ) {

            # Increment to collect correct read 2 from %infile
            $paired_end_tracker++;
            push @fastq_files, $temp_infile_paths[$paired_end_tracker];
        }

        # If second read direction is present
        if ( $sequence_run_mode ne q{paired-end} ) {

            salmon_quant(
                {
                    FILEHANDLE        => $FILEHANDLE,
                    gc_bias           => 1,
                    index_path        => $referencefile_dir_path,
                    outfile_path      => $file_path_prefix,
                    read_1_fastq_path => $fastq_files[0],
                },
            );
            say {$FILEHANDLE} $NEWLINE;
        }
        else {
            salmon_quant(
                {
                    FILEHANDLE        => $FILEHANDLE,
                    gc_bias           => 1,
                    index_path        => $referencefile_dir_path,
                    outfile_path      => $file_path_prefix,
                    read_1_fastq_path => $fastq_files[0],
                    read_2_fastq_path => $fastq_files[1],
                },
            );
            say {$FILEHANDLE} $NEWLINE;
        }

        ## Increment paired end tracker
        $paired_end_tracker++;

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $file_path_prefix . $ASTERISK,
                outfile_path => $outdir_path,
            }
        );
        say {$FILEHANDLE} q{wait}, $NEWLINE;

        ## Close FILEHANDLES
        close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

        if ( $recipe_mode == 1 ) {

            ## Collect QC metadata info for later use
            set_recipe_outfile_to_sample_info(
                {
                    infile           => $outfile_name_prefix,
                    path             => $outfile_path,
                    recipe_name      => $recipe_name,
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );

            submit_recipe(
                {
                    dependency_method       => q{sample_to_sample_parallel},
                    case_id                 => $case_id,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    log                     => $log,
                    job_id_chain            => $job_id_chain,
                    recipe_file_path        => $recipe_file_path,
                    recipe_files_tracker    => $infile_index,
                    sample_id               => $sample_id,
                    submission_profile => $active_parameter_href->{submission_profile},
                }
            );
        }
    }
    return;
}

1;
