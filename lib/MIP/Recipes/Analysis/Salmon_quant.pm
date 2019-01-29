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
use List::MoreUtils qw{ uniq };
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.07;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_salmon_quant };

}

## Constants
Readonly my $NEWLINE    => qq{\n};

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
    use MIP::IO::Files qw{ migrate_files };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Program::Variantcalling::Salmon qw{ salmon_quant };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::QC::Sample_info qw{ set_recipe_outfile_in_sample_info };
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
    my @infile_names      = @{ $io{in}{file_names} };
    my $indir_path        = $io{in}{dir_path};
    my @temp_infile_paths = @{ $io{temp}{file_paths} };
    my $recipe_mode       = $active_parameter_href->{$recipe_name};
    my $referencefile_dir_path =
        $active_parameter_href->{salmon_quant_reference_genome}
      . $file_info_href->{salmon_quant_reference_genome}[0];
    my %rec_atr = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $job_id_chain = $rec_atr{chain};
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set outfile
    my $recipe_dir =
      catdir( $active_parameter_href->{outdata_dir}, $sample_id, $recipe_name );
    my $file_path = catfile( $recipe_dir, $rec_atr{file_tag} . $rec_atr{outfile_suffix} );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id       => $job_id_chain,
                id             => $sample_id,
                file_info_href => $file_info_href,
                file_paths_ref => [$file_path],
                outdata_dir    => $active_parameter_href->{outdata_dir},
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
                temp_directory => $temp_directory,
            }
        )
    );
    my $outdir_path  = $io{out}{dir_path};
    my $outfile_name = ${ $io{out}{file_names} }[0];
    my $outfile_path = ${ $io{out}{file_paths} }[0];

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Perform per single-end or read pair
    my @infile_prefixes = @{ $infile_lane_prefix_href->{$sample_id} };

    ## Collect paired-end or single-end sequence run mode
    my @sequence_run_modes;
  INFILE_PREFIX:
    foreach my $infile_prefix (@infile_prefixes) {
        push @sequence_run_modes,
          $sample_info_href->{sample}{$sample_id}{file}{$infile_prefix}
          {sequence_run_type};
    }
    @sequence_run_modes = uniq @sequence_run_modes;

    # Fail on mixed sequencing modes - THIS CHECK WILL BE MOVED IN A LATER PR
    if ( scalar @sequence_run_modes > 1 ) {
        $log->fatal(
            q{Salmon quant cannot operate on mixed single-end and paired-end fastq files}
        );
        return;
    }

    my $sequence_run_mode = $sequence_run_modes[0];

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
    migrate_files(
        {
            core_number  => $core_number,
            FILEHANDLE   => $FILEHANDLE,
            infiles_ref  => \@infile_names,
            indirectory  => $indir_path,
            outfile_path => $temp_directory,
        }
    );

    ## Salmon quant
    say {$FILEHANDLE} q{## Quantifying transcripts using } . $recipe_name;

    ## For paired_end, split first and second reads
    if ( $sequence_run_mode eq q{paired-end} ) {
        my @read_1_fastq_paths =
          @temp_infile_paths[ grep { !( $_ % 2 ) } 0 .. $#temp_infile_paths ];
        my @read_2_fastq_paths =
          @temp_infile_paths[ grep { ( $_ % 2 ) } 0 .. $#temp_infile_paths ];

        salmon_quant(
            {
                FILEHANDLE             => $FILEHANDLE,
                gc_bias                => 1,
                index_path             => $referencefile_dir_path,
                outdir_path            => $outdir_path,
                read_1_fastq_paths_ref => \@read_1_fastq_paths,
                read_2_fastq_paths_ref => \@read_2_fastq_paths,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }
    ## For single end
    else {
        salmon_quant(
            {
                FILEHANDLE        => $FILEHANDLE,
                gc_bias           => 1,
                index_path        => $referencefile_dir_path,
                outdir_path       => $outdir_path,
                read_1_fastq_path => \@temp_infile_paths,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Close FILEHANDLES
    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => $outfile_name,
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                dependency_method       => q{sample_to_sample},
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                job_id_chain            => $job_id_chain,
                recipe_file_path        => $recipe_file_path,
                sample_id               => $sample_id,
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return;
}

1;
