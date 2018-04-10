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
    our $VERSION = q{1.0.0};

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_expansionhunter };

}

## Constants
Readonly my $ASTERISK => q{*};
Readonly my $NEWLINE  => qq{\n};

sub analysis_expansionhunter {

## Function : Call expansions of STR using Expansion Hunter
## Returns  :
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => The file_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $insample_directory      => In sample directory
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outsample_directory     => Out sample directory
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $reference_dir           => MIP reference directory
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $sample_id               => Sample id
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
    my $reference_dir;
    my $temp_directory;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $FILEHANDLE;
    my $infile_lane_prefix_href;
    my $insample_directory;
    my $job_id_href;
    my $outsample_directory;
    my $parameter_href;
    my $program_name;
    my $sample_id;
    my $sample_info_href;

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
        insample_directory => {
            defined     => 1,
            required    => 1,
            store       => \$insample_directory,
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
        outsample_directory => {
            defined     => 1,
            required    => 1,
            store       => \$outsample_directory,
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
        sample_id => {
            default     => 1,
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

    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::Program::Variantcalling::Expansionhunter qw{ expansionhunter };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };

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
    my $program_outdirectory_name =
      $parameter_href->{$mip_program_name}{outdir_name};
    my $human_genome_reference_ref =
      $arg_href->{active_parameter_href}{human_genome_reference};

    ## Filehandles
    # Create anonymous filehandle
    $FILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $sample_id,
            program_name          => $program_name,
            program_directory =>
              catfile( $outaligner_dir, $program_outdirectory_name ),
            core_number                     => $core_number,
            process_time                    => $time,
            source_environment_commands_ref => [$source_environment_cmd],
            temp_directory                  => $temp_directory,
        }
    );

    #Used downstream
    $parameter_href->{$mip_program_name}{$sample_id}{indirectory} =
      $outsample_directory;

    ## Files
    my $infile = $file_info_href->{$sample_id}{merged_infile};
    my $infile_tag =
      $file_info_href->{$sample_id}{pgatk_baserecalibration}{file_tag};
    my $outfile_tag =
      $file_info_href->{$sample_id}{$mip_program_name}{file_tag};
    my $infile_prefix  = $infile . $infile_tag;
    my $outfile_prefix = $infile . $outfile_tag;

    ## Paths
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ### Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
            jobid_chain    => $parameter_href->{pgatk_baserecalibration}{chain},
        }
    );

    ## Set file suffix for next module within jobid chain
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
    my $infile_path = catfile( $insample_directory,
        $infile_prefix . substr( $infile_suffix, 0, 2 ) . $ASTERISK );
    migrate_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $infile_path,
            outfile_path => $temp_directory,
        }
    );
    say {$FILEHANDLE} q{wait} . $NEWLINE;

    ## Rename the bam file index file so that Expansion Hunter can find it
    say {$FILEHANDLE} q{## Rename index file};
    gnu_mv(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $file_path_prefix . q{.bai},
            outfile_path => $file_path_prefix . q{.bam.bai},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Get path to repeat_specs directory
    my $reference_genome_path =
      $active_parameter_href->{human_genome_reference};
    my $repeat_specs_dir_path =
      $active_parameter_href->{expansionhunter_repeat_specs_dir};

    ## Try to get default directory if varible is unset
    if ( not $repeat_specs_dir_path ) {
        $repeat_specs_dir_path = _get_default_repeat_specs_dir_path(
            {
                log                   => $log,
                reference_genome_path => $reference_genome_path,
            }
        );
    }

    if ( not -d $repeat_specs_dir_path ) {
        $log->fatal(
q{Can't find repeat specification directory for Expansion Hunter. Please set "expansionhunter_repeat_specs_dir" in config file or on command line.}
        );
        exit 1;
    }

    ## Run Expansion Hunter
    say {$FILEHANDLE} q{## Run ExpansionHunter};
    my $sample_sex = $sample_info_href->{sample}{$sample_id}{sex};
    expansionhunter(
        {
            FILEHANDLE            => $FILEHANDLE,
            infile_path           => $file_path_prefix . $infile_suffix,
            json_outfile_path     => $outfile_path_prefix . q{.json},
            log_outfile_path      => $outfile_path_prefix . q{.log},
            reference_genome_path => $reference_genome_path,
            repeat_specs_dir_path => $repeat_specs_dir_path,
            sex                   => $sample_sex,
            vcf_outfile_path      => $outfile_path_prefix . $outfile_suffix,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            infile_path  => $outfile_path_prefix . $ASTERISK,
            outfile_path => $outsample_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait} . $NEWLINE;

    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {

        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => q{expansionhunter},
                path             => catfile(
                    $outsample_directory, $outfile_prefix . $outfile_suffix
                ),
            }
        );

        slurm_submit_job_sample_id_dependency_add_to_sample(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                family_id               => $family_id,
                sample_id               => $sample_id,
                path                    => $job_id_chain,
                log                     => $log,
                sbatch_file_name        => $file_path
            }
        );
    }
    return;
}

sub _get_default_repeat_specs_dir_path {

## Function : Return the path to the repeat specs directory in the Expansionhunter directory
## Returns  : $repeat_specs_dir_path
## Arguments: $log                   => Log
##          : $reference_genome_path => Path to the reference genome used

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $reference_genome_path;

    my $tmpl = {
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        reference_genome_path => {
            defined     => 1,
            required    => 1,
            store       => \$reference_genome_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Cwd qw{ abs_path };
    use File::Basename qw{ fileparse };
    use File::Find::Rule;
    use IPC::Cmd qw{ can_run };

    Readonly my $MINUS_ONE => -1;
    Readonly my $MINUS_TWO => -2;

    ## Path to set
    my $repeat_specs_dir_path;

    ## Get path to binary
    my $expansionhunter_bin_path = can_run(q{ExpansionHunter});

    ## Follow potential link
    $expansionhunter_bin_path = abs_path($expansionhunter_bin_path);

    ## Get the path to the repeat specs dirs
    my @expansionhunter_dirs = File::Spec->splitdir($expansionhunter_bin_path);
    splice @expansionhunter_dirs, $MINUS_TWO;
    my $parent_repeat_specs_dir_path =
      catdir( @expansionhunter_dirs, qw{ data repeat-specs } );

    ## Get list of genome version directories
    my @repeat_specs_dir_paths =
      File::Find::Rule->directory->in($parent_repeat_specs_dir_path);

    ## Remove top directory
    @repeat_specs_dir_paths =
      grep { !/^$parent_repeat_specs_dir_path$/xms } @repeat_specs_dir_paths;

    # Test that some directories has been found
    if ( scalar @repeat_specs_dir_paths == 0 ) {
        $log->fatal(
q{Can't find repeat specification directory for Expansion Hunter. Please set "expansionhunter_repeat_specs_dir" in config file or on command line.}
        );
        exit 1;
    }

    ## Find correct repeat spec folder
    my $genome_reference = fileparse($reference_genome_path);
  REPEAT_SPECS_VERSION:
    foreach my $repeat_specs_version (@repeat_specs_dir_paths) {

        ## Get version
        my @genome_version_dirs = File::Spec->splitdir($repeat_specs_version);
        my $genome_version_dir = splice @genome_version_dirs, $MINUS_ONE;

        ## Match version to reference used
        if ( $genome_reference =~ / $genome_version_dir /ixms ) {
            $repeat_specs_dir_path = $repeat_specs_version;
            last;
        }
    }
    return $repeat_specs_dir_path;
}
1;
