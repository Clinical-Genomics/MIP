package MIP::Recipes::Tiddit;

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
use File::Spec::Functions qw{ catdir catfile devnull };

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_tiddit };

}

## Constants
Readonly my $UNDERSCORE => q{_};
Readonly my $SPACE      => q{ };
Readonly my $ASTERISK   => q{*};
Readonly my $NEWLINE    => qq{\n};

sub analysis_tiddit {

## analysis_tiddit

## Function : Call structural variants using tiddit
## Returns  : ""
## Arguments: $parameter_href, $active_parameter_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $sample_id, family_id, $temp_directory, $reference_dir, $outaligner_dir, $call_type, $outfamily_directory
##          : $parameter_href             => Parameter hash {REF}
##          : $active_parameter_href      => Active parameters for this analysis hash {REF}
##          : $sample_info_href           => Info on samples and family hash {REF}
##          : $file_info_href             => The file_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href                => Job id hash {REF}
##          : $family_id                  => Family id
##          : $temp_directory             => Temporary directory
##          : $reference_dir              => MIP reference directory
##          : $outaligner_dir             => Outaligner_dir used in the analysis
##          : $call_type                  => The variant call type
##          : $outfamily_directory     => Out family directory

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id;
    my $temp_directory;
    my $reference_dir;
    my $outaligner_dir;
    my $call_type;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;
    my $outfamily_directory;

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
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory
        },
        reference_dir => {
            default     => $arg_href->{active_parameter_href}{reference_dir},
            strict_type => 1,
            store       => \$reference_dir
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir
        },
        call_type => {
            default     => $UNDERSCORE . q{SV},
            strict_type => 1,
            store       => \$call_type
        },
        outfamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfamily_directory
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Processes qw{ print_wait };
    use MIP::Cluster qw{ get_core_number };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Set::File qw{ set_file_suffix };
    use MIP::Program::Variantcalling::Tiddit qw{ tiddit_sv };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};
    my $modifier_core_number =
      scalar( @{ $active_parameter_href->{sample_ids} } );
    my $max_cores_per_node = $active_parameter_href->{max_cores_per_node};
    my $program_outdirectory_name =
      $parameter_href->{$mip_program_name}{outdir_name};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    $core_number = get_core_number(
        {
            module_core_number   => $core_number,
            modifier_core_number => $modifier_core_number,
            max_cores_per_node   => $max_cores_per_node,
        }
    );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $family_id,
            program_name          => $program_name,
            program_directory =>
              catfile( $outaligner_dir, $program_outdirectory_name ),
            core_number    => $core_number,
            process_time   => $time,
            temp_directory => $temp_directory,
        }
    );

    # Used downstream
    $parameter_href->{$mip_program_name}{indirectory} = $outfamily_directory;

    ## Assign file_tags
    my %file_path_prefix;
    my $outfile_tag =
      $file_info_href->{$family_id}{ q{p} . $program_name }{file_tag};
    my $outfile_prefix = $family_id . $outfile_tag . $call_type;
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ## Assign suffix
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

    my $process_batches_count = 1;

    ## Collect infiles for all sample_ids to enable migration to temporary directory
    while ( my ( $sample_id_index, $sample_id ) =
        each @{ $active_parameter_href->{sample_ids} } )
    {

        ## Assign directories
        my $insample_directory = catdir( $active_parameter_href->{outdata_dir},
            $sample_id, $outaligner_dir );

        ## Add merged infile name after merging all BAM files per sample_id
        my $infile = $file_info_href->{$sample_id}{merge_infile};

        ## Assign file_tags
        my $infile_tag =
          $file_info_href->{$sample_id}{pgatk_baserecalibration}{file_tag};
        my $infile_prefix         = $infile . $infile_tag;
        my $sample_outfile_prefix = $infile . $outfile_tag;

        #q{.bam} -> ".b*" for getting index as well
        my $infile_path = catfile( $insample_directory,
            $infile_prefix . substr( $infile_suffix, 0, 2 ) . $ASTERISK );

        $file_path_prefix{$sample_id}{in} =
          catfile( $temp_directory, $infile_prefix );
        $file_path_prefix{$sample_id}{out} =
          catfile( $temp_directory, $sample_outfile_prefix );

        $process_batches_count = print_wait(
            {
                process_counter       => $sample_id_index,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                FILEHANDLE            => $FILEHANDLE,
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
                process_counter       => $sample_id_index,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                FILEHANDLE            => $FILEHANDLE,
            }
        );

        ## Tiddit
        tiddit_sv(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $file_path_prefix{$sample_id}{in}
                  . $infile_suffix,
                outfile_path_prefix => $file_path_prefix{$sample_id}{out},
                minimum_number_supporting_pairs => $active_parameter_href
                  ->{tiddit_minimum_number_supporting_pairs},
            }
        );
        say {$FILEHANDLE} q{&} . $SPACE . $NEWLINE;
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Get parameters
    ## Tiddit sample outfiles
    my @infile_paths = map { $file_path_prefix{$_}{out} . $outfile_suffix }
      ( keys %file_path_prefix );

    Program::Variantcalling::Svdb::merge(
        {
            infile_paths_ref => \@infile_paths,
            outfile_path     => $outfile_path_prefix . $outfile_suffix,
            FILEHANDLE       => $FILEHANDLE,
            notag            => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            infile_path  => $outfile_path_prefix . $outfile_suffix . $ASTERISK,
            outfile_path => $outfamily_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {

        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => q{tiddit},
                outdirectory     => $outfamily_directory,
                outfile          => $outfile_prefix . $outfile_suffix,
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

1;
