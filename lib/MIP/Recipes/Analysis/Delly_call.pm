package MIP::Recipes::Analysis::Delly_call;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use List::MoreUtils qw { any };
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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_delly_call };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_delly_call {

## Function : Call structural variants using delly
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $insample_directory      => In sample directory
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outsample_directory     => Out sample directory
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $reference_dir           => MIP reference directory
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $insample_directory;
    my $job_id_href;
    my $outsample_directory;
    my $parameter_href;
    my $program_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
    my $reference_dir;
    my $temp_directory;
    my $xargs_file_counter;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        insample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$insample_directory,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
        },
        outsample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outsample_directory,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        reference_dir => {
            default     => $arg_href->{active_parameter_href}{reference_dir},
            strict_type => 1,
            store       => \$reference_dir,
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory,
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ### SET IMPORT IN ALPHABETIC ORDER
    use MIP::Delete::List qw{ delete_contig_elements };
    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Program::Variantcalling::Delly qw{ delly_call };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $program_outdirectory_name =
      $parameter_href->{$mip_program_name}{outdir_name};
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};

    my $program_directory =
      catfile( $outaligner_dir, $program_outdirectory_name );

    ## Filehandles
    # Create anonymous filehandles
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $sample_id,
            program_name          => $program_name,
            program_directory     => $program_directory,
            core_number           => $core_number,
            process_time          => $time,
            temp_directory        => $temp_directory,
        }
    );

    # Used downstream
    $parameter_href->{$mip_program_name}{$sample_id}{indirectory} =
      $outsample_directory;

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Files
    my $infile_tag =
      $file_info_href->{$sample_id}{pgatk_baserecalibration}{file_tag};
    my $outfile_tag =
      $file_info_href->{$sample_id}{$mip_program_name}{file_tag};
    my $infile_prefix  = $merged_infile_prefix . $infile_tag;
    my $outfile_prefix = $merged_infile_prefix . $outfile_tag;

    ## Paths
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ## Assign suffix
    # Get infile_suffix from baserecalibration jobid chain
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
            jobid_chain    => $parameter_href->{pgatk_baserecalibration}{chain},
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

    ### Update contigs
    ## Removes an element from array and return new array while leaving orginal elements_ref untouched
    # Skip contig Y along with MT throughout since sometimes there are no variants particularly for INS
    my @contigs = delete_contig_elements(
        {
            elements_ref       => \@{ $file_info_href->{contigs_size_ordered} },
            remove_contigs_ref => [qw{ MT M Y }],
        }
    );

    # If element is part of array
    if ( ( any { $_ eq q{TRA} } @{ $active_parameter_href->{delly_types} } ) ) {

        ## Required for processing complete file (TRA)
        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        migrate_file(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => catfile(
                    $insample_directory,
                    $infile_prefix . substr( $infile_suffix, 0, 2 ) . $ASTERISK
                ),
                outfile_path => $active_parameter_href->{temp_directory}
            }
        );
    }

    if (
        (
            any { /DEL|DUP|INS|INV/ }
            @{ $active_parameter_href->{delly_types} }
        )
      )
    {
        # If element is part of array
        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        ($xargs_file_counter) = xargs_migrate_contig_files(
            {
                FILEHANDLE        => $FILEHANDLE,
                XARGSFILEHANDLE   => $XARGSFILEHANDLE,
                contigs_ref       => \@contigs,
                file_path         => $file_path,
                program_info_path => $program_info_path,
                core_number       => ( $core_number - 1 )
                ,    #Compensate for cp of entire BAM (INS, TRA), see above
                xargs_file_counter => $xargs_file_counter,
                infile             => $infile_prefix,
                indirectory        => $insample_directory,
                file_ending    => substr( $infile_suffix, 0, 2 ) . $ASTERISK,
                temp_directory => $temp_directory,
            }
        );
        say {$FILEHANDLE} q{wait}, $NEWLINE;
    }

    ## delly
    say {$FILEHANDLE} q{## delly};

    my $xargs_file_path_prefix;

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

    foreach my $sv_type ( @{ $active_parameter_href->{delly_types} } ) {

        if ( $sv_type ne q{TRA} ) {

            ## Process per contig
            foreach my $contig (@contigs) {

                delly_call(
                    {
                        infile_path => $file_path_prefix
                          . $UNDERSCORE
                          . $contig
                          . $infile_suffix,
                        outfile_path => $outfile_path_prefix
                          . $UNDERSCORE
                          . $contig
                          . $UNDERSCORE
                          . $sv_type
                          . $outfile_suffix,
                        stdoutfile_path => $xargs_file_path_prefix
                          . $DOT
                          . $contig
                          . $DOT
                          . $sv_type
                          . $DOT
                          . q{stdout.txt},
                        stderrfile_path => $xargs_file_path_prefix
                          . $DOT
                          . $contig
                          . $DOT
                          . $sv_type
                          . $DOT
                          . q{stderr.txt},
                        sv_type => $sv_type,
                        exclude_file_path =>
                          $active_parameter_href->{delly_exclude_file},
                        referencefile_path =>
                          $active_parameter_href->{human_genome_reference},
                        FILEHANDLE => $XARGSFILEHANDLE,
                    }
                );
                say {$XARGSFILEHANDLE} $NEWLINE;
            }
        }
        else {

            delly_call(
                {
                    infile_path  => $file_path_prefix . $infile_suffix,
                    outfile_path => $outfile_path_prefix
                      . $UNDERSCORE
                      . $sv_type
                      . $outfile_suffix,
                    stdoutfile_path => $xargs_file_path_prefix
                      . $DOT
                      . $sv_type
                      . $DOT
                      . q{stdout.txt},
                    stderrfile_path => $xargs_file_path_prefix
                      . $DOT
                      . $sv_type
                      . $DOT
                      . q{stderr.txt},
                    sv_type => $sv_type,
                    exclude_file_path =>
                      $active_parameter_href->{delly_exclude_file},
                    referencefile_path =>
                      $active_parameter_href->{human_genome_reference},
                    FILEHANDLE => $XARGSFILEHANDLE,
                }
            );
            say {$XARGSFILEHANDLE} $NEWLINE;
        }
    }

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            infile_path => $outfile_path_prefix
              . $ASTERISK
              . $outfile_suffix
              . $ASTERISK,
            outfile_path => $outsample_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});
    close $XARGSFILEHANDLE
      or $log->logcroak(q{Could not close $XARGSFILEHANDLE});

    if ( $mip_program_mode == 1 ) {

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

1;
