package MIP::Recipes::Analysis::Vt;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile devnull splitpath };
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
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_vt analysis_vt_rio };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $DASH       => q{-};
Readonly my $DOT        => q{.};
Readonly my $EMPTY_STR  => q{};
Readonly my $NEWLINE    => qq{\n};
Readonly my $PIPE       => q{|};
Readonly my $SEMICOLON  => q{;};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_vt {

## Function : Split multi allelic records into single records and normalize
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => Variant call type
##          : $family_id               => Family id
##          : $file_path               => File path
##          : $file_info_href          => File_info hash {REF}
##          : $infamily_directory      => In family directory
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outfamily_directory     => Out family directory
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $program_info_path       => The program info path
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $stderr_path             => The stderr path of the block script
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_path;
    my $file_info_href;
    my $infamily_directory;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $outfamily_directory;
    my $parameter_href;
    my $program_name;
    my $program_info_path;
    my $sample_info_href;
    my $stderr_path;

    ## Default(s)
    my $call_type;
    my $family_id;
    my $outaligner_dir;
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
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
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
        file_path          => { strict_type => 1, store => \$file_path, },
        infamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infamily_directory,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
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
        outfamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfamily_directory,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        program_name      => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        stderr_path    => { strict_type => 1, store => \$stderr_path },
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

    use MIP::Cluster qw{ get_core_number };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Vt_core qw{ analysis_vt_core_rio};
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $program_mode = $active_parameter_href->{$program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$program_name}{chain};

    my ( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name      => $program_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            module_core_number   => $core_number,
            modifier_core_number => scalar @{ $file_info_href->{contigs} },
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            call_type                       => $call_type,
            core_number                     => $core_number,
            directory_id                    => $family_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            program_directory               => $outaligner_dir,
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );
    $stderr_path = $program_info_path . $DOT . q{stderr.txt};

    # Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) = splitpath($stderr_path);

    # Used downstream
    $parameter_href->{$program_name}{indirectory} = $outfamily_directory;

    ## Tags
    my $infile_tag = $file_info_href->{$family_id}{rhocall}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$program_name}{file_tag};

    ## Files
    my $infile_prefix  = $family_id . $infile_tag . $call_type;
    my $outfile_prefix = $family_id . $outfile_tag . $call_type;

    ## Paths
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain

    my $infile_suffix = get_file_suffix(
        {
            jobid_chain    => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );
    my $outfile_suffix = set_file_suffix(
        {
            file_suffix => $parameter_href->{$program_name}{outfile_suffix},
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    ($xargs_file_counter) = xargs_migrate_contig_files(
        {
            contigs_ref        => \@{ $file_info_href->{contigs_size_ordered} },
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $file_path,
            indirectory        => $infamily_directory,
            infile             => $infile_prefix,
            program_info_path  => $program_info_path,
            temp_directory     => $temp_directory,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    say {$FILEHANDLE}
q{## vt - Decompose (split multi allelic records into single records) and/or normalize variants};

    my $xargs_file_path_prefix;

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Split vcf into contigs
    while ( my ( $contig_index, $contig ) =
        each @{ $file_info_href->{contigs_size_ordered} } )
    {

        ## vt - Split multi allelic records into single records and normalize
        analysis_vt_core_rio(
            {
                active_parameter_href => $active_parameter_href,
                cmd_break             => $SEMICOLON,
                contig                => $contig,
                decompose             => $active_parameter_href->{vt_decompose},
                FILEHANDLE            => $XARGSFILEHANDLE,
                gnu_sed               => 1,
                infile_path           => $file_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $infile_suffix,
                instream     => 0,
                normalize    => $active_parameter_href->{vt_normalize},
                outfile_path => $outfile_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $outfile_suffix,
                uniq                   => $active_parameter_href->{vt_uniq},
                xargs_file_path_prefix => $xargs_file_path_prefix,
            }
        );

        if (   $contig_index == 0
            && $program_mode == 1 )
        {

            ## Split to enable submission to &SampleInfoQC later
            my ( $volume_xargs, $directory_xargs, $stderr_file_xargs ) =
              splitpath($xargs_file_path_prefix);

            ## Collect QC metadata info for later use
            my $qc_vt_outfile =
              $stderr_file_xargs . $DOT . $contig . $DOT . q{stderr.txt};
            add_program_outfile_to_sample_info(
                {
                    path             => catfile( $directory, $qc_vt_outfile ),
                    program_name     => q{vt},
                    sample_info_href => $sample_info_href,
                }
            );
        }

        my $alt_file_tag = $EMPTY_STR;

        ## VEP does not annotate '*' since the alt allele does not exist, this is captured in the upstream indel and SNV record associated with '*'
        ## Remove decomposed '*' entries
        if ( $active_parameter_href->{vt_missing_alt_allele} ) {

            # Update file tag
            $alt_file_tag .= $UNDERSCORE . q{nostar};

            _remove_decomposed_asterisk_entries(
                {
                    alt_file_tag           => $UNDERSCORE . q{nostar},
                    contig                 => $contig,
                    outfile_prefix         => $outfile_prefix,
                    outfile_suffix         => $outfile_suffix,
                    outfile_path_prefix    => $outfile_path_prefix,
                    temp_directory         => $temp_directory,
                    XARGSFILEHANDLE        => $XARGSFILEHANDLE,
                    xargs_file_path_prefix => $xargs_file_path_prefix,
                }
            );
        }

        gnu_mv(
            {
                FILEHANDLE  => $XARGSFILEHANDLE,
                infile_path => $outfile_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $alt_file_tag
                  . $outfile_suffix,
                outfile_path => $outfile_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $outfile_suffix,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $outfile_path_prefix
              . $UNDERSCORE
              . $ASTERISK
              . $outfile_suffix
              . $ASTERISK,
            outfile_path => $outfamily_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;
    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $program_mode == 1 ) {

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

sub analysis_vt_rio {

## Function : Split multi allelic records into single records and normalize
## Returns  : |xargs_file_counter
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => Variant call type
##          : $family_id               => Family id
##          : $FILEHANDLE              => Filehandle to write to
##          : $file_path               => File path
##          : $file_info_href          => File_info hash {REF}
##          : $infamily_directory      => In family directory
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outfamily_directory     => Out family directory
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $program_info_path       => The program info path
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $stderr_path             => The stderr path of the block script
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $FILEHANDLE;
    my $file_path;
    my $file_info_href;
    my $infamily_directory;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $outfamily_directory;
    my $parameter_href;
    my $program_name;
    my $program_info_path;
    my $sample_info_href;
    my $stderr_path;

    ## Default(s)
    my $call_type;
    my $family_id;
    my $outaligner_dir;
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
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        FILEHANDLE     => { store => \$FILEHANDLE, },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        file_path          => { strict_type => 1, store => \$file_path, },
        infamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infamily_directory,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
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
        outfamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfamily_directory,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        program_name      => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        stderr_path    => { strict_type => 1, store => \$stderr_path },
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

    use MIP::Cluster qw{ get_core_number };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Vt_core qw{ analysis_vt_core_rio};
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Set::File qw{ set_file_suffix };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $program_mode = $active_parameter_href->{$program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$program_name}{chain};

    my ( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name      => $program_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            module_core_number   => $core_number,
            modifier_core_number => scalar @{ $file_info_href->{contigs} },
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    # Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) = splitpath($stderr_path);

    # Used downstream
    $parameter_href->{$program_name}{indirectory} = $outfamily_directory;

    ## Tags
    my $infile_tag = $file_info_href->{$family_id}{rhocall}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$program_name}{file_tag};

    ## Files
    my $infile_prefix  = $family_id . $infile_tag . $call_type;
    my $outfile_prefix = $family_id . $outfile_tag . $call_type;

    ## Paths
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain

    my $infile_suffix = get_file_suffix(
        {
            jobid_chain    => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );
    my $outfile_suffix = set_file_suffix(
        {
            file_suffix => $parameter_href->{$program_name}{outfile_suffix},
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

    say {$FILEHANDLE}
q{## vt - Decompose (split multi allelic records into single records) and/or normalize variants};

    my $xargs_file_path_prefix;

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Split vcf into contigs
    while ( my ( $contig_index, $contig ) =
        each @{ $file_info_href->{contigs_size_ordered} } )
    {

        ## vt - Split multi allelic records into single records and normalize
        analysis_vt_core_rio(
            {
                active_parameter_href => $active_parameter_href,
                cmd_break             => $SEMICOLON,
                contig                => $contig,
                decompose             => $active_parameter_href->{vt_decompose},
                FILEHANDLE            => $XARGSFILEHANDLE,
                gnu_sed               => 1,
                infile_path           => $file_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $infile_suffix,
                instream     => 0,
                normalize    => $active_parameter_href->{vt_normalize},
                outfile_path => $outfile_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $outfile_suffix,
                uniq                   => $active_parameter_href->{vt_uniq},
                xargs_file_path_prefix => $xargs_file_path_prefix,
            }
        );

        if (   $contig_index == 0
            && $program_mode == 1 )
        {

            #Split to enable submission to &SampleInfoQC later
            my ( $volume_xargs, $directory_xargs, $stderr_file_xargs ) =
              splitpath($xargs_file_path_prefix);

            ## Collect QC metadata info for later use
            my $qc_vt_outfile =
              $stderr_file_xargs . $DOT . $contig . $DOT . q{stderr.txt};
            add_program_outfile_to_sample_info(
                {
                    path             => catfile( $directory, $qc_vt_outfile ),
                    program_name     => q{vt},
                    sample_info_href => $sample_info_href,
                }
            );
        }

        my $alt_file_tag = $EMPTY_STR;

        ## VEP does not annotate '*' since the alt allele does not exist, this is captured in the upstream indel and SNV record associated with '*'
        ## Remove decomposed '*' entries
        if ( $active_parameter_href->{vt_missing_alt_allele} ) {

            # Update file tag
            $alt_file_tag .= $UNDERSCORE . q{nostar};

            _remove_decomposed_asterisk_entries(
                {
                    alt_file_tag           => $UNDERSCORE . q{nostar},
                    contig                 => $contig,
                    outfile_prefix         => $outfile_prefix,
                    outfile_suffix         => $outfile_suffix,
                    outfile_path_prefix    => $outfile_path_prefix,
                    temp_directory         => $temp_directory,
                    XARGSFILEHANDLE        => $XARGSFILEHANDLE,
                    xargs_file_path_prefix => $xargs_file_path_prefix,
                }
            );
        }

        gnu_mv(
            {
                FILEHANDLE  => $XARGSFILEHANDLE,
                infile_path => $outfile_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $alt_file_tag
                  . $outfile_suffix,
                outfile_path => $outfile_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $outfile_suffix,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    # Track the number of created xargs scripts per module for Block algorithm
    return $xargs_file_counter;
}

sub _remove_decomposed_asterisk_entries {

## Function : Remove decomposed '*' entries
## Returns  :
## Arguments: $alt_file_tag           => Alt. file tag
##          : $contig                 => contig
##          : $outfile_prefix         => Outfile prefix
##          : $outfile_suffix         => Outfile suffix
##          : $outfile_path_prefix    => Outfile path prefix
##          : $temp_directory         => Temporary directory
##          : $XARGSFILEHANDLE        => XARGS file handle
##          : $xargs_file_path_prefix => Xargs file path prefix

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $alt_file_tag;
    my $contig;
    my $outfile_path_prefix;
    my $outfile_prefix;
    my $outfile_suffix;
    my $temp_directory;
    my $XARGSFILEHANDLE;
    my $xargs_file_path_prefix;

    my $tmpl = {
        alt_file_tag => {
            required => 1,
            defined  => 1,
            store    => \$alt_file_tag,
        },
        contig => {
            required => 1,
            defined  => 1,
            store    => \$contig,
        },
        outfile_path_prefix => {
            required => 1,
            defined  => 1,
            store    => \$outfile_path_prefix,
        },
        outfile_prefix => {
            required => 1,
            defined  => 1,
            store    => \$outfile_prefix,
        },
        outfile_suffix => {
            required => 1,
            defined  => 1,
            store    => \$outfile_suffix,
        },
        temp_directory => {
            required => 1,
            defined  => 1,
            store    => \$temp_directory,
        },
        XARGSFILEHANDLE => {
            required => 1,
            defined  => 1,
            store    => \$XARGSFILEHANDLE,
        },
        xargs_file_path_prefix => {
            required => 1,
            defined  => 1,
            store    => \$xargs_file_path_prefix,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Execute perl
    my $remove_star_regexp = q?perl -nae \'?;

    ## Print if line does not contain asterisk
    $remove_star_regexp .= q?unless\($F\[4\] eq \"\*\") \{print $_\}\' ?;

    $alt_file_tag = $UNDERSCORE . q{nostar};

    print {$XARGSFILEHANDLE} catfile( $remove_star_regexp . $temp_directory,
        $outfile_prefix . $UNDERSCORE . $contig . $outfile_suffix )
      . $SPACE;

    print {$XARGSFILEHANDLE} q{>}
      . $SPACE
      . $outfile_path_prefix
      . $UNDERSCORE
      . $contig
      . $alt_file_tag
      . $outfile_suffix
      . $SPACE;

    print {$XARGSFILEHANDLE} q{2>>}
      . $SPACE
      . $xargs_file_path_prefix
      . $DOT
      . $contig
      . $DOT
      . q{stderr.txt}
      . $SPACE;

    # Redirect xargs output to program specific stderr file
    print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

    return;
}

1;
