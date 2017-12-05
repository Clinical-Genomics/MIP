package MIP::Recipes::Analysis::Vt;

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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_vt };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $DASH       => q{-};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $PIPE       => q{|};
Readonly my $SEMICOLON  => q{;};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_vt {

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
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Genmod
      qw{ genmod_annotate genmod_filter };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Unpack parameters
    my $job_id_chain  = $parameter_href->{$mip_program_name}{chain};
    my $reduce_io_ref = \$active_parameter_href->{reduce_io};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};

    ## Filehandles
    if ( not defined $FILEHANDLE ) {

        #Create anonymous filehandle
        $FILEHANDLE = IO::Handle->new();
    }

    # Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    my $core_number = get_core_number(
        {
            module_core_number =>
              $active_parameter_href->{module_core_number}{$mip_program_name},
            modifier_core_number => scalar( @{ $file_info_href->{contigs} } ),
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    #Run as individual sbatch script
    if ( !${$reduce_io_ref} ) {

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        ( $file_path, $program_info_path ) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                job_id_href           => $job_id_href,
                FILEHANDLE            => $FILEHANDLE,
                directory_id          => $family_id,
                program_name          => $program_name,
                program_directory     => $outaligner_dir,
                call_type             => $call_type,
                core_number           => $core_number,
                process_time          => $time,
                temp_directory        => $temp_directory,
            }
        );
        $stderr_path = $program_info_path . $DOT . q{stderr.txt};
    }

    #Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) = splitpath($stderr_path);

    #Used downstream
    $parameter_href->{ $mip_program_name }{indirectory} =
      $outfamily_directory;

    ## Tags
    my $infile_tag = $file_info_href->{$family_id}{prhocall}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$mip_program_name}{file_tag};

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
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
            jobid_chain    => $job_id_chain,
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

    if ( !${$reduce_io_ref} ) {    #Run as individual sbatch script

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        ($xargs_file_counter) = xargs_migrate_contig_files(
            {
                FILEHANDLE      => $FILEHANDLE,
                XARGSFILEHANDLE => $XARGSFILEHANDLE,
                contigs_ref => \@{ $file_info_href->{contigs_size_ordered} },
                file_path   => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
                xargs_file_counter => $xargs_file_counter,
                infile             => $infile_prefix,
                indirectory        => $infamily_directory,
                temp_directory     => $temp_directory,
            }
        );
    }

    say {$FILEHANDLE}
q{## vt - Decompose (split multi allelic records into single records) and/or normalize variants};

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

#VEP does not annotate '*' since the alt allele does not exist, this is captured in the upstream indel and SNV record associated with '*'
    my $remove_star_regexp =
      q?perl -nae \'unless\($F\[4\] eq \"\*\") \{print $_\}\' ?;

    ## Split vcf into contigs
    while ( my ( $contig_index, $contig ) =
        each @{ $file_info_href->{contigs_size_ordered} } )
    {

        ## vt - Split multi allelic records into single records and normalize
        analysis_vt_core_rio(
            {
                active_parameter_href => $active_parameter_href,
                FILEHANDLE            => $XARGSFILEHANDLE,
                infile_path           => $file_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $infile_suffix,
                outfile_path => $outfile_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $outfile_suffix,
                decompose => $active_parameter_href->{vt_decompose},
                normalize => $active_parameter_href->{vt_normalize},
                uniq      => $active_parameter_href->{vt_uniq},
                gnu_sed   => 1,
                instream  => 0,
                cmd_break => $SEMICOLON,
                xargs_file_path_prefix => $xargs_file_path_prefix,
                contig                 => $contig,
            }
        );

        if (   ( $contig_index == 0 )
            && ( $active_parameter_href->{$mip_program_name} == 1 ) )
        {

            my ( $volume, $directory, $stderr_file ) =
              splitpath($xargs_file_path_prefix)
              ;    #Split to enable submission to &SampleInfoQC later

            ## Collect QC metadata info for later use
            my $qc_vt_outfile =
              $stderr_file . $DOT . $contig . $DOT . q{stderr.txt};
            add_program_outfile_to_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    program_name     => q{vt},
                    outdirectory     => $directory,
                    outfile          => $qc_vt_outfile,
                }
            );
        }

        my $alt_file_tag = q{};

        ## Remove decomposed '*' entries
        if ( $active_parameter_href->{vt_missing_alt_allele} ) {

            $alt_file_tag = $UNDERSCORE . q{nostar};
            print {$XARGSFILEHANDLE}
              catfile( $remove_star_regexp . $temp_directory,
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
              . $SPACE;   #Redirect xargs output to program specific stderr file
            print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;
        }

        ## Remove common variants
        if ( $active_parameter_href->{vt_genmod_filter} ) {

            genmod_annotate(
                {
                    infile_path => $outfile_path_prefix
                      . $UNDERSCORE
                      . $contig
                      . $alt_file_tag
                      . $outfile_suffix,
                    outfile_path => catfile( dirname( devnull() ), q{stdout} ),
                    stderrfile_path => $xargs_file_path_prefix
                      . $DOT
                      . $contig
                      . $DOT
                      . q{stderr.txt},
                    verbosity           => q{v},
                    temp_directory_path => $temp_directory,
                    thousand_g_file_path =>
                      $active_parameter_href->{vt_genmod_filter_1000g},
                    max_af => $active_parameter_href->{vt_genmod_filter_max_af},
                    FILEHANDLE => $XARGSFILEHANDLE,
                }
            );
            print {$XARGSFILEHANDLE} $PIPE . $SPACE;

            #Update file tag
            $alt_file_tag .= $UNDERSCORE . q{genmod_filter};

            genmod_filter(
                {
                    infile_path  => $DASH,
                    outfile_path => $outfile_path_prefix
                      . $UNDERSCORE
                      . $contig
                      . $alt_file_tag
                      . $outfile_suffix,
                    stderrfile_path_append => $xargs_file_path_prefix
                      . $DOT
                      . $contig
                      . $DOT
                      . q{stderr.txt},
                    verbosity => q{v},
                    threshold =>
                      $active_parameter_href->{sv_genmod_filter_threshold},
                    FILEHANDLE => $XARGSFILEHANDLE,
                }
            );
            print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;
        }

        gnu_mv(
            {
                infile_path => $outfile_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $alt_file_tag
                  . $outfile_suffix,
                outfile_path => $outfile_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $outfile_suffix,
                FILEHANDLE => $XARGSFILEHANDLE,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
        close $XARGSFILEHANDLE
          or $log->logcroak(q{Could not close $XARGSFILEHANDLE});
    }

    #Run as individual sbatch script
    if ( !${$reduce_io_ref} ) {

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                infile_path => $outfile_path_prefix
                  . $UNDERSCORE
                  . $ASTERISK
                  . $outfile_suffix
                  . $ASTERISK,
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} q{wait}, $NEWLINE;
        close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});
    }

    if ( $mip_program_mode == 1 ) {

        #Run as individual sbatch script
        if ( !${$reduce_io_ref} ) {

            slurm_submit_job_sample_id_dependency_add_to_family(
                {
                    job_id_href             => $job_id_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    sample_ids_ref =>
                      \@{ $active_parameter_href->{sample_ids} },
                    family_id        => $family_id,
                    path             => $job_id_chain,
                    log              => $log,
                    sbatch_file_name => $file_path,
                }
            );
        }
    }
    if ( ${$reduce_io_ref} ) {

       #Track the number of created xargs scripts per module for Block algorithm
        return $xargs_file_counter;
    }
}

1;
