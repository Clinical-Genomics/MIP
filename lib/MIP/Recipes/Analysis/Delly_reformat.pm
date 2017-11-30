package MIP::Recipes::Analysis::Delly_reformat;

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
    our @EXPORT_OK = qw{ analysis_delly_reformat };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SEMICOLON  => q{;};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_delly_reformat {

## Function : Merge, regenotype, and filter using Delly
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => Variant call type
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outfamily_directory     => Out family directory
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $reference_dir           => MIP reference directory
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

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
    my $call_type;
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
        call_type =>
          { default => q{SV}, strict_type => 1, store => \$call_type, },
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
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        reference_dir_ref => {
            default     => $arg_href->{active_parameter_href}{reference_dir},
            strict_type => 1,
            store       => \$reference_dir,
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
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_merge bcftools_index bcftools_concat };
    use MIP::Program::Variantcalling::Delly
      qw{ delly_call delly_merge delly_filter };
    use MIP::Program::Variantcalling::Picardtools qw{ picardtools_sortvcf };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};
    my $program_outdirectory_name =
      $parameter_href->{$mip_program_name}{outdir_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};

    ## Filehandles
    # Create anonymous filehandles
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    my $program_directory =
      catfile( $outaligner_dir, $program_outdirectory_name );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $family_id,
            program_name          => $program_name,
            program_directory     => $program_directory,
            core_number           => $core_number,
            process_time          => $time,
            temp_directory        => $temp_directory,
        }
    );

    # Used downstream
    $parameter_href->{$mip_program_name}{indirectory} = $outfamily_directory;

    ## Tags
    my $outfile_tag =
      $file_info_href->{$family_id}{$mip_program_name}{file_tag};

    ## Files
    my $outfile_prefix = $family_id . $outfile_tag . $UNDERSCORE . $call_type;
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );
    my $file_suffix = $parameter_href->{$mip_program_name}{outfile_suffix};

    ### Assign suffix
    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
            job_id_chain   => $job_id_chain,
            file_suffix    => $file_suffix,
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

    ## Collect files and suffix for all sample_ids
    my %infile_path_prefix;
    my %file_path_prefix;
    my %suffix;
    my @program_tag_keys = (qw{ pgatk_baserecalibration pdelly_call });

    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

      say {$FILEHANDLE} q{@@@@@@@@@@@ Sample ---> } . $sample_id;

        ## Assign directories
        my $insample_directory_bam =
          catdir( $active_parameter_href->{outdata_dir},
            $sample_id, $outaligner_dir );

        my $insample_directory_bcf =
          catdir( $active_parameter_href->{outdata_dir},
            $sample_id, $outaligner_dir,
            $parameter_href->{pdelly_call}{outdir_name} );

        foreach my $infile_tag_key (@program_tag_keys) {

            ## Add merged infile name prefix after merging all BAM files per sample_id
            my $merged_infile_prefix = get_merged_infile_prefix(
                {
                    file_info_href => $file_info_href,
                    sample_id      => $sample_id,
                }
            );

            ## Tags
            my $infile_tag =
              $file_info_href->{$sample_id}{$infile_tag_key}{file_tag};
            my $infile_prefix = $merged_infile_prefix . $infile_tag;

            # Used downstream
            $infile_path_prefix{$sample_id}{$infile_tag_key} =
              catfile( $temp_directory, $infile_prefix );
            $file_path_prefix{$sample_id} =
              catfile( $temp_directory, $merged_infile_prefix . $outfile_tag );

            # BCFs
            if ( $infile_tag_key eq q{pdelly_call} ) {

                ## Assign suffix
                $suffix{$infile_tag_key} = get_file_suffix(
                    {
                        parameter_href => $parameter_href,
                        suffix_key     => q{outfile_suffix},
                        program_name   => $infile_tag_key,
                    }
                );

              SV_TYPE:
                foreach
                  my $sv_type ( @{ $active_parameter_href->{delly_types} } )
                {

                    my $file_ending =
                        $UNDERSCORE
                      . $sv_type
                      . substr( $suffix{$infile_tag_key}, 0, 2 )
                      . $ASTERISK;

                    if ( $sv_type ne q{TRA} ) {

                        ## Copy file(s) to temporary directory
                        say {$FILEHANDLE}
                          q{## Copy file(s) to temporary directory};
                        ($xargs_file_counter) = xargs_migrate_contig_files(
                            {
                                FILEHANDLE         => $FILEHANDLE,
                                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                                contigs_ref        => \@contigs,
                                file_path          => $file_path,
                                program_info_path  => $program_info_path,
                                core_number        => $core_number,
                                xargs_file_counter => $xargs_file_counter,
                                infile => $merged_infile_prefix . $infile_tag,
                                indirectory    => $insample_directory_bcf,
                                file_ending    => $file_ending,
                                temp_directory => $temp_directory,
                            }
                        );
                    }
                    else {

                        my $infile_path = catfile( $insample_directory_bcf,
                                $merged_infile_prefix
                              . $infile_tag
                              . $file_ending );
                        say {$FILEHANDLE}
                          q{## Copy file(s) to temporary directory};
                        migrate_file(
                            {
                                FILEHANDLE  => $FILEHANDLE,
                                infile_path => $infile_path,
                                outfile_path =>
                                  $active_parameter_href->{temp_directory},
                            }
                        );
                    }
                }
            }
            else {    #BAMs

                ## Assign suffix
                $suffix{$infile_tag_key} = get_file_suffix(
                    {
                        parameter_href => $parameter_href,
                        suffix_key     => q{alignment_file_suffix},
                        jobid_chain =>
                          $parameter_href->{$infile_tag_key}{chain},
                    }
                );

                my $file_ending =
                  substr( $suffix{$infile_tag_key}, 0, 2 ) . $ASTERISK;

                ## Copy file(s) to temporary directory
                say {$FILEHANDLE} q{## Copy file(s) to temporary directory};

                my $infile_path = catfile( $insample_directory_bam,
                    $merged_infile_prefix . $infile_tag . $file_ending );

                migrate_file(
                    {
                        FILEHANDLE  => $FILEHANDLE,
                        infile_path => $infile_path,
                        outfile_path =>
                          $active_parameter_href->{temp_directory},
                    }
                );

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
                        ,    # Compensate for cp of entire BAM TRA, see above
                        xargs_file_counter => $xargs_file_counter,
                        infile         => $merged_infile_prefix . $infile_tag,
                        indirectory    => $insample_directory_bam,
                        file_ending    => $file_ending,
                        temp_directory => $temp_directory,
                    }
                );
            }
            say {$FILEHANDLE} q{wait}, $NEWLINE;
        }
    }

    if ( scalar @{ $active_parameter_href->{sample_ids} } > 1 ) {

        say {$FILEHANDLE} q{@@@@@@@@@@@ ---> IM IN IF SCALAR};

        ### Delly merge
        say {$FILEHANDLE} q{## delly merge} . $NEWLINE;

        say {$FILEHANDLE}
          q{## Fix locale bug using old centosOS and Boost library};
        say {$FILEHANDLE} q?LC_ALL="C"; export LC_ALL ?, $NEWLINE . $NEWLINE;

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

      SV_TYPE:
        foreach my $sv_type ( @{ $active_parameter_href->{delly_types} } ) {

            say {$FILEHANDLE} q{@@@@@@@@@@@ LINE 439, SV TYPE ---> } . $sv_type;

            if ( $sv_type ne q{TRA} ) {

              CONTIG:
                foreach my $contig (@contigs) {

                    say {$FILEHANDLE} q{@@@@@@@@@@@ LINE 446, CONTIG ---> } . $contig;

                    ## Assemble file paths by adding file ending
                    my @file_paths = map {
                            $infile_path_prefix{$_}{pdelly_call}
                          . $UNDERSCORE
                          . $contig
                          . $UNDERSCORE
                          . $sv_type
                          . $suffix{pdelly_call}
                    } @{ $active_parameter_href->{sample_ids} };

                    delly_merge(
                        {
                            infile_paths_ref => \@file_paths,
                            outfile_path     => $outfile_path_prefix
                              . $UNDERSCORE
                              . $contig
                              . $UNDERSCORE
                              . $sv_type
                              . $DOT . q{bcf},
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
                            sv_type    => $sv_type,
                            min_size   => 0,
                            max_size   => 100000000,
                            FILEHANDLE => $XARGSFILEHANDLE,
                        }
                    );
                    say {$XARGSFILEHANDLE} $NEWLINE;
                }
            }
            else {

                ## Assemble file paths by adding file ending
                my @file_paths = map {
                        $infile_path_prefix{$_}{pdelly_call}
                      . $UNDERSCORE
                      . $sv_type
                      . $suffix{pdelly_call}
                } @{ $active_parameter_href->{sample_ids} };

                delly_merge(
                    {
                        infile_paths_ref => \@file_paths,
                        outfile_path     => $outfile_path_prefix
                          . $UNDERSCORE
                          . $sv_type
                          . $DOT . q{bcf},
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
                        sv_type    => $sv_type,
                        min_size   => 0,
                        max_size   => 100000000,
                        FILEHANDLE => $XARGSFILEHANDLE,
                    }
                );
                say {$XARGSFILEHANDLE} $NEWLINE;
            }
        }

        ## Delly call regenotype
        say {$FILEHANDLE} q{## delly call regenotype};

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

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
          SV_TYPE:
            foreach my $sv_type ( @{ $active_parameter_href->{delly_types} } ) {

                if ( $sv_type ne q{TRA} ) {

                  CONTIG:
                    foreach my $contig (@contigs) {

                        ## Assemble file path
                        my $alignment_sample_file_path =
                          $infile_path_prefix{$sample_id}
                          {pgatk_baserecalibration}
                          . $UNDERSCORE
                          . $contig
                          . $suffix{pgatk_baserecalibration};

                        delly_call(
                            {
                                infile_path => $alignment_sample_file_path,
                                genotypefile_path => $outfile_path_prefix
                                  . $UNDERSCORE
                                  . $contig
                                  . $UNDERSCORE
                                  . $sv_type
                                  . $suffix{pdelly_call},
                                outfile_path => $file_path_prefix{$sample_id}
                                  . $UNDERSCORE
                                  . $contig
                                  . $UNDERSCORE
                                  . $sv_type
                                  . $UNDERSCORE . q{geno}
                                  . $suffix{pdelly_call},
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
                                referencefile_path => $referencefile_path,
                                FILEHANDLE         => $XARGSFILEHANDLE,
                            }
                        );
                        say {$XARGSFILEHANDLE} $NEWLINE;
                    }
                }
                else {

                    ## Assemble file path
                    my $alignment_sample_file_path =
                        $infile_path_prefix{$sample_id}{pgatk_baserecalibration}
                      . $suffix{pgatk_baserecalibration};

                    delly_call(
                        {
                            infile_path       => $alignment_sample_file_path,
                            genotypefile_path => $outfile_path_prefix
                              . $UNDERSCORE
                              . $sv_type
                              . $suffix{pdelly_call},
                            outfile_path => $file_path_prefix{$sample_id}
                              . $UNDERSCORE
                              . $sv_type
                              . $UNDERSCORE . q{geno}
                              . $suffix{pdelly_call},
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
                            referencefile_path => $referencefile_path,
                            FILEHANDLE         => $XARGSFILEHANDLE,
                        }
                    );
                    say {$XARGSFILEHANDLE} $NEWLINE;
                }
            }

            ### Merge calls
            say {$FILEHANDLE} q{## bcftools merge};

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

          SV_TYPE:
            foreach my $sv_type ( @{ $active_parameter_href->{delly_types} } ) {

                if ( $sv_type ne q{TRA} ) {

                  CONTIG:
                    foreach my $contig (@contigs) {

                        ## Assemble file paths by adding file ending
                        my @file_paths = map {
                                $file_path_prefix{$_}
                              . $UNDERSCORE
                              . $contig
                              . $UNDERSCORE
                              . $sv_type
                              . $UNDERSCORE . q{geno}
                              . $suffix{pdelly_call}
                        } @{ $active_parameter_href->{sample_ids} };

                        bcftools_merge(
                            {
                                infile_paths_ref => \@file_paths,
                                outfile_path     => $outfile_path_prefix
                                  . $UNDERSCORE
                                  . $contig
                                  . $UNDERSCORE
                                  . $sv_type
                                  . $suffix{pdelly_call},
                                output_type     => q{b},
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
                                FILEHANDLE => $XARGSFILEHANDLE,
                            }
                        );
                        print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

                        bcftools_index(
                            {
                                infile_path => $outfile_path_prefix
                                  . $UNDERSCORE
                                  . $contig
                                  . $UNDERSCORE
                                  . $sv_type
                                  . $suffix{pdelly_call},
                                stderrfile_path => $xargs_file_path_prefix
                                  . $DOT
                                  . $contig
                                  . $DOT
                                  . $sv_type
                                  . $UNDERSCORE
                                  . q{index.stderr.txt},
                                output_type => q{csi},
                                FILEHANDLE  => $XARGSFILEHANDLE,
                            }
                        );
                        say {$XARGSFILEHANDLE} $NEWLINE;
                    }
                }
                else {

                    ## Assemble file paths by adding file ending
                    my @file_paths = map {
                            $file_path_prefix{$_}
                          . $UNDERSCORE
                          . $sv_type
                          . $UNDERSCORE . q{geno}
                          . $suffix{pdelly_call}
                    } @{ $active_parameter_href->{sample_ids} };

                    bcftools_merge(
                        {
                            infile_paths_ref => \@file_paths,
                            outfile_path     => $outfile_path_prefix
                              . $UNDERSCORE
                              . $sv_type
                              . $suffix{pdelly_call},
                            output_type     => q{b},
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
                            FILEHANDLE => $XARGSFILEHANDLE,
                        }
                    );
                    print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

                    bcftools_index(
                        {
                            infile_path => $outfile_path_prefix
                              . $UNDERSCORE
                              . $sv_type
                              . $suffix{pdelly_call},
                            stderrfile_path => $xargs_file_path_prefix
                              . $DOT
                              . $sv_type
                              . $UNDERSCORE
                              . q{index.stderr.txt},
                            output_type => q{csi},
                            FILEHANDLE  => $XARGSFILEHANDLE,
                        }
                    );
                    say {$XARGSFILEHANDLE} $NEWLINE;
                }
            }

            ### Concatenate SV types
            say {$FILEHANDLE}
              q{## bcftools concat - concatenate SV type per contigs};

            ## Assemble file paths by adding file ending
            my @file_paths;

          SV_TYPE:
            foreach my $sv_type ( @{ $active_parameter_href->{delly_types} } ) {

                if ( $sv_type ne q{TRA} ) {

                    push @file_paths, map {
                            $outfile_path_prefix
                          . $UNDERSCORE
                          . $_
                          . $UNDERSCORE
                          . $sv_type
                          . $suffix{pdelly_call}
                    } @contigs;
                }
                else {

                    push @file_paths,
                        $outfile_path_prefix
                      . $UNDERSCORE
                      . $sv_type
                      . $suffix{pdelly_call};
                }

                bcftools_concat(
                    {
                        infile_paths_ref => \@file_paths,
                        outfile_path     => $outfile_path_prefix
                          . $UNDERSCORE
                          . $sv_type
                          . $UNDERSCORE
                          . q{concat}
                          . $suffix{pdelly_call},
                        output_type     => q{b},
                        stderrfile_path => $program_info_path
                          . $UNDERSCORE
                          . $sv_type
                          . $UNDERSCORE
                          . q{concat.stderr.txt},
                        allow_overlaps => 1,
                        rm_dups        => q{all},
                        FILEHANDLE     => $FILEHANDLE,
                    }
                );
                say {$FILEHANDLE} "\n";

                bcftools_index(
                    {
                        infile_path => $outfile_path_prefix
                          . $UNDERSCORE
                          . $sv_type
                          . $UNDERSCORE
                          . q{concat}
                          . $suffix{pdelly_call},
                        stderrfile_path => $xargs_file_path_prefix
                          . $DOT
                          . $sv_type
                          . $UNDERSCORE
                          . q{index.stderr.txt},
                        output_type => q{csi},
                        FILEHANDLE  => $FILEHANDLE,
                    }
                );
                say {$FILEHANDLE} $NEWLINE;
            }
            ### Filter calls
            say {$FILEHANDLE} q{## Delly filter};

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
          SV_TYPE:
            foreach my $sv_type ( @{ $active_parameter_href->{delly_types} } ) {

                if ( $sv_type ne q{TRA} ) {

                    delly_filter(
                        {
                            infile_path => $outfile_path_prefix
                              . $UNDERSCORE
                              . $sv_type
                              . $UNDERSCORE
                              . q{concat}
                              . $suffix{pdelly_call},
                            outfile_path => $outfile_path_prefix
                              . $UNDERSCORE
                              . $sv_type
                              . $UNDERSCORE
                              . q{filtered}
                              . $suffix{pdelly_call},
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
                            sv_type     => $sv_type,
                            filter_mode => q{germline},
                            FILEHANDLE  => $XARGSFILEHANDLE,
                        }
                    );
                    say {$XARGSFILEHANDLE} $NEWLINE;
                }
                else {

                    delly_filter(
                        {
                            infile_path => $outfile_path_prefix
                              . $UNDERSCORE
                              . $sv_type
                              . $suffix{pdelly_call},
                            outfile_path => $outfile_path_prefix
                              . $UNDERSCORE
                              . $sv_type
                              . $UNDERSCORE
                              . q{filtered}
                              . $suffix{pdelly_call},
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
                            sv_type     => $sv_type,
                            filter_mode => q{germline},
                            FILEHANDLE  => $XARGSFILEHANDLE,
                        }
                    );
                    say {$XARGSFILEHANDLE} $NEWLINE;
                }
            }

        }

        ## Assemble filepaths
        my @file_paths;

        ### Concatenate SV types
        if ( scalar( @{ $active_parameter_href->{sample_ids} } ) > 1 ) {

            say {$FILEHANDLE} q{## bcftools concat - merge all SV types};

            @file_paths = map {
                    $outfile_path_prefix
                  . $UNDERSCORE
                  . $_
                  . $UNDERSCORE
                  . q{filtered}
                  . $suffix{pdelly_call}
            } @{ $active_parameter_href->{delly_types} };
        }
        else {    #Only one sample

            say {$FILEHANDLE}
              q{## Only one sample - skip merging and regenotyping};
            say {$FILEHANDLE}
              q{## bcftools concat - merge all SV types and contigs};

          SV_TYPE:
            foreach my $sv_type ( @{ $active_parameter_href->{delly_types} } ) {

                if ( $sv_type ne q{TRA} ) {

                  CONTIG:
                    foreach my $contig (@contigs) {

                        ## Assemble file paths by adding file ending
                        push @file_paths, map {
                                $infile_path_prefix{$_}{pdelly_call}
                              . $UNDERSCORE
                              . $contig
                              . $UNDERSCORE
                              . $sv_type
                              . $suffix{pdelly_call}
                        } @{ $active_parameter_href->{sample_ids} };
                    }
                }
                else {

                    push @file_paths, map {
                            $infile_path_prefix{$_}{pdelly_call}
                          . $UNDERSCORE
                          . $sv_type
                          . $suffix{pdelly_call}
                    } @{ $active_parameter_href->{sample_ids} };
                }
            }
        }
        bcftools_concat(
            {
                infile_paths_ref => \@file_paths,
                outfile_path     => $outfile_path_prefix
                  . $UNDERSCORE
                  . q{concat}
                  . $outfile_suffix,
                output_type     => q{v},
                stderrfile_path => $program_info_path
                  . $UNDERSCORE
                  . q{concat.stderr.txt},
                allow_overlaps => 1,
                rm_dups        => q{all},
                FILEHANDLE     => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Writes sbatch code to supplied filehandle to sort variants in vcf format
        say {$FILEHANDLE} q{## Picard SortVcf};
        picardtools_sortvcf(
            {
                FILEHANDLE       => $FILEHANDLE,
                infile_paths_ref => [
                        $outfile_path_prefix
                      . $UNDERSCORE
                      . q{concat}
                      . $outfile_suffix
                ],
                java_jar => catfile(
                    $active_parameter_href->{picardtools_path},
                    q{picard.jar}
                ),
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                outfile_path        => $outfile_path_prefix . $outfile_suffix,
                memory_allocation   => q{Xmx2g},
                referencefile_path  => $referencefile_path,
                sequence_dictionary => catfile(
                    $reference_dir,
                    $file_info_href->{human_genome_reference_name_prefix}
                      . $DOT . q{dict}
                ),
                temp_directory => $temp_directory,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Copies file from temporary directory.
        say {$FILEHANDLE} $NEWLINE . q{## Copy file from temporary directory};
        migrate_file(
            {
                infile_path  => $outfile_path_prefix . $outfile_suffix,
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} q{wait}, $NEWLINE;

        if ( $active_parameter_href->{$mip_program_name} == 1 ) {

            add_program_outfile_to_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    program_name     => q{delly},
                    outdirectory     => $outfamily_directory,
                    outfile          => $outfile_prefix . $outfile_suffix,
                }
            );

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

        close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});
        close $XARGSFILEHANDLE
          or $log->logcroak(q{Could not close $XARGSFILEHANDLE});

    }

    return;
}

1;
