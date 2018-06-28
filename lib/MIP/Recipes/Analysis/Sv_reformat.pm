package MIP::Recipes::Analysis::Sv_reformat;

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
use List::MoreUtils qw{ first_index };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_sv_reformat };

}

## Constants
Readonly my $AMPERSAND  => q{&};
Readonly my $ASTERIX    => q{*};
Readonly my $DOT        => q{.};
Readonly my $EMPTY_STR  => q{};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_sv_reformat {

## Function : Concatenate and sort contig files. Optionally remove variants from genelist
## Returns  : |$xargs_file_counter
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => Variant call type
##          : $family_id               => Family id
##          : $FILEHANDLE              => Sbatch filehandle to write to
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $reference_dir           => MIP reference directory
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $sample_info_href;
    my $parameter_href;
    my $program_name;

    ## Default(s)
    my $call_type;
    my $family_id;
    my $outaligner_dir;
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
        call_type => {
            default     => q{SV},
            store       => \$call_type,
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
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir,
            strict_type => 1,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir,
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

    use MIP::Delete::List qw{ delete_male_contig };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Gnu::Software::Gnu_grep qw{ gnu_grep };
    use MIP::IO::Files qw{ migrate_files };
    use MIP::Program::Utility::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Picardtools qw{ sort_vcf };
    use MIP::QC::Record
      qw{ add_most_complete_vcf add_program_metafile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constant
    Readonly my $VCFPARSER_OUTFILE_COUNT =>
      $active_parameter_href->{vcfparser_outfile_count} - 1;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $program_mode = $active_parameter_href->{$program_name};

    ## Unpack parameters
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $job_id_chain = $parameter_href->{$program_name}{chain};
    my ( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name      => $program_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

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
            program_directory               => catfile($outaligner_dir),
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir );
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir );
    my $outfamily_file_directory =
      catfile( $active_parameter_href->{outdata_dir}, $family_id );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$family_id}{sv_rankvariant}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$program_name}{file_tag};

    ## Files
    my $infile_prefix  = $family_id . $infile_tag . $call_type;
    my $outfile_prefix = $family_id . $outfile_tag . $call_type;

    ## Paths
    my $file_path_prefix = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix =
      catfile( $temp_directory, $family_id . $outfile_tag . $call_type );
    my $final_path_prefix = catfile( $outfamily_directory, $outfile_prefix );

    ## Assign suffix
    my $file_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
            jobid_chain    => $job_id_chain,
        }
    );

    ## Set analysis types and infiles
    my @vcfparser_analysis_types = ( $EMPTY_STR, q{.selected} );
    splice @vcfparser_analysis_types, $VCFPARSER_OUTFILE_COUNT + 1;
    my @vcfparser_infiles =
      map { $infile_prefix . $_ . $file_suffix } @vcfparser_analysis_types;

    ## Copy file(s) to temporary directory
    migrate_files(
        {
            core_number  => scalar @vcfparser_infiles,
            FILEHANDLE   => $FILEHANDLE,
            infiles_ref  => \@vcfparser_infiles,
            indirectory  => $infamily_directory,
            outfile_path => $temp_directory,
        }
    );

    ## Sort vcf;
  VCFPARSER_ANLAYSIS_TYPE:
    foreach my $vcfparser_analysis_type (@vcfparser_analysis_types) {

        ## Writes sbatch code to supplied filehandle to sort variants in vcf format
        sort_vcf(
            {
                active_parameter_href => $active_parameter_href,
                FILEHANDLE            => $FILEHANDLE,
                sequence_dict_file    => catfile(
                    $reference_dir,
                    $file_info_href->{human_genome_reference_name_prefix}
                      . $DOT . q{dict}
                ),
                infile_paths_ref => [
                    $file_path_prefix . $vcfparser_analysis_type . $file_suffix
                ],
                outfile => $outfile_path_prefix
                  . $vcfparser_analysis_type
                  . $file_suffix,
            }
        );
        print {$FILEHANDLE} $AMPERSAND . $NEWLINE;
    }
    say {$FILEHANDLE} q{wait} . $NEWLINE;

    ## Remove variants in hgnc_id list from vcf
    if ( $active_parameter_href->{sv_reformat_remove_genes_file} ) {
      VCFPARSER_ANLAYSIS_TYPE:
        foreach my $vcfparser_analysis_type (@vcfparser_analysis_types) {

            ## Removes contig_names from contigs array if no male or other found
            gnu_grep(
                {
                    FILEHANDLE       => $FILEHANDLE,
                    filter_file_path => catfile(
                        $reference_dir,
                        $active_parameter_href->{sv_reformat_remove_genes_file}
                    ),
                    infile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $file_suffix,
                    invert_match => 1,
                    outfile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $UNDERSCORE
                      . q{filtered}
                      . $file_suffix,
                }
            );
            say {$FILEHANDLE} $AMPERSAND . $NEWLINE;

            ## Save filtered file
            my $sv_reformat_remove_genes_path =
                $final_path_prefix
              . $vcfparser_analysis_type
              . $UNDERSCORE
              . q{filtered}
              . $file_suffix;

            if ( $vcfparser_analysis_type eq $vcfparser_analysis_types[1] ) {

                ## Save filtered file
                add_program_metafile_to_sample_info(
                    {
                        program_name => $program_name,
                        metafile_tag =>
                          q{sv_reformat_remove_genes_file_clinical},
                        path             => $sv_reformat_remove_genes_path,
                        sample_info_href => $sample_info_href,
                    }
                );
            }
            else {

                ## Save filtered file
                add_program_metafile_to_sample_info(
                    {
                        program_name => $program_name,
                        metafile_tag =>
                          q{sv_reformat_remove_genes_file_research},
                        path             => $sv_reformat_remove_genes_path,
                        sample_info_href => $sample_info_href,
                    }
                );
            }
        }
        say {$FILEHANDLE} q{wait};

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        my @filtered_files =
          map { $infile_prefix . $_ . $UNDERSCORE . q{filtered} . $file_suffix }
          @vcfparser_analysis_types;
        migrate_files(
            {
                FILEHANDLE   => $FILEHANDLE,
                core_number  => scalar @filtered_files,
                indirectory  => $temp_directory,
                infiles_ref  => \@filtered_files,
                outfile_path => $outfamily_directory,
            }
        );
    }

    if ( $active_parameter_href->{sv_rankvariant_binary_file} ) {

        say {$FILEHANDLE} q{## Compress};
      VCFPARSER_ANLAYSIS_TYPE:
        foreach my $vcfparser_analysis_type (@vcfparser_analysis_types) {
            ## Compress or decompress original file or stream to outfile (if supplied)
            htslib_bgzip(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $file_suffix,
                    stdoutfile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $file_suffix
                      . $DOT . q{gz},
                    write_to_stdout => 1,
                }
            );
            say {$FILEHANDLE} $AMPERSAND;
        }
        say {$FILEHANDLE} q{wait} . $NEWLINE;

        say {$FILEHANDLE} q{## Index};
      VCFPARSER_ANLAYSIS_TYPE:
        foreach my $vcfparser_analysis_type (@vcfparser_analysis_types) {
            ## Index file using tabix
            htslib_tabix(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    force       => 1,
                    infile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $file_suffix
                      . $DOT . q{gz},
                    preset => substr $file_suffix,
                    1,
                }
            );
            say {$FILEHANDLE} $AMPERSAND;
        }
        say {$FILEHANDLE} q{wait} . $NEWLINE;
    }

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file(s) from temporary directory};
    my @outfiles = map { $outfile_prefix . $_ . $file_suffix . $ASTERIX }
      @vcfparser_analysis_types;
    migrate_files(
        {
            FILEHANDLE   => $FILEHANDLE,
            core_number  => scalar @outfiles,
            indirectory  => $temp_directory,
            infiles_ref  => \@outfiles,
            outfile_path => $outfamily_directory,
        }
    );

  VCFPARSER_ANLAYSIS_TYPE:
    foreach my $vcfparser_analysis_type (@vcfparser_analysis_types) {

        ## Adds the most complete vcf file to sample_info
        my $vcfparser_outfile_counter =
          first_index { $_ eq $vcfparser_analysis_type }
        @vcfparser_analysis_types;
        add_most_complete_vcf(
            {
                active_parameter_href => $active_parameter_href,
                path                  => $final_path_prefix
                  . $vcfparser_analysis_type
                  . $file_suffix,
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
                vcf_file_key     => q{sv}
                  . $UNDERSCORE
                  . substr( $file_suffix, 1 )
                  . $UNDERSCORE . q{file},
                vcfparser_outfile_counter => $vcfparser_outfile_counter,
            }
        );

        if ( $program_mode == 1 ) {

            if ( $vcfparser_outfile_counter == 1 ) {

                # Save clinical candidate list path
                my $clinical_candidate_path =
                  $final_path_prefix . $vcfparser_analysis_type . $file_suffix;
                add_program_metafile_to_sample_info(
                    {
                        metafile_tag     => q{clinical},
                        path             => $clinical_candidate_path,
                        program_name     => $program_name,
                        sample_info_href => $sample_info_href,
                    }
                );

                if ( $active_parameter_href->{sv_rankvariant_binary_file} ) {

                    my $sv_rankvariant_binary_file_path =
                        $final_path_prefix
                      . $vcfparser_analysis_type
                      . $file_suffix
                      . $DOT . q{gz};
                    add_most_complete_vcf(
                        {
                            active_parameter_href => $active_parameter_href,
                            path         => $sv_rankvariant_binary_file_path,
                            program_name => $program_name,
                            sample_info_href => $sample_info_href,
                            vcf_file_key     => q{sv_vcf_binary_file},
                            vcfparser_outfile_counter =>
                              $vcfparser_outfile_counter,
                        }
                    );
                }
            }
            else {

                # Save research candidate list path
                my $research_candidate_path =
                  $final_path_prefix . $vcfparser_analysis_type . $file_suffix;
                add_program_metafile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        program_name     => $program_name,
                        metafile_tag     => q{research},
                        path             => $research_candidate_path,
                    }
                );

                if ( $active_parameter_href->{sv_rankvariant_binary_file} ) {

                    my $sv_rankvariant_binary_file_path =
                        $final_path_prefix
                      . $vcfparser_analysis_type
                      . $file_suffix
                      . $DOT . q{gz};
                    add_most_complete_vcf(
                        {
                            active_parameter_href => $active_parameter_href,
                            path         => $sv_rankvariant_binary_file_path,
                            program_name => $program_name,
                            sample_info_href => $sample_info_href,
                            vcf_file_key     => q{sv_vcf_binary_file},
                            vcfparser_outfile_counter =>
                              $vcfparser_outfile_counter,
                        }
                    );
                }
            }
        }
    }
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

1;
