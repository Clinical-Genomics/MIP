package MIP::Recipes::Analysis::Endvariantannotationblock;

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
    our @EXPORT_OK =
      qw{ analysis_endvariantannotationblock analysis_endvariantannotationblock_rio };

}

## Constants
Readonly my $ASTERIX    => q{*};
Readonly my $DOT        => q{.};
Readonly my $EMPTY_STR  => q{};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_endvariantannotationblock {

## Function : Concatenate ouput from variant annotation block.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => Variant call type
##          : $family_id               => Family id
##          : $file_info_href          => File info hash {REF}
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $parameter_href          => Parameter hash {REF}
##          : $program_info_path       => The program info path
##          : $program_name            => Program name
##          : $reference_dir           => MIP reference directory
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_path;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_info_path;
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
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        file_path      => { strict_type => 1, store => \$file_path },
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
        reference_dir => {
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

    use MIP::Delete::List qw{ delete_contig_elements };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Gnu::Software::Gnu_grep qw{ gnu_grep };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Utility::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Program::Variantcalling::Gatk qw{ gatk_concatenate_variants };
    use MIP::QC::Record
      qw{ add_most_complete_vcf add_program_metafile_to_sample_info };
    use MIP::Set::File qw{ set_file_suffix };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constant
    Readonly my $VCFPARSER_OUTFILE_COUNT =>
      $active_parameter_href->{vcfparser_outfile_count} - 1;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Unpack parameters
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $job_id_chain  = $parameter_href->{$mip_program_name}{chain};
    my $reduce_io_ref = \$active_parameter_href->{reduce_io};
    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    my $vcfparser_analysis_type = $EMPTY_STR;

    ## Set default contigs
    my $contigs_size_ordered_ref =
      \@{ $file_info_href->{contigs_size_ordered} };
    my @contigs = @{ $file_info_href->{contigs} };

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            FILEHANDLE                      => $FILEHANDLE,
            directory_id                    => $family_id,
            job_id_href                     => $job_id_href,
            process_time                    => $time,
            program_directory               => catfile($outaligner_dir),
            program_name                    => $program_name,
            source_environment_commands_ref => [$source_environment_cmd],
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
    my $infile_tag = $file_info_href->{$family_id}{prankvariant}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{prankvariant}{file_tag};

    ## Files
    my $infile_prefix = $family_id . $infile_tag . $call_type;
    my $file_path_prefix = catfile( $temp_directory, $infile_prefix );

    ## Paths
    my $outfile_prefix = $family_id . $outfile_tag . $call_type;
    my $outfile_path_prefix =
      catfile( $temp_directory, $family_id . $outfile_tag . $call_type );
    my $final_path_prefix = catfile( $outfamily_directory, $outfile_prefix );

    ## Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            jobid_chain    => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

    my $outfile_suffix = set_file_suffix(
        {
            file_suffix => $parameter_href->{$mip_program_name}{outfile_suffix},
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

    ## Determined by vcfparser output
    for my $vcfparser_outfile_counter ( 0 .. $VCFPARSER_OUTFILE_COUNT ) {

        if ( $vcfparser_outfile_counter == 1 ) {

            ## SelectFile variants
            $vcfparser_analysis_type = $DOT . q{selected};

            ## Selectfile contigs
            $contigs_size_ordered_ref =
              \@{ $file_info_href->{sorted_select_file_contigs} };
            @contigs = @{ $file_info_href->{select_file_contigs} };

            ## Remove MT|M since no exome kit so far has mitochondrial probes
            if ( $consensus_analysis_type eq q{wes} ) {

                ## Removes an element from array and return new array while leaving orginal elements_ref untouched
                @contigs = delete_contig_elements(
                    {
                        elements_ref =>
                          \@{ $file_info_href->{select_file_contigs} },
                        remove_contigs_ref => [qw{ MT M }],
                    }
                );
            }
        }

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        $xargs_file_counter = xargs_migrate_contig_files(
            {
                contigs_ref        => $contigs_size_ordered_ref,
                core_number        => $core_number,
                FILEHANDLE         => $FILEHANDLE,
                file_path          => $file_path,
                infile             => $infile_prefix,
                program_info_path  => $program_info_path,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                xargs_file_counter => $xargs_file_counter,
                file_ending        => $vcfparser_analysis_type
                  . $infile_suffix
                  . $ASTERIX,
                indirectory    => $infamily_directory,
                temp_directory => $active_parameter_href->{temp_directory},
            }
        );

        ## Writes sbatch code to supplied filehandle to concatenate variants in vcf format. Each array element is combined with the infile prefix and postfix.
        gatk_concatenate_variants(
            {
                active_parameter_href => $active_parameter_href,
                elements_ref          => \@contigs,
                FILEHANDLE            => $FILEHANDLE,
                infile_prefix         => $file_path_prefix . $UNDERSCORE,
                infile_postfix => $vcfparser_analysis_type . $infile_suffix,
                outfile_path_prefix => $outfile_path_prefix
                  . $vcfparser_analysis_type,
                outfile_suffix => $outfile_suffix,
            }
        );

        ## Remove variants in hgnc_id list from vcf
        if ( $active_parameter_href
            ->{endvariantannotationblock_remove_genes_file} )
        {

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
                      . $outfile_suffix,
                    outfile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $UNDERSCORE
                      . q{filtered}
                      . $outfile_suffix,
                    invert_match => 1,
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            if ( $vcfparser_outfile_counter == 1 ) {

                ## Save filtered file
                $sample_info_href->{program}{$program_name}
                  {reformat_remove_genes_file}{clinical}{path} =
                    $final_path_prefix
                  . $vcfparser_analysis_type
                  . $UNDERSCORE
                  . q{filtered}
                  . $outfile_suffix;
            }
            else {

                ## Save filtered file
                $sample_info_href->{program}{$program_name}
                  {reformat_remove_genes_file}{research}{path} =
                    $final_path_prefix
                  . $vcfparser_analysis_type
                  . $UNDERSCORE
                  . q{filtered}
                  . $outfile_suffix;
            }

            ## Copies file from temporary directory.
            say {$FILEHANDLE} q{## Copy file from temporary directory};
            migrate_file(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . q{_filtered}
                      . $outfile_suffix,
                    outfile_path => $outfamily_directory,
                }
            );
            say {$FILEHANDLE} q{wait}, $NEWLINE;
        }

        if ( $active_parameter_href->{rankvariant_binary_file} ) {

            ## Compress or decompress original file or stream to outfile (if supplied)
            htslib_bgzip(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $outfile_suffix,
                    stdoutfile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $outfile_suffix
                      . $DOT . q{gz},
                    write_to_stdout => 1,
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            ## Index file using tabix
            htslib_tabix(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    force       => 1,
                    infile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $outfile_suffix
                      . $DOT . q{gz},
                    preset => substr $outfile_suffix,
                    1,
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                infile_path => $outfile_path_prefix
                  . $vcfparser_analysis_type
                  . $outfile_suffix
                  . $ASTERIX,
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} q{wait}, $NEWLINE;

        ## Adds the most complete vcf file to sample_info
        add_most_complete_vcf(
            {
                active_parameter_href => $active_parameter_href,
                path                  => $final_path_prefix
                  . $vcfparser_analysis_type
                  . $outfile_suffix,
                program_name              => $program_name,
                sample_info_href          => $sample_info_href,
                vcfparser_outfile_counter => $vcfparser_outfile_counter,
            }
        );

        if ( $mip_program_mode == 1 ) {

            if ( $vcfparser_outfile_counter == 1 ) {

                # Save clinical candidate list path
                my $clinical_candidate_path =
                    $final_path_prefix
                  . $vcfparser_analysis_type
                  . $outfile_suffix;
                add_program_metafile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        program_name     => $program_name,
                        metafile_tag     => q{clinical},
                        path             => $clinical_candidate_path,
                    }
                );

                if ( $active_parameter_href->{rankvariant_binary_file} ) {

                    $sample_info_href->{vcf_binary_file}{clinical}{path} =
                        $final_path_prefix
                      . $vcfparser_analysis_type
                      . $outfile_suffix
                      . $DOT . q{gz};
                }
            }
            else {

                # Save research candidate list path
                my $research_candidate_path =
                    $final_path_prefix
                  . $vcfparser_analysis_type
                  . $outfile_suffix;
                add_program_metafile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        program_name     => $program_name,
                        metafile_tag     => q{research},
                        path             => $research_candidate_path,
                    }
                );

                if ( $active_parameter_href->{rankvariant_binary_file} ) {

                    $sample_info_href->{vcf_binary_file}{research}{path} =
                        $final_path_prefix
                      . $vcfparser_analysis_type
                      . $outfile_suffix
                      . $DOT . q{gz};
                }
            }
        }
    }

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $mip_program_mode == 1 ) {

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

sub analysis_endvariantannotationblock_rio {

## Function : Concatenate ouput from variant annotation block.
## Returns  : $xargs_file_counter
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => Variant call type
##          : $family_id               => Family id
##          : $FILEHANDLE              => Sbatch filehandle to write to
##          : $file_info_href          => File info hash {REF}
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $parameter_href          => Parameter hash {REF}
##          : $program_info_path       => The program info path
##          : $program_name            => Program name
##          : $reference_dir           => MIP reference directory
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $FILEHANDLE;
    my $file_path;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_info_path;
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
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        FILEHANDLE => { store       => \$FILEHANDLE, },
        file_path  => { strict_type => 1, store => \$file_path },
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
        reference_dir => {
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

    use MIP::Delete::List qw{ delete_contig_elements };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Gnu::Software::Gnu_grep qw{ gnu_grep };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Utility::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Program::Variantcalling::Gatk qw{ gatk_concatenate_variants };
    use MIP::QC::Record
      qw{ add_most_complete_vcf add_program_metafile_to_sample_info };
    use MIP::Set::File qw{ set_file_suffix };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constant
    Readonly my $VCFPARSER_OUTFILE_COUNT =>
      $active_parameter_href->{vcfparser_outfile_count} - 1;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Unpack parameters
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $job_id_chain  = $parameter_href->{$mip_program_name}{chain};
    my $reduce_io_ref = \$active_parameter_href->{reduce_io};
    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    my $vcfparser_analysis_type = $EMPTY_STR;

    ## Set default contigs
    my $contigs_size_ordered_ref =
      \@{ $file_info_href->{contigs_size_ordered} };
    my @contigs = @{ $file_info_href->{contigs} };

    ## Filehandles
    # Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir );
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir );
    my $outfamily_file_directory =
      catfile( $active_parameter_href->{outdata_dir}, $family_id );

    ## Assign file_tags
    my $infile_tag = $file_info_href->{$family_id}{prankvariant}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{prankvariant}{file_tag};

    ## Files
    my $infile_prefix = $family_id . $infile_tag . $call_type;
    my $file_path_prefix = catfile( $temp_directory, $infile_prefix );

    ## Paths
    my $outfile_prefix = $family_id . $outfile_tag . $call_type;
    my $outfile_path_prefix =
      catfile( $temp_directory, $family_id . $outfile_tag . $call_type );
    my $final_path_prefix = catfile( $outfamily_directory, $outfile_prefix );

    ## Assign suffix
    my $infile_suffix = get_file_suffix(
        {
            jobid_chain    => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

    my $outfile_suffix = set_file_suffix(
        {
            file_suffix => $parameter_href->{$mip_program_name}{outfile_suffix},
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

    ## Determined by vcfparser output
    for my $vcfparser_outfile_counter ( 0 .. $VCFPARSER_OUTFILE_COUNT ) {

        if ( $vcfparser_outfile_counter == 1 ) {

            ## SelectFile variants
            $vcfparser_analysis_type = $DOT . q{selected};

            ## Selectfile contigs
            $contigs_size_ordered_ref =
              \@{ $file_info_href->{sorted_select_file_contigs} };
            @contigs = @{ $file_info_href->{select_file_contigs} };

            ## Remove MT|M since no exome kit so far has mitochondrial probes
            if ( $consensus_analysis_type eq q{wes} ) {

                ## Removes an element from array and return new array while leaving orginal elements_ref untouched
                @contigs = delete_contig_elements(
                    {
                        elements_ref =>
                          \@{ $file_info_href->{select_file_contigs} },
                        remove_contigs_ref => [qw{ MT M }],
                    }
                );
            }
        }

        ## Writes sbatch code to supplied filehandle to concatenate variants in vcf format. Each array element is combined with the infile prefix and postfix.
        gatk_concatenate_variants(
            {
                active_parameter_href => $active_parameter_href,
                elements_ref          => \@contigs,
                FILEHANDLE            => $FILEHANDLE,
                infile_prefix         => $file_path_prefix . $UNDERSCORE,
                infile_postfix => $vcfparser_analysis_type . $infile_suffix,
                outfile_path_prefix => $outfile_path_prefix
                  . $vcfparser_analysis_type,
                outfile_suffix => $outfile_suffix,
            }
        );

        ## Remove variants in hgnc_id list from vcf
        if ( $active_parameter_href
            ->{endvariantannotationblock_remove_genes_file} )
        {

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
                      . $outfile_suffix,
                    outfile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $UNDERSCORE
                      . q{filtered}
                      . $outfile_suffix,
                    invert_match => 1,
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            if ( $vcfparser_outfile_counter == 1 ) {

                ## Save filtered file
                $sample_info_href->{program}{$program_name}
                  {reformat_remove_genes_file}{clinical}{path} =
                    $final_path_prefix
                  . $vcfparser_analysis_type
                  . $UNDERSCORE
                  . q{filtered}
                  . $outfile_suffix;
            }
            else {

                ## Save filtered file
                $sample_info_href->{program}{$program_name}
                  {reformat_remove_genes_file}{research}{path} =
                    $final_path_prefix
                  . $vcfparser_analysis_type
                  . $UNDERSCORE
                  . q{filtered}
                  . $outfile_suffix;
            }

            ## Copies file from temporary directory.
            say {$FILEHANDLE} q{## Copy file from temporary directory};
            migrate_file(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . q{_filtered}
                      . $outfile_suffix,
                    outfile_path => $outfamily_directory,
                }
            );
            say {$FILEHANDLE} q{wait}, $NEWLINE;
        }

        if ( $active_parameter_href->{rankvariant_binary_file} ) {

            ## Compress or decompress original file or stream to outfile (if supplied)
            htslib_bgzip(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $outfile_suffix,
                    stdoutfile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $outfile_suffix
                      . $DOT . q{gz},
                    write_to_stdout => 1,
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            ## Index file using tabix
            htslib_tabix(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    force       => 1,
                    infile_path => $outfile_path_prefix
                      . $vcfparser_analysis_type
                      . $outfile_suffix
                      . $DOT . q{gz},
                    preset => substr $outfile_suffix,
                    1,
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                infile_path => $outfile_path_prefix
                  . $vcfparser_analysis_type
                  . $outfile_suffix
                  . $ASTERIX,
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} q{wait}, $NEWLINE;

        ## Adds the most complete vcf file to sample_info
        add_most_complete_vcf(
            {
                active_parameter_href => $active_parameter_href,
                path                  => $final_path_prefix
                  . $vcfparser_analysis_type
                  . $outfile_suffix,
                program_name              => $program_name,
                sample_info_href          => $sample_info_href,
                vcfparser_outfile_counter => $vcfparser_outfile_counter,
            }
        );

        if ( $mip_program_mode == 1 ) {

            if ( $vcfparser_outfile_counter == 1 ) {

                # Save clinical candidate list path
                my $clinical_candidate_path =
                    $final_path_prefix
                  . $vcfparser_analysis_type
                  . $outfile_suffix;
                add_program_metafile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        program_name     => $program_name,
                        metafile_tag     => q{clinical},
                        path             => $clinical_candidate_path,
                    }
                );

                if ( $active_parameter_href->{rankvariant_binary_file} ) {

                    $sample_info_href->{vcf_binary_file}{clinical}{path} =
                        $final_path_prefix
                      . $vcfparser_analysis_type
                      . $outfile_suffix
                      . $DOT . q{gz};
                }
            }
            else {

                # Save research candidate list path
                my $research_candidate_path =
                    $final_path_prefix
                  . $vcfparser_analysis_type
                  . $outfile_suffix;
                add_program_metafile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        program_name     => $program_name,
                        metafile_tag     => q{research},
                        path             => $research_candidate_path,
                    }
                );

                if ( $active_parameter_href->{rankvariant_binary_file} ) {

                    $sample_info_href->{vcf_binary_file}{research}{path} =
                        $final_path_prefix
                      . $vcfparser_analysis_type
                      . $outfile_suffix
                      . $DOT . q{gz};
                }
            }
        }
    }

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $mip_program_mode == 1 ) {

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
    ## Track the number of created xargs scripts per module for Block algorithm
    return $xargs_file_counter;
}

1;
