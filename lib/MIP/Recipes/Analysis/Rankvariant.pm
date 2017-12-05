package MIP::Recipes::Analysis::Rankvariant;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile devnull };
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
      qw{ analysis_rankvariant analysis_rankvariant_rio analysis_rankvariant_rio_unaffected analysis_rankvariant_unaffected };

}

## Constants
Readonly my $ASTERIX    => q{*};
Readonly my $DASH       => q{-};
Readonly my $DOT        => q{.};
Readonly my $EMPTY_STR  => q{};
Readonly my $NEWLINE    => qq{\n};
Readonly my $PIPE       => q{|};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_rankvariant {

## Function : Annotate and score variants depending on mendelian inheritance, frequency and phenotype etc.
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
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $file_path;
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
        file_path               => { strict_type => 1, store => \$file_path },
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

    use MIP::Check::Cluster qw{ check_max_core_number };
    use MIP::Cluster qw{ get_core_number };
    use MIP::Delete::List qw{ delete_contig_elements };
    use MIP::File::Format::Pedigree qw{ create_fam_file };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Genmod
      qw{ genmod_annotate genmod_compound genmod_models genmod_score };
    use MIP::QC::Record
      qw{ add_program_metafile_to_sample_info add_program_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Set::File qw{ set_file_suffix };
    use MIP::Script::Setup_script qw(setup_script);

    ## Constant
    Readonly my $CORE_NUMBER_REQUESTED => 16;
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
    my $job_id_chain            = $parameter_href->{$mip_program_name}{chain};
    my $reduce_io_ref           = \$active_parameter_href->{reduce_io};
    my $vcfparser_analysis_type = $EMPTY_STR;
    my $xargs_file_path_prefix;
    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    ## Set default contigs
    my @contigs_size_ordered = @{ $file_info_href->{contigs_size_ordered} };
    my @contigs              = @{ $file_info_href->{contigs} };

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            module_core_number   => $core_number,
            modifier_core_number => scalar(@contigs),
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    ### Detect the number of cores to use per genmod process.
    ## Limit number of cores requested to the maximum number of cores available per node
    my $genmod_core_number = check_max_core_number(
        {
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
            core_number_requested => $CORE_NUMBER_REQUESTED,
        }
    );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $family_id,
            FILEHANDLE                      => $FILEHANDLE,
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
    my $infile_tag = $file_info_href->{$family_id}{psnpeff}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$mip_program_name}{file_tag};

    ## Files
    my $infile_prefix  = $family_id . $infile_tag . $call_type;
    my $outfile_prefix = $family_id . $outfile_tag . $call_type;

    ## Paths
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ## Assign suffix
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

    my $family_file =
      catfile( $outfamily_file_directory, $family_id . $DOT . q{fam} );

    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            active_parameter_href => $active_parameter_href,
            fam_file_path         => $family_file,
            FILEHANDLE            => $FILEHANDLE,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Determined by vcfparser output
    for my $vcfparser_outfile_counter ( 0 .. $VCFPARSER_OUTFILE_COUNT ) {

        if ( $vcfparser_outfile_counter == 1 ) {

            ## SelectFile variants
            $vcfparser_analysis_type = $DOT . q{selected};
            ## Selectfile contigs
            @contigs_size_ordered =
              @{ $file_info_href->{sorted_select_file_contigs} };
            @contigs = @{ $file_info_href->{select_file_contigs} };

            ## Remove MT|M since no exome kit so far has mitochondrial probes
            if ( $consensus_analysis_type eq q{wes} ) {

                ## Removes an element from array and return new array while leaving orginal elements_ref untouched
                @contigs = delete_contig_elements(
                    {
                        elements_ref =>
                          \@{ $file_info_href->{select_file_contigs} },
                        remove_contigs_ref => [q{ MT M }],
                    }
                );

                ## Removes an element from array and return new array while leaving orginal elements_ref untouched
                @contigs_size_ordered = delete_contig_elements(
                    {
                        elements_ref =>
                          \@{ $file_info_href->{sorted_select_file_contigs} },
                        remove_contigs_ref => [qw{ MT M }],
                    }
                );
            }
        }

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        $xargs_file_counter = xargs_migrate_contig_files(
            {
                contigs_ref => \@contigs_size_ordered,
                core_number => $core_number,
                FILEHANDLE  => $FILEHANDLE,
                file_ending => $vcfparser_analysis_type
                  . $infile_suffix
                  . $ASTERIX,
                file_path          => $file_path,
                indirectory        => $infamily_directory,
                infile             => $infile_prefix,
                program_info_path  => $program_info_path,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                temp_directory     => $active_parameter_href->{temp_directory},
                xargs_file_counter => $xargs_file_counter,
            }
        );

        ## Genmod
        say {$FILEHANDLE} q{## Genmod};

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                core_number        => $genmod_core_number,
                FILEHANDLE         => $FILEHANDLE,
                file_path          => $file_path,
                program_info_path  => $program_info_path,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                xargs_file_counter => $xargs_file_counter,
            }
        );

        ## Track which genmod modules has been processed
        my $genmod_module = $EMPTY_STR;

        ## Process per contig
        while ( my ( $contig_index, $contig ) = each @contigs_size_ordered ) {

            ## Get parameters
            # Restart for next contig
            $genmod_module = $EMPTY_STR;

            ## Infile
            my $genmod_indata =
                $file_path_prefix
              . $UNDERSCORE
              . $contig
              . $vcfparser_analysis_type
              . $infile_suffix;

            ## Output stream
            my $genmod_outfile_path =
              catfile( dirname( devnull() ), q{stdout} );

            ## Genmod Annotate
            $genmod_module = q{_annotate};

            genmod_annotate(
                {
                    annotate_region =>
                      $active_parameter_href->{genmod_annotate_regions},
                    cadd_file_paths_ref =>
                      \@{ $active_parameter_href->{genmod_annotate_cadd_files}
                      },
                    FILEHANDLE   => $XARGSFILEHANDLE,
                    infile_path  => $genmod_indata,
                    outfile_path => $genmod_outfile_path,
                    spidex_file_path =>
                      $active_parameter_href->{genmod_annotate_spidex_file},
                    stderrfile_path => $xargs_file_path_prefix
                      . $DOT
                      . $contig
                      . $genmod_module
                      . $DOT
                      . q{stderr.txt},
                    temp_directory_path => $temp_directory,
                    verbosity           => q{v},
                }
            );

            ## Pipe
            print {$XARGSFILEHANDLE} $PIPE . $SPACE;

            # Preparation for next module
            $genmod_indata = $DASH;

            ## Genmod Models
            $genmod_module .= q{_models};

            my $use_vep;
            ## Use VEP annotations in compound models
            if ( $active_parameter_href->{pvarianteffectpredictor}
                and not $active_parameter_href->{genmod_annotate_regions} )
            {

                $use_vep = 1;
            }
            genmod_models(
                {
                    FILEHANDLE  => $XARGSFILEHANDLE,
                    family_file => $family_file,
                    family_type =>
                      $active_parameter_href->{genmod_models_family_type},
                    infile_path  => $genmod_indata,
                    outfile_path => catfile( dirname( devnull() ), q{stdout} ),
                    reduced_penetrance_file_path => $active_parameter_href
                      ->{genmod_models_reduced_penetrance_file},
                    stderrfile_path => $xargs_file_path_prefix
                      . $DOT
                      . $contig
                      . $genmod_module
                      . $DOT
                      . q{stderr.txt},
                    temp_directory_path => $temp_directory,
                    thread_number       => 4,
                    vep                 => $use_vep,
                    verbosity           => q{v},
                    whole_gene =>
                      $active_parameter_href->{genmod_models_whole_gene},
                }
            );
            ## Pipe
            print {$XARGSFILEHANDLE} $PIPE . $SPACE;

            ## Genmod Score
            $genmod_module .= q{_score};
            genmod_score(
                {
                    family_file => $family_file,
                    family_type =>
                      $active_parameter_href->{genmod_models_family_type},
                    FILEHANDLE   => $XARGSFILEHANDLE,
                    infile_path  => $genmod_indata,
                    outfile_path => catfile( dirname( devnull() ), q{stdout} ),
                    rank_result  => 1,
                    rank_model_file_path =>
                      $active_parameter_href->{rank_model_file},
                    stderrfile_path => $xargs_file_path_prefix
                      . $DOT
                      . $contig
                      . $genmod_module
                      . $DOT
                      . q{stderr.txt},
                    verbosity => q{v},
                }
            );
            ## Pipe
            print {$XARGSFILEHANDLE} $PIPE . $SPACE;

            ##Genmod Compound
            $genmod_module .= q{_compound};

            genmod_compound(
                {
                    FILEHANDLE   => $XARGSFILEHANDLE,
                    infile_path  => $genmod_indata,
                    outfile_path => $outfile_path_prefix
                      . $UNDERSCORE
                      . $contig
                      . $vcfparser_analysis_type
                      . $infile_suffix,
                    stderrfile_path => $xargs_file_path_prefix
                      . $DOT
                      . $contig
                      . $genmod_module
                      . $DOT
                      . q{stderr.txt},
                    temp_directory_path => $temp_directory,
                    verbosity           => q{v},
                    vep                 => $use_vep,
                }
            );
            say {$XARGSFILEHANDLE} $NEWLINE;
        }

        ## Copies file from temporary directory. Per contig
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        ($xargs_file_counter) = xargs_migrate_contig_files(
            {
                contigs_ref => \@contigs_size_ordered,
                core_number => $active_parameter_href->{max_cores_per_node},
                FILEHANDLE  => $FILEHANDLE,
                file_path   => $file_path,
                file_ending => $vcfparser_analysis_type
                  . $outfile_suffix
                  . $ASTERIX,
                outdirectory       => $outfamily_directory,
                outfile            => $outfile_prefix,
                program_info_path  => $program_info_path,
                temp_directory     => $temp_directory,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                xargs_file_counter => $xargs_file_counter,
            }
        );
    }

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $mip_program_mode == 1 ) {

        my $qc_genmod_outfile =
            $outfile_prefix
          . $UNDERSCORE
          . $file_info_href->{contigs_size_ordered}[0]
          . $vcfparser_analysis_type
          . $outfile_suffix;
        add_program_outfile_to_sample_info(
            {
                outdirectory     => $outfamily_directory,
                outfile          => $qc_genmod_outfile,
                program_name     => q{genmod},
                sample_info_href => $sample_info_href,
            }
        );

        # Add to Sample_info
        if ( defined $active_parameter_href->{rank_model_file} ) {

            my $rank_model_version;
            if ( $active_parameter_href->{rank_model_file} =~
                / v(\d+[.]\d+[.]\d+ | \d+[.]\d+) /sxm )
            {

                $rank_model_version = $1;
            }
            add_program_metafile_to_sample_info(
                {
                    file =>
                      basename( $active_parameter_href->{rank_model_file} ),
                    metafile_tag => q{rank_model},
                    path         => $active_parameter_href->{rank_model_file},
                    program_name => q{genmod},
                    sample_info_href => $sample_info_href,
                    version          => $rank_model_version,
                }
            );
        }
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

sub analysis_rankvariant_rio {

## Function : Annotate and score variants depending on mendelian inheritance, frequency and phenotype etc.
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
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $FILEHANDLE;
    my $file_info_href;
    my $file_path;
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
        file_path               => { strict_type => 1, store => \$file_path },
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

    use MIP::Check::Cluster qw{ check_max_core_number };
    use MIP::Cluster qw{ get_core_number };
    use MIP::Delete::List qw{ delete_contig_elements };
    use MIP::File::Format::Pedigree qw{ create_fam_file };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Genmod
      qw{ genmod_annotate genmod_compound genmod_models genmod_score };
    use MIP::QC::Record
      qw{ add_program_metafile_to_sample_info add_program_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Set::File qw{ set_file_suffix };

    ## Constant
    Readonly my $CORE_NUMBER_REQUESTED => 16;
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
    my $job_id_chain            = $parameter_href->{$mip_program_name}{chain};
    my $reduce_io_ref           = \$active_parameter_href->{reduce_io};
    my $vcfparser_analysis_type = $EMPTY_STR;
    my $xargs_file_path_prefix;
    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    ## Set default contigs
    my @contigs_size_ordered = @{ $file_info_href->{contigs_size_ordered} };
    my @contigs              = @{ $file_info_href->{contigs} };

    ## Filehandles
    # Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            module_core_number   => $core_number,
            modifier_core_number => scalar(@contigs),
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    ### Detect the number of cores to use per genmod process.
    ## Limit number of cores requested to the maximum number of cores available per node
    my $genmod_core_number = check_max_core_number(
        {
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
            core_number_requested => $CORE_NUMBER_REQUESTED,
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
    my $infile_tag = $file_info_href->{$family_id}{psnpeff}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$mip_program_name}{file_tag};

    ## Files
    my $infile_prefix  = $family_id . $infile_tag . $call_type;
    my $outfile_prefix = $family_id . $outfile_tag . $call_type;

    ## Paths
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ## Assign suffix
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

    my $family_file =
      catfile( $outfamily_file_directory, $family_id . $DOT . q{fam} );

    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            active_parameter_href => $active_parameter_href,
            fam_file_path         => $family_file,
            FILEHANDLE            => $FILEHANDLE,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Determined by vcfparser output
    for my $vcfparser_outfile_counter ( 0 .. $VCFPARSER_OUTFILE_COUNT ) {

        if ( $vcfparser_outfile_counter == 1 ) {

            ## SelectFile variants
            $vcfparser_analysis_type = $DOT . q{selected};
            ## Selectfile contigs
            @contigs_size_ordered =
              @{ $file_info_href->{sorted_select_file_contigs} };
            @contigs = @{ $file_info_href->{select_file_contigs} };

            ## Remove MT|M since no exome kit so far has mitochondrial probes
            if ( $consensus_analysis_type eq q{wes} ) {

                ## Removes an element from array and return new array while leaving orginal elements_ref untouched
                @contigs = delete_contig_elements(
                    {
                        elements_ref =>
                          \@{ $file_info_href->{select_file_contigs} },
                        remove_contigs_ref => [q{ MT M }],
                    }
                );

                ## Removes an element from array and return new array while leaving orginal elements_ref untouched
                @contigs_size_ordered = delete_contig_elements(
                    {
                        elements_ref =>
                          \@{ $file_info_href->{sorted_select_file_contigs} },
                        remove_contigs_ref => [qw{ MT M }],
                    }
                );
            }
        }

        ## Genmod
        say {$FILEHANDLE} q{## Genmod};

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                core_number        => $genmod_core_number,
                FILEHANDLE         => $FILEHANDLE,
                file_path          => $file_path,
                program_info_path  => $program_info_path,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                xargs_file_counter => $xargs_file_counter,
            }
        );

        ## Track which genmod modules has been processed
        my $genmod_module = $EMPTY_STR;

        ## Process per contig
        while ( my ( $contig_index, $contig ) = each @contigs_size_ordered ) {

            ## Get parameters
            # Restart for next contig
            $genmod_module = $EMPTY_STR;

            ## InFile
            my $genmod_indata =
                $file_path_prefix
              . $UNDERSCORE
              . $contig
              . $vcfparser_analysis_type
              . $infile_suffix;

            ## Output stream
            my $genmod_outfile_path =
              catfile( dirname( devnull() ), q{stdout} );

            ## Genmod Annotate
            $genmod_module = q{_annotate};

            genmod_annotate(
                {
                    annotate_region =>
                      $active_parameter_href->{genmod_annotate_regions},
                    cadd_file_paths_ref =>
                      \@{ $active_parameter_href->{genmod_annotate_cadd_files}
                      },
                    FILEHANDLE   => $XARGSFILEHANDLE,
                    infile_path  => $genmod_indata,
                    outfile_path => $genmod_outfile_path,
                    spidex_file_path =>
                      $active_parameter_href->{genmod_annotate_spidex_file},
                    stderrfile_path => $xargs_file_path_prefix
                      . $DOT
                      . $contig
                      . $genmod_module
                      . $DOT
                      . q{stderr.txt},
                    temp_directory_path => $temp_directory,
                    verbosity           => q{v},
                }
            );

            ## Pipe
            print {$XARGSFILEHANDLE} $PIPE . $SPACE;

            # Preparation for next module
            $genmod_indata = $DASH;

            ## Genmod Models
            $genmod_module .= q{_models};

            my $use_vep;
            ## Use VEP annotations in compound models
            if ( $active_parameter_href->{pvarianteffectpredictor}
                and not $active_parameter_href->{genmod_annotate_regions} )
            {

                $use_vep = 1;
            }
            genmod_models(
                {
                    FILEHANDLE  => $XARGSFILEHANDLE,
                    family_file => $family_file,
                    family_type =>
                      $active_parameter_href->{genmod_models_family_type},
                    infile_path  => $genmod_indata,
                    outfile_path => catfile( dirname( devnull() ), q{stdout} ),
                    reduced_penetrance_file_path => $active_parameter_href
                      ->{genmod_models_reduced_penetrance_file},
                    stderrfile_path => $xargs_file_path_prefix
                      . $DOT
                      . $contig
                      . $genmod_module
                      . $DOT
                      . q{stderr.txt},
                    temp_directory_path => $temp_directory,
                    thread_number       => 4,
                    vep                 => $use_vep,
                    verbosity           => q{v},
                    whole_gene =>
                      $active_parameter_href->{genmod_models_whole_gene},
                }
            );
            ## Pipe
            print {$XARGSFILEHANDLE} $PIPE . $SPACE;

            ## Genmod Score
            $genmod_module .= q{_score};
            genmod_score(
                {
                    family_file => $family_file,
                    family_type =>
                      $active_parameter_href->{genmod_models_family_type},
                    FILEHANDLE   => $XARGSFILEHANDLE,
                    infile_path  => $genmod_indata,
                    outfile_path => catfile( dirname( devnull() ), q{stdout} ),
                    rank_result  => 1,
                    rank_model_file_path =>
                      $active_parameter_href->{rank_model_file},
                    stderrfile_path => $xargs_file_path_prefix
                      . $DOT
                      . $contig
                      . $genmod_module
                      . $DOT
                      . q{stderr.txt},
                    verbosity => q{v},
                }
            );
            ## Pipe
            print {$XARGSFILEHANDLE} $PIPE . $SPACE;

            ##Genmod Compound
            $genmod_module .= q{_compound};

            genmod_compound(
                {
                    FILEHANDLE   => $XARGSFILEHANDLE,
                    infile_path  => $genmod_indata,
                    outfile_path => $outfile_path_prefix
                      . $UNDERSCORE
                      . $contig
                      . $vcfparser_analysis_type
                      . $infile_suffix,
                    stderrfile_path => $xargs_file_path_prefix
                      . $DOT
                      . $contig
                      . $genmod_module
                      . $DOT
                      . q{stderr.txt},
                    temp_directory_path => $temp_directory,
                    verbosity           => q{v},
                    vep                 => $use_vep,
                }
            );
            say {$XARGSFILEHANDLE} $NEWLINE;
        }

        ## QC Data File(s)
        migrate_file(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $outfile_path_prefix
                  . $UNDERSCORE
                  . $file_info_href->{contigs_size_ordered}[0]
                  . $vcfparser_analysis_type
                  . $outfile_suffix,
                outfile_path => $outfamily_directory,
            }
        );
        say {$FILEHANDLE} q{wait}, $NEWLINE;
    }

    if ( $mip_program_mode == 1 ) {

        my $qc_genmod_outfile =
            $outfile_prefix
          . $UNDERSCORE
          . $file_info_href->{contigs_size_ordered}[0]
          . $vcfparser_analysis_type
          . $outfile_suffix;
        add_program_outfile_to_sample_info(
            {
                outdirectory     => $outfamily_directory,
                outfile          => $qc_genmod_outfile,
                program_name     => q{genmod},
                sample_info_href => $sample_info_href,
            }
        );

        # Add to Sample_info
        if ( defined $active_parameter_href->{rank_model_file} ) {

            my $rank_model_version;
            if ( $active_parameter_href->{rank_model_file} =~
                / v(\d+[.]\d+[.]\d+ | \d+[.]\d+) /sxm )
            {

                $rank_model_version = $1;
            }
            add_program_metafile_to_sample_info(
                {
                    file =>
                      basename( $active_parameter_href->{rank_model_file} ),
                    metafile_tag => q{rank_model},
                    path         => $active_parameter_href->{rank_model_file},
                    program_name => q{genmod},
                    sample_info_href => $sample_info_href,
                    version          => $rank_model_version,
                }
            );
        }
    }

    ## Track the number of created xargs scripts per module for Block algorithm
    return $xargs_file_counter;
}

sub analysis_rankvariant_rio_unaffected {

## Function : Annotate and score variants depending on mendelian inheritance, frequency and phenotype etc.
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
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $FILEHANDLE;
    my $file_info_href;
    my $file_path;
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
        file_path               => { strict_type => 1, store => \$file_path },
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

    use MIP::Check::Cluster qw{ check_max_core_number };
    use MIP::Cluster qw{ get_core_number };
    use MIP::Delete::List qw{ delete_contig_elements };
    use MIP::File::Format::Pedigree qw{ create_fam_file };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Genmod
      qw{ genmod_annotate genmod_compound genmod_models genmod_score };
    use MIP::QC::Record
      qw{ add_program_metafile_to_sample_info add_program_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Set::File qw{ set_file_suffix };

    ## Constant
    Readonly my $CORE_NUMBER_REQUESTED => 16;
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
    my $job_id_chain            = $parameter_href->{$mip_program_name}{chain};
    my $reduce_io_ref           = \$active_parameter_href->{reduce_io};
    my $vcfparser_analysis_type = $EMPTY_STR;
    my $xargs_file_path_prefix;
    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    ## Set default contigs
    my @contigs_size_ordered = @{ $file_info_href->{contigs_size_ordered} };
    my @contigs              = @{ $file_info_href->{contigs} };

    ## Filehandles
    # Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            module_core_number   => $core_number,
            modifier_core_number => scalar(@contigs),
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    ### Detect the number of cores to use per genmod process.
    ## Limit number of cores requested to the maximum number of cores available per node
    my $genmod_core_number = check_max_core_number(
        {
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
            core_number_requested => $CORE_NUMBER_REQUESTED,
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
    my $infile_tag = $file_info_href->{$family_id}{psnpeff}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$mip_program_name}{file_tag};

    ## Files
    my $infile_prefix  = $family_id . $infile_tag . $call_type;
    my $outfile_prefix = $family_id . $outfile_tag . $call_type;

    ## Paths
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ## Assign suffix
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

    my $family_file =
      catfile( $outfamily_file_directory, $family_id . $DOT . q{fam} );

    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            active_parameter_href => $active_parameter_href,
            fam_file_path         => $family_file,
            FILEHANDLE            => $FILEHANDLE,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Determined by vcfparser output
    for my $vcfparser_outfile_counter ( 0 .. $VCFPARSER_OUTFILE_COUNT ) {

        if ( $vcfparser_outfile_counter == 1 ) {

            ## SelectFile variants
            $vcfparser_analysis_type = $DOT . q{selected};
            ## Selectfile contigs
            @contigs_size_ordered =
              @{ $file_info_href->{sorted_select_file_contigs} };
            @contigs = @{ $file_info_href->{select_file_contigs} };

            ## Remove MT|M since no exome kit so far has mitochondrial probes
            if ( $consensus_analysis_type eq q{wes} ) {

                ## Removes an element from array and return new array while leaving orginal elements_ref untouched
                @contigs = delete_contig_elements(
                    {
                        elements_ref =>
                          \@{ $file_info_href->{select_file_contigs} },
                        remove_contigs_ref => [q{ MT M }],
                    }
                );

                ## Removes an element from array and return new array while leaving orginal elements_ref untouched
                @contigs_size_ordered = delete_contig_elements(
                    {
                        elements_ref =>
                          \@{ $file_info_href->{sorted_select_file_contigs} },
                        remove_contigs_ref => [qw{ MT M }],
                    }
                );
            }
        }

        ## Genmod
        say {$FILEHANDLE} q{## Genmod};

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                core_number        => $genmod_core_number,
                FILEHANDLE         => $FILEHANDLE,
                file_path          => $file_path,
                program_info_path  => $program_info_path,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                xargs_file_counter => $xargs_file_counter,
            }
        );

        ## Track which genmod modules has been processed
        my $genmod_module = $EMPTY_STR;

        ## Process per contig
        while ( my ( $contig_index, $contig ) = each @contigs_size_ordered ) {

            ## Get parameters
            # Restart for next contig
            $genmod_module = $EMPTY_STR;

            ## InFile
            my $genmod_indata =
                $file_path_prefix
              . $UNDERSCORE
              . $contig
              . $vcfparser_analysis_type
              . $infile_suffix;

            ## Output file
            my $genmod_outfile_path =
                $outfile_path_prefix
              . $UNDERSCORE
              . $contig
              . $vcfparser_analysis_type
              . $infile_suffix;

            ## Genmod Annotate
            $genmod_module = q{_annotate};

            genmod_annotate(
                {
                    annotate_region =>
                      $active_parameter_href->{genmod_annotate_regions},
                    cadd_file_paths_ref =>
                      \@{ $active_parameter_href->{genmod_annotate_cadd_files}
                      },
                    FILEHANDLE   => $XARGSFILEHANDLE,
                    infile_path  => $genmod_indata,
                    outfile_path => $genmod_outfile_path,
                    spidex_file_path =>
                      $active_parameter_href->{genmod_annotate_spidex_file},
                    stderrfile_path => $xargs_file_path_prefix
                      . $DOT
                      . $contig
                      . $genmod_module
                      . $DOT
                      . q{stderr.txt},
                    temp_directory_path => $temp_directory,
                    verbosity           => q{v},
                }
            );
            say {$XARGSFILEHANDLE} $NEWLINE;
        }

        ## QC Data File(s)
        migrate_file(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $outfile_path_prefix
                  . $UNDERSCORE
                  . $file_info_href->{contigs_size_ordered}[0]
                  . $vcfparser_analysis_type
                  . $outfile_suffix,
                outfile_path => $outfamily_directory,
            }
        );
        say {$FILEHANDLE} q{wait}, $NEWLINE;
    }

    if ( $mip_program_mode == 1 ) {

        my $qc_genmod_outfile =
            $outfile_prefix
          . $UNDERSCORE
          . $file_info_href->{contigs_size_ordered}[0]
          . $vcfparser_analysis_type
          . $outfile_suffix;
        add_program_outfile_to_sample_info(
            {
                outdirectory     => $outfamily_directory,
                outfile          => $qc_genmod_outfile,
                program_name     => q{genmod},
                sample_info_href => $sample_info_href,
            }
        );

        # Add to Sample_info
        if ( defined $active_parameter_href->{rank_model_file} ) {

            my $rank_model_version;
            if ( $active_parameter_href->{rank_model_file} =~
                / v(\d+[.]\d+[.]\d+ | \d+[.]\d+) /sxm )
            {

                $rank_model_version = $1;
            }
            add_program_metafile_to_sample_info(
                {
                    file =>
                      basename( $active_parameter_href->{rank_model_file} ),
                    metafile_tag => q{rank_model},
                    path         => $active_parameter_href->{rank_model_file},
                    program_name => q{genmod},
                    sample_info_href => $sample_info_href,
                    version          => $rank_model_version,
                }
            );
        }
    }

    ## Track the number of created xargs scripts per module for Block algorithm
    return $xargs_file_counter;
}

sub analysis_rankvariant_unaffected {

## Function : Annotate variants but do not score.
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
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $file_path;
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
        file_path               => { strict_type => 1, store => \$file_path },
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

    use MIP::Check::Cluster qw{ check_max_core_number };
    use MIP::Cluster qw{ get_core_number };
    use MIP::Delete::List qw{ delete_contig_elements };
    use MIP::File::Format::Pedigree qw{ create_fam_file };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Genmod
      qw{ genmod_annotate genmod_compound genmod_models genmod_score };
    use MIP::QC::Record
      qw{ add_program_metafile_to_sample_info add_program_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Set::File qw{ set_file_suffix };
    use MIP::Script::Setup_script qw(setup_script);

    ## Constant
    Readonly my $CORE_NUMBER_REQUESTED => 16;
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
    my $job_id_chain            = $parameter_href->{$mip_program_name}{chain};
    my $reduce_io_ref           = \$active_parameter_href->{reduce_io};
    my $vcfparser_analysis_type = $EMPTY_STR;
    my $xargs_file_path_prefix;
    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    ## Set default contigs
    my @contigs_size_ordered = @{ $file_info_href->{contigs_size_ordered} };
    my @contigs              = @{ $file_info_href->{contigs} };

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            module_core_number   => $core_number,
            modifier_core_number => scalar(@contigs),
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    ### Detect the number of cores to use per genmod process.
    ## Limit number of cores requested to the maximum number of cores available per node
    my $genmod_core_number = check_max_core_number(
        {
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
            core_number_requested => $CORE_NUMBER_REQUESTED,
        }
    );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $family_id,
            FILEHANDLE                      => $FILEHANDLE,
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
    my $infile_tag = $file_info_href->{$family_id}{psnpeff}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$mip_program_name}{file_tag};

    ## Files
    my $infile_prefix  = $family_id . $infile_tag . $call_type;
    my $outfile_prefix = $family_id . $outfile_tag . $call_type;

    ## Paths
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ## Assign suffix
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

    my $family_file =
      catfile( $outfamily_file_directory, $family_id . $DOT . q{fam} );

    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            active_parameter_href => $active_parameter_href,
            fam_file_path         => $family_file,
            FILEHANDLE            => $FILEHANDLE,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Determined by vcfparser output
    for my $vcfparser_outfile_counter ( 0 .. $VCFPARSER_OUTFILE_COUNT ) {

        if ( $vcfparser_outfile_counter == 1 ) {

            ## SelectFile variants
            $vcfparser_analysis_type = $DOT . q{selected};
            ## Selectfile contigs
            @contigs_size_ordered =
              @{ $file_info_href->{sorted_select_file_contigs} };
            @contigs = @{ $file_info_href->{select_file_contigs} };

            ## Remove MT|M since no exome kit so far has mitochondrial probes
            if ( $consensus_analysis_type eq q{wes} ) {

                ## Removes an element from array and return new array while leaving orginal elements_ref untouched
                @contigs = delete_contig_elements(
                    {
                        elements_ref =>
                          \@{ $file_info_href->{select_file_contigs} },
                        remove_contigs_ref => [q{ MT M }],
                    }
                );

                ## Removes an element from array and return new array while leaving orginal elements_ref untouched
                @contigs_size_ordered = delete_contig_elements(
                    {
                        elements_ref =>
                          \@{ $file_info_href->{sorted_select_file_contigs} },
                        remove_contigs_ref => [qw{ MT M }],
                    }
                );
            }
        }

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        $xargs_file_counter = xargs_migrate_contig_files(
            {
                contigs_ref => \@contigs_size_ordered,
                core_number => $core_number,
                FILEHANDLE  => $FILEHANDLE,
                file_ending => $vcfparser_analysis_type
                  . $infile_suffix
                  . $ASTERIX,
                file_path          => $file_path,
                indirectory        => $infamily_directory,
                infile             => $infile_prefix,
                program_info_path  => $program_info_path,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                temp_directory     => $active_parameter_href->{temp_directory},
                xargs_file_counter => $xargs_file_counter,
            }
        );

        ## Genmod
        say {$FILEHANDLE} q{## Genmod};

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                core_number        => $genmod_core_number,
                FILEHANDLE         => $FILEHANDLE,
                file_path          => $file_path,
                program_info_path  => $program_info_path,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                xargs_file_counter => $xargs_file_counter,
            }
        );

        ## Track which genmod modules has been processed
        my $genmod_module = $EMPTY_STR;

        ## Process per contig
        while ( my ( $contig_index, $contig ) = each @contigs_size_ordered ) {

            ## Get parameters
            # Restart for next contig
            $genmod_module = $EMPTY_STR;

            ## Infile
            my $genmod_indata =
                $file_path_prefix
              . $UNDERSCORE
              . $contig
              . $vcfparser_analysis_type
              . $infile_suffix;

            ## Output file
            my $genmod_outfile_path =
                $outfile_path_prefix
              . $UNDERSCORE
              . $contig
              . $vcfparser_analysis_type
              . $infile_suffix;

            ## Genmod Annotate
            $genmod_module = q{_annotate};

            genmod_annotate(
                {
                    annotate_region =>
                      $active_parameter_href->{genmod_annotate_regions},
                    cadd_file_paths_ref =>
                      \@{ $active_parameter_href->{genmod_annotate_cadd_files}
                      },
                    FILEHANDLE   => $XARGSFILEHANDLE,
                    infile_path  => $genmod_indata,
                    outfile_path => $genmod_outfile_path,
                    spidex_file_path =>
                      $active_parameter_href->{genmod_annotate_spidex_file},
                    stderrfile_path => $xargs_file_path_prefix
                      . $DOT
                      . $contig
                      . $genmod_module
                      . $DOT
                      . q{stderr.txt},
                    temp_directory_path => $temp_directory,
                    verbosity           => q{v},
                }
            );
            say {$XARGSFILEHANDLE} $NEWLINE;
        }

        ## Copies file from temporary directory. Per contig
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        ($xargs_file_counter) = xargs_migrate_contig_files(
            {
                contigs_ref => \@contigs_size_ordered,
                core_number => $active_parameter_href->{max_cores_per_node},
                FILEHANDLE  => $FILEHANDLE,
                file_path   => $file_path,
                file_ending => $vcfparser_analysis_type
                  . $outfile_suffix
                  . $ASTERIX,
                outdirectory       => $outfamily_directory,
                outfile            => $outfile_prefix,
                program_info_path  => $program_info_path,
                temp_directory     => $temp_directory,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                xargs_file_counter => $xargs_file_counter,
            }
        );
    }

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $mip_program_mode == 1 ) {

        my $qc_genmod_outfile =
            $outfile_prefix
          . $UNDERSCORE
          . $file_info_href->{contigs_size_ordered}[0]
          . $vcfparser_analysis_type
          . $outfile_suffix;
        add_program_outfile_to_sample_info(
            {
                outdirectory     => $outfamily_directory,
                outfile          => $qc_genmod_outfile,
                program_name     => q{genmod},
                sample_info_href => $sample_info_href,
            }
        );

        # Add to Sample_info
        if ( defined $active_parameter_href->{rank_model_file} ) {

            my $rank_model_version;
            if ( $active_parameter_href->{rank_model_file} =~
                / v(\d+[.]\d+[.]\d+ | \d+[.]\d+) /sxm )
            {

                $rank_model_version = $1;
            }
            add_program_metafile_to_sample_info(
                {
                    file =>
                      basename( $active_parameter_href->{rank_model_file} ),
                    metafile_tag => q{rank_model},
                    path         => $active_parameter_href->{rank_model_file},
                    program_name => q{genmod},
                    sample_info_href => $sample_info_href,
                    version          => $rank_model_version,
                }
            );
        }
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
