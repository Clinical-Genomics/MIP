package MIP::Recipes::Analysis::Snpeff;

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
    our @EXPORT_OK = qw{ analysis_snpeff analysis_snpeff_rio };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $DOT        => q{.};
Readonly my $EMPTY_STR  => q{};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_snpeff {

## Function : snpeff annotates variants from different sources.
## Returns  : |$xargs_file_counter
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => Variant call type
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF
##          : $file_path               => File path
##          : $infamily_directory      => In family directory
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outfamily_directory     => Out family directory
##          : $parameter_href          => Parameter hash {REF}
##          : $program_info_path       => The program info path
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $file_path;
    my $infamily_directory;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $outfamily_directory;
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
        program_info_path =>
          { strict_type => 1, store => \$program_info_path, },
        program_name => {
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

    use MIP::Cluster qw{ get_core_number };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Mip_vcfparser qw{ mip_vcfparser };
    use MIP::Program::Variantcalling::Snpeff qw{ snpeff_ann };
    use MIP::Program::Variantcalling::Snpsift
      qw{ snpsift_annotate snpsift_dbnsfp };
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
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};

    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
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
            process_time                    => $time,
            program_directory               => $outaligner_dir,
            program_name                    => $program_name,
            source_environment_commands_ref => [$source_environment_cmd],
            temp_directory                  => $temp_directory,
        }
    );

    # Used downstream
    $parameter_href->{$mip_program_name}{indirectory} = $outfamily_directory;

    ## Tags
    my $infile_tag = $file_info_href->{$family_id}{pvcfparser}{file_tag};
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

    my $vcfparser_analysis_type = $EMPTY_STR;

    # Set default
    my $vcfparser_contigs_ref = \@{ $file_info_href->{contigs_size_ordered} };

    for (
        my $vcfparser_outfile_counter = 0 ;
        $vcfparser_outfile_counter <
        $active_parameter_href->{vcfparser_outfile_count} ;
        $vcfparser_outfile_counter++
      )
    {

        if ( $vcfparser_outfile_counter == 1 ) {

            # SelectFile variants
            $vcfparser_analysis_type = $DOT . q{selected};

            # Selectfile contigs
            $vcfparser_contigs_ref =
              \@{ $file_info_href->{sorted_select_file_contigs} };
        }

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        ($xargs_file_counter) = xargs_migrate_contig_files(
            {
                contigs_ref => $vcfparser_contigs_ref,
                core_number => $core_number,
                FILEHANDLE  => $FILEHANDLE,
                file_ending => $vcfparser_analysis_type
                  . $infile_suffix
                  . $ASTERISK,
                file_path          => $file_path,
                indirectory        => $infamily_directory,
                infile             => $infile_prefix,
                program_info_path  => $program_info_path,
                temp_directory     => $temp_directory,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                xargs_file_counter => $xargs_file_counter,
            }
        );

        ## SnpSift Annotation
        say {$FILEHANDLE} q{## SnpSift Annotation};

        my $annotation_file_counter = 0;
        my $xargs_file_path_prefix;

        # Annotate using snpeff
        if ( $active_parameter_href->{snpeff_ann} == 1 ) {
            ## Create file commands for xargs
            ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
                {
                    core_number   => $core_number,
                    first_command => q{java},
                    FILEHANDLE    => $FILEHANDLE,
                    file_path     => $file_path,
                    java_jar      => catfile(
                        $active_parameter_href->{snpeff_path},
                        q{snpEff.jar}
                    ),
                    java_use_large_pages =>
                      $active_parameter_href->{java_use_large_pages},
                    memory_allocation  => q{Xmx4g -XX:-UseConcMarkSweepGC},
                    program_info_path  => $program_info_path,
                    temp_directory     => $temp_directory,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    xargs_file_counter => $xargs_file_counter,
                }
            );

            foreach my $contig ( @{$vcfparser_contigs_ref} ) {

                snpeff_ann(
                    {
                        config_file_path => catfile(
                            $active_parameter_href->{snpeff_path},
                            q{snpEff.config}
                        ),
                        FILEHANDLE => $XARGSFILEHANDLE,
                        genome_build_version =>
                          $active_parameter_href->{snpeff_genome_build_version},
                        infile_path => $file_path_prefix
                          . $UNDERSCORE
                          . $contig
                          . $vcfparser_analysis_type
                          . $infile_suffix,
                        stderrfile_path => $xargs_file_path_prefix
                          . $UNDERSCORE
                          . $contig
                          . $DOT
                          . q{stderr.txt},
                        stdoutfile_path => $file_path_prefix
                          . $UNDERSCORE
                          . $contig
                          . $vcfparser_analysis_type
                          . $infile_suffix
                          . $DOT
                          . $xargs_file_counter,
                        verbosity => q{v},

                    }
                );
                say {$XARGSFILEHANDLE} $NEWLINE;
            }
            $annotation_file_counter = $xargs_file_counter;
        }

        while ( my ( $annotation_file, $annotation_info_key ) =
            each %{ $active_parameter_href->{snpsift_annotation_files} } )
        {

            ## Create file commands for xargs
            ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
                {
                    core_number   => $core_number,
                    FILEHANDLE    => $FILEHANDLE,
                    file_path     => $file_path,
                    first_command => q{java},
                    java_jar      => catfile(
                        $active_parameter_href->{snpeff_path},
                        q{SnpSift.jar}
                    ),
                    java_use_large_pages =>
                      $active_parameter_href->{java_use_large_pages},
                    memory_allocation  => q{Xmx2g -XX:-UseConcMarkSweepGC},
                    program_info_path  => $program_info_path,
                    temp_directory     => $temp_directory,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    xargs_file_counter => $xargs_file_counter,
                }
            );
            ## Get parameters
            my $name_prefix;
            my $info_key;
            if ( defined $annotation_info_key ) {

                ## Apply specific INFO field output key for easier downstream processing
                if (
                    defined(
                        $active_parameter_href->{snpsift_annotation_outinfo_key}
                          {$annotation_file}
                    )
                  )
                {

                    $name_prefix =
                      $active_parameter_href->{snpsift_annotation_outinfo_key}
                      {$annotation_file};
                }
                $info_key = $annotation_info_key;    #Database
            }

            foreach my $contig ( @{$vcfparser_contigs_ref} ) {

                ## Get contig specific parameters
                my $infile_path;

                # First file per contig
                if ( !$annotation_file_counter ) {

                    $infile_path =
                        $file_path_prefix
                      . $UNDERSCORE
                      . $contig
                      . $vcfparser_analysis_type
                      . $infile_suffix;
                }
                else {

                    my $annotation_infile_number = $xargs_file_counter - 1;
                    $infile_path =
                        $file_path_prefix
                      . $UNDERSCORE
                      . $contig
                      . $vcfparser_analysis_type
                      . $infile_suffix
                      . $DOT
                      . $annotation_infile_number;
                }
                snpsift_annotate(
                    {
                        config_file_path => catfile(
                            $active_parameter_href->{snpeff_path},
                            q{snpEff.config}
                        ),
                        database_path   => $annotation_file,
                        FILEHANDLE      => $XARGSFILEHANDLE,
                        infile_path     => $infile_path,
                        info            => $info_key,
                        name_prefix     => $name_prefix,
                        stderrfile_path => $xargs_file_path_prefix
                          . $DOT
                          . $contig
                          . $DOT
                          . q{stderr.txt},
                        stderrfile_path_append => $xargs_file_path_prefix
                          . $DOT
                          . $contig
                          . $DOT
                          . q{stderr.txt},
                        stdoutfile_path => $file_path_prefix
                          . $UNDERSCORE
                          . $contig
                          . $vcfparser_analysis_type
                          . $infile_suffix
                          . $DOT
                          . $xargs_file_counter,
                        verbosity => q{v},
                    }
                );
                say {$XARGSFILEHANDLE} $NEWLINE;
            }

            # Increment counter
            $annotation_file_counter++;
            close $XARGSFILEHANDLE;
        }

        if ( @{ $active_parameter_href->{snpsift_dbnsfp_annotations} } ) {

            ## SnpSiftDbNSFP Annotation
            say {$FILEHANDLE} q{## SnpSiftDnNSFP Annotation};

            ## Create file commands for xargs
            ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
                {
                    core_number       => $core_number,
                    FILEHANDLE        => $FILEHANDLE,
                    file_path         => $file_path,
                    first_command     => q{java},
                    memory_allocation => q{Xmx2g -XX:-UseConcMarkSweepGC},
                    java_jar          => catfile(
                        $active_parameter_href->{snpeff_path},
                        q{SnpSift.jar}
                    ),
                    java_use_large_pages =>
                      $active_parameter_href->{java_use_large_pages},
                    program_info_path  => $program_info_path,
                    temp_directory     => $temp_directory,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    xargs_file_counter => $xargs_file_counter,
                }
            );

            my $annotation_infile_number = $xargs_file_counter - 1;

            foreach my $contig ( @{$vcfparser_contigs_ref} ) {

                snpsift_dbnsfp(
                    {
                        annotate_fields_ref => \@{
                            $active_parameter_href->{snpsift_dbnsfp_annotations}
                        },
                        config_file_path => catfile(
                            $active_parameter_href->{snpeff_path},
                            q{snpEff.config}
                        ),
                        database_path =>
                          $active_parameter_href->{snpsift_dbnsfp_file},
                        FILEHANDLE  => $XARGSFILEHANDLE,
                        infile_path => $file_path_prefix
                          . $UNDERSCORE
                          . $contig
                          . $vcfparser_analysis_type
                          . $infile_suffix
                          . $DOT
                          . $annotation_infile_number,
                        stderrfile_path => $xargs_file_path_prefix
                          . $DOT
                          . $contig
                          . $DOT
                          . q{stderr.txt},
                        stderrfile_path_append => $xargs_file_path_prefix
                          . $DOT
                          . $contig
                          . $DOT
                          . q{stderr.txt},
                        stdoutfile_path => $file_path_prefix
                          . $UNDERSCORE
                          . $contig
                          . $vcfparser_analysis_type
                          . $infile_suffix
                          . $DOT
                          . $xargs_file_counter,
                        verbosity => q{v},
                    }
                );
                say {$XARGSFILEHANDLE} $NEWLINE;
            }
            close $XARGSFILEHANDLE;
        }

        ## Add INFO headers and FIX_INFO for annotations using vcfparser
        say {$FILEHANDLE}
          q{## Add INFO headers and FIX_INFO for annotations using vcfparser};

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

        my $annotation_infile_number = $xargs_file_counter - 1;

        foreach my $contig ( @{$vcfparser_contigs_ref} ) {

            mip_vcfparser(
                {
                    FILEHANDLE  => $XARGSFILEHANDLE,
                    infile_path => $file_path_prefix
                      . $UNDERSCORE
                      . $contig
                      . $vcfparser_analysis_type
                      . $infile_suffix
                      . $DOT
                      . $annotation_infile_number,
                    stderrfile_path => $xargs_file_path_prefix
                      . $DOT
                      . $contig
                      . $DOT
                      . q{stderr.txt}
                      . $SPACE,
                    stderrfile_path_append => $xargs_file_path_prefix
                      . $DOT
                      . $contig
                      . $DOT
                      . q{stderr.txt}
                      . $SPACE,
                    stdoutfile_path => $outfile_path_prefix
                      . $UNDERSCORE
                      . $contig
                      . $vcfparser_analysis_type
                      . $outfile_suffix,
                }
            );
            say {$XARGSFILEHANDLE} $NEWLINE;
        }

        ## Copies file from temporary directory. Per contig
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        ($xargs_file_counter) = xargs_migrate_contig_files(
            {
                contigs_ref => $vcfparser_contigs_ref,
                core_number => $active_parameter_href->{max_cores_per_node},
                FILEHANDLE  => $FILEHANDLE,
                file_ending => $vcfparser_analysis_type
                  . $outfile_suffix
                  . $ASTERISK,
                file_path          => $file_path,
                outfile            => $outfile_prefix,
                outdirectory       => $outfamily_directory,
                program_info_path  => $program_info_path,
                temp_directory     => $temp_directory,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                xargs_file_counter => $xargs_file_counter,
            }
        );
    }

    if ( $mip_program_mode == 1 ) {

        ## Collect QC metadata info for later use
        my $qc_snpeff_outfile =
            $outfile_prefix
          . $UNDERSCORE
          . $file_info_href->{contigs_size_ordered}[0]
          . $vcfparser_analysis_type
          . $outfile_suffix;
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => $program_name,
                outdirectory     => $outfamily_directory,
                outfile          => $qc_snpeff_outfile,
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

sub analysis_snpeff_rio {

## Function : snpeff annotates variants from different sources.
## Returns  : |$xargs_file_counter
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => Variant call type
##          : $family_id               => Family id
##          : $FILEHANDLE              => Filehandle to write to
##          : $file_info_href          => File_info hash {REF
##          : $file_path               => File path
##          : $infamily_directory      => In family directory
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outfamily_directory     => Out family directory
##          : $parameter_href          => Parameter hash {REF}
##          : $program_info_path       => The program info path
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $FILEHANDLE;
    my $file_info_href;
    my $file_path;
    my $infamily_directory;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $outfamily_directory;
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
        program_info_path =>
          { strict_type => 1, store => \$program_info_path, },
        program_name => {
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

    use MIP::Cluster qw{ get_core_number };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Mip_vcfparser qw{ mip_vcfparser };
    use MIP::Program::Variantcalling::Snpeff qw{ snpeff_ann };
    use MIP::Program::Variantcalling::Snpsift
      qw{ snpsift_annotate snpsift_dbnsfp };
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
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};

    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
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

    # Used downstream
    $parameter_href->{$mip_program_name}{indirectory} = $outfamily_directory;

    ## Tags
    my $infile_tag = $file_info_href->{$family_id}{pvcfparser}{file_tag};
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

    my $vcfparser_analysis_type = $EMPTY_STR;

    # Set default
    my $vcfparser_contigs_ref = \@{ $file_info_href->{contigs_size_ordered} };

    for (
        my $vcfparser_outfile_counter = 0 ;
        $vcfparser_outfile_counter <
        $active_parameter_href->{vcfparser_outfile_count} ;
        $vcfparser_outfile_counter++
      )
    {

        if ( $vcfparser_outfile_counter == 1 ) {

            # SelectFile variants
            $vcfparser_analysis_type = $DOT . q{selected};

            # Selectfile contigs
            $vcfparser_contigs_ref =
              \@{ $file_info_href->{sorted_select_file_contigs} };
        }

        ## SnpSift Annotation
        say {$FILEHANDLE} q{## SnpSift Annotation};

        my $annotation_file_counter = 0;
        my $xargs_file_path_prefix;

        # Annotate using snpeff
        if ( $active_parameter_href->{snpeff_ann} == 1 ) {
            ## Create file commands for xargs
            ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
                {
                    core_number   => $core_number,
                    first_command => q{java},
                    FILEHANDLE    => $FILEHANDLE,
                    file_path     => $file_path,
                    java_jar      => catfile(
                        $active_parameter_href->{snpeff_path},
                        q{snpEff.jar}
                    ),
                    java_use_large_pages =>
                      $active_parameter_href->{java_use_large_pages},
                    memory_allocation  => q{Xmx4g -XX:-UseConcMarkSweepGC},
                    program_info_path  => $program_info_path,
                    temp_directory     => $temp_directory,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    xargs_file_counter => $xargs_file_counter,
                }
            );

            foreach my $contig ( @{$vcfparser_contigs_ref} ) {

                snpeff_ann(
                    {
                        config_file_path => catfile(
                            $active_parameter_href->{snpeff_path},
                            q{snpEff.config}
                        ),
                        FILEHANDLE => $XARGSFILEHANDLE,
                        genome_build_version =>
                          $active_parameter_href->{snpeff_genome_build_version},
                        infile_path => $file_path_prefix
                          . $UNDERSCORE
                          . $contig
                          . $vcfparser_analysis_type
                          . $infile_suffix,
                        stderrfile_path => $xargs_file_path_prefix
                          . $UNDERSCORE
                          . $contig
                          . $DOT
                          . q{stderr.txt},
                        stdoutfile_path => $file_path_prefix
                          . $UNDERSCORE
                          . $contig
                          . $vcfparser_analysis_type
                          . $infile_suffix
                          . $DOT
                          . $xargs_file_counter,
                        verbosity => q{v},

                    }
                );
                say {$XARGSFILEHANDLE} $NEWLINE;
            }
            $annotation_file_counter = $xargs_file_counter;
        }

        while ( my ( $annotation_file, $annotation_info_key ) =
            each %{ $active_parameter_href->{snpsift_annotation_files} } )
        {

            ## Create file commands for xargs
            ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
                {
                    core_number   => $core_number,
                    FILEHANDLE    => $FILEHANDLE,
                    file_path     => $file_path,
                    first_command => q{java},
                    java_jar      => catfile(
                        $active_parameter_href->{snpeff_path},
                        q{SnpSift.jar}
                    ),
                    java_use_large_pages =>
                      $active_parameter_href->{java_use_large_pages},
                    memory_allocation  => q{Xmx2g -XX:-UseConcMarkSweepGC},
                    program_info_path  => $program_info_path,
                    temp_directory     => $temp_directory,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    xargs_file_counter => $xargs_file_counter,
                }
            );
            ## Get parameters
            my $name_prefix;
            my $info_key;
            if ( defined $annotation_info_key ) {

                ## Apply specific INFO field output key for easier downstream processing
                if (
                    defined(
                        $active_parameter_href->{snpsift_annotation_outinfo_key}
                          {$annotation_file}
                    )
                  )
                {

                    $name_prefix =
                      $active_parameter_href->{snpsift_annotation_outinfo_key}
                      {$annotation_file};
                }
                $info_key = $annotation_info_key;    #Database
            }

            foreach my $contig ( @{$vcfparser_contigs_ref} ) {

                ## Get contig specific parameters
                my $infile_path;

                # First file per contig
                if ( !$annotation_file_counter ) {

                    $infile_path =
                        $file_path_prefix
                      . $UNDERSCORE
                      . $contig
                      . $vcfparser_analysis_type
                      . $infile_suffix;
                }
                else {

                    my $annotation_infile_number = $xargs_file_counter - 1;
                    $infile_path =
                        $file_path_prefix
                      . $UNDERSCORE
                      . $contig
                      . $vcfparser_analysis_type
                      . $infile_suffix
                      . $DOT
                      . $annotation_infile_number;
                }
                snpsift_annotate(
                    {
                        config_file_path => catfile(
                            $active_parameter_href->{snpeff_path},
                            q{snpEff.config}
                        ),
                        database_path   => $annotation_file,
                        FILEHANDLE      => $XARGSFILEHANDLE,
                        infile_path     => $infile_path,
                        info            => $info_key,
                        name_prefix     => $name_prefix,
                        stderrfile_path => $xargs_file_path_prefix
                          . $DOT
                          . $contig
                          . $DOT
                          . q{stderr.txt},
                        stderrfile_path_append => $xargs_file_path_prefix
                          . $DOT
                          . $contig
                          . $DOT
                          . q{stderr.txt},
                        stdoutfile_path => $file_path_prefix
                          . $UNDERSCORE
                          . $contig
                          . $vcfparser_analysis_type
                          . $infile_suffix
                          . $DOT
                          . $xargs_file_counter,
                        verbosity => q{v},
                    }
                );
                say {$XARGSFILEHANDLE} $NEWLINE;
            }

            # Increment counter
            $annotation_file_counter++;
            close $XARGSFILEHANDLE;
        }

        if ( @{ $active_parameter_href->{snpsift_dbnsfp_annotations} } ) {

            ## SnpSiftDbNSFP Annotation
            say {$FILEHANDLE} q{## SnpSiftDnNSFP Annotation};

            ## Create file commands for xargs
            ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
                {
                    core_number       => $core_number,
                    FILEHANDLE        => $FILEHANDLE,
                    file_path         => $file_path,
                    first_command     => q{java},
                    memory_allocation => q{Xmx2g -XX:-UseConcMarkSweepGC},
                    java_jar          => catfile(
                        $active_parameter_href->{snpeff_path},
                        q{SnpSift.jar}
                    ),
                    java_use_large_pages =>
                      $active_parameter_href->{java_use_large_pages},
                    program_info_path  => $program_info_path,
                    temp_directory     => $temp_directory,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    xargs_file_counter => $xargs_file_counter,
                }
            );

            my $annotation_infile_number = $xargs_file_counter - 1;

            foreach my $contig ( @{$vcfparser_contigs_ref} ) {

                snpsift_dbnsfp(
                    {
                        annotate_fields_ref => \@{
                            $active_parameter_href->{snpsift_dbnsfp_annotations}
                        },
                        config_file_path => catfile(
                            $active_parameter_href->{snpeff_path},
                            q{snpEff.config}
                        ),
                        database_path =>
                          $active_parameter_href->{snpsift_dbnsfp_file},
                        FILEHANDLE  => $XARGSFILEHANDLE,
                        infile_path => $file_path_prefix
                          . $UNDERSCORE
                          . $contig
                          . $vcfparser_analysis_type
                          . $infile_suffix
                          . $DOT
                          . $annotation_infile_number,
                        stderrfile_path => $xargs_file_path_prefix
                          . $DOT
                          . $contig
                          . $DOT
                          . q{stderr.txt},
                        stderrfile_path_append => $xargs_file_path_prefix
                          . $DOT
                          . $contig
                          . $DOT
                          . q{stderr.txt},
                        stdoutfile_path => $file_path_prefix
                          . $UNDERSCORE
                          . $contig
                          . $vcfparser_analysis_type
                          . $infile_suffix
                          . $DOT
                          . $xargs_file_counter,
                        verbosity => q{v},
                    }
                );
                say {$XARGSFILEHANDLE} $NEWLINE;
            }
            close $XARGSFILEHANDLE;
        }

        ## Add INFO headers and FIX_INFO for annotations using vcfparser
        say {$FILEHANDLE}
          q{## Add INFO headers and FIX_INFO for annotations using vcfparser};

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

        my $annotation_infile_number = $xargs_file_counter - 1;

        foreach my $contig ( @{$vcfparser_contigs_ref} ) {

            mip_vcfparser(
                {
                    FILEHANDLE  => $XARGSFILEHANDLE,
                    infile_path => $file_path_prefix
                      . $UNDERSCORE
                      . $contig
                      . $vcfparser_analysis_type
                      . $infile_suffix
                      . $DOT
                      . $annotation_infile_number,
                    stderrfile_path => $xargs_file_path_prefix
                      . $DOT
                      . $contig
                      . $DOT
                      . q{stderr.txt}
                      . $SPACE,
                    stderrfile_path_append => $xargs_file_path_prefix
                      . $DOT
                      . $contig
                      . $DOT
                      . q{stderr.txt}
                      . $SPACE,
                    stdoutfile_path => $outfile_path_prefix
                      . $UNDERSCORE
                      . $contig
                      . $vcfparser_analysis_type
                      . $outfile_suffix,
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

        ## Collect QC metadata info for later use
        my $qc_snpeff_outfile =
            $outfile_prefix
          . $UNDERSCORE
          . $file_info_href->{contigs_size_ordered}[0]
          . $vcfparser_analysis_type
          . $outfile_suffix;
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => $program_name,
                outdirectory     => $outfamily_directory,
                outfile          => $qc_snpeff_outfile,
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );
    }

    # Track the number of created xargs scripts per module for Block algorithm
    return $xargs_file_counter;
}
1;
