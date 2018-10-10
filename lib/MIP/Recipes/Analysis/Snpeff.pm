package MIP::Recipes::Analysis::Snpeff;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename };
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
    our $VERSION = 1.02;

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

## Function : Snpeff annotates variants from different sources.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
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
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_info_path;
    my $program_name;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $temp_directory;
    my $xargs_file_counter;

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
        file_path               => { store => \$file_path, strict_type => 1, },
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
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        program_info_path =>
          { store => \$program_info_path, strict_type => 1, },
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
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
        xargs_file_counter => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw{ get_core_number };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_module_parameters get_program_attributes };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Mip_vcfparser qw{ mip_vcfparser };
    use MIP::Program::Variantcalling::Snpeff qw{ snpeff_ann };
    use MIP::Program::Variantcalling::Snpsift
      qw{ snpsift_annotate snpsift_dbnsfp };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $family_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            program_name   => $program_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my %infile_path        = %{ $io{in}{file_path_href} };
    my $job_id_chain       = get_program_attributes(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
            attribute      => q{chain},
        }
    );
    my $program_mode = $active_parameter_href->{$program_name};
    my %snpsift_annotation_outinfo_key =
      %{ $active_parameter_href->{snpsift_annotation_outinfo_key} };
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $family_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => [ ( keys %infile_path ) ],
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                program_name     => $program_name,
                temp_directory   => $temp_directory,
            }
        )
    );

    my %outfile_path  = %{ $io{out}{file_path_href} };
    my @outfile_paths = @{ $io{out}{file_paths} };

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
            modifier_core_number => scalar @{ $file_info_href->{contigs} },
            module_core_number   => $core_number,
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
            log                             => $log,
            process_time                    => $time,
            program_directory               => $program_name,
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    ## SnpSift Ann Annotation
    say {$FILEHANDLE} q{## SnpSift Ann Annotation};

    my $annotation_file_counter = 0;
    my $xargs_file_path_prefix;
    my $snpeff_config_file_path =
      catfile( $active_parameter_href->{snpeff_path}, q{snpEff.config} );
    my $java_snpeff_jar =
      catfile( $active_parameter_href->{snpeff_path}, q{snpEff.jar} );

    # Annotate using snpeff
    if ( $active_parameter_href->{snpeff_ann} ) {

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                core_number   => $core_number,
                first_command => q{java},
                FILEHANDLE    => $FILEHANDLE,
                file_path     => $file_path,
                java_jar      => $java_snpeff_jar,
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                memory_allocation  => q{Xmx4g -XX:-UseConcMarkSweepGC},
                program_info_path  => $program_info_path,
                temp_directory     => $temp_directory,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                xargs_file_counter => $xargs_file_counter,
            }
        );

      CONTIG:
        foreach my $contig ( keys %infile_path ) {

            my $ann_outfile_path =
              $infile_path{$contig} . $DOT . $xargs_file_counter;
            snpeff_ann(
                {
                    config_file_path => $snpeff_config_file_path,
                    FILEHANDLE       => $XARGSFILEHANDLE,
                    genome_build_version =>
                      $active_parameter_href->{snpeff_genome_build_version},
                    infile_path     => $infile_path{$contig},
                    stderrfile_path => $xargs_file_path_prefix
                      . $UNDERSCORE
                      . $contig
                      . $DOT
                      . q{stderr.txt},
                    stdoutfile_path => $ann_outfile_path,
                    verbosity       => q{v},
                }
            );
            say {$XARGSFILEHANDLE} $NEWLINE;
        }
        $annotation_file_counter = $xargs_file_counter;
    }

  ANNOTATION:
    while ( my ( $annotation_file, $annotation_info_key ) =
        each %{ $active_parameter_href->{snpsift_annotation_files} } )
    {

        say {$FILEHANDLE} q{## SnpSift Annotate } . basename($annotation_file);

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
                defined $snpsift_annotation_outinfo_key{$annotation_file}

              )
            {

                $name_prefix =
                  $snpsift_annotation_outinfo_key{$annotation_file};
            }

            # Database
            $info_key = $annotation_info_key;
        }

      CONTIG:
        foreach my $contig ( keys %infile_path ) {

            ## Get contig specific parameters
            my $annotate_infile_path;
            my $annotate_outfile_path =
              $outfile_path{$contig} . $DOT . $xargs_file_counter;

            # First file per contig
            if ( not $annotation_file_counter ) {

                $annotate_infile_path = $infile_path{$contig};
            }
            else {

                ## Decrement to make last output file as input
                my $annotation_infile_number = $xargs_file_counter - 1;
                $annotate_infile_path =
                  $infile_path{$contig} . $DOT . $annotation_infile_number;
            }
            snpsift_annotate(
                {
                    config_file_path => $snpeff_config_file_path,
                    database_path    => $annotation_file,
                    FILEHANDLE       => $XARGSFILEHANDLE,
                    infile_path      => $annotate_infile_path,
                    info             => $info_key,
                    name_prefix      => $name_prefix,
                    stderrfile_path  => $xargs_file_path_prefix
                      . $DOT
                      . $contig
                      . $DOT
                      . q{stderr.txt},
                    stdoutfile_path => $annotate_outfile_path,
                    verbosity       => q{v},
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
        say {$FILEHANDLE} q{## SnpSift DnNSFP Annotation};

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

        ## Decrement to make last output file as input to dbnsfp
        my $annotation_infile_number = $xargs_file_counter - 1;

      CONTIG:
        foreach my $contig ( keys %infile_path ) {

            my $dbnsfp_outfile_path =
              $outfile_path{$contig} . $DOT . $xargs_file_counter;
            snpsift_dbnsfp(
                {
                    annotate_fields_ref =>
                      \@{ $active_parameter_href->{snpsift_dbnsfp_annotations}
                      },
                    config_file_path => $snpeff_config_file_path,
                    database_path =>
                      $active_parameter_href->{snpsift_dbnsfp_file},
                    FILEHANDLE  => $XARGSFILEHANDLE,
                    infile_path => $outfile_path{$contig}
                      . $DOT
                      . $annotation_infile_number,
                    stderrfile_path => $xargs_file_path_prefix
                      . $DOT
                      . $contig
                      . $DOT
                      . q{stderr.txt},
                    stdoutfile_path => $dbnsfp_outfile_path,
                    verbosity       => q{v},
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

    ## Decrement to make last output file as input to vcfparser
    my $annotation_infile_number = $xargs_file_counter - 1;

  CONTIG:
    foreach my $contig ( keys %infile_path ) {

        mip_vcfparser(
            {
                FILEHANDLE  => $XARGSFILEHANDLE,
                infile_path => $outfile_path{$contig}
                  . $DOT
                  . $annotation_infile_number,
                stderrfile_path => $xargs_file_path_prefix
                  . $DOT
                  . $contig
                  . $DOT
                  . q{stderr.txt},
                stdoutfile_path => $outfile_path{$contig}
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    if ( $program_mode == 1 ) {

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                path             => $outfile_paths[0],
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );

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

sub analysis_snpeff_rio {

## Function : Snpeff annotates variants from different sources.
## Returns  : |$xargs_file_counter
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => Variant call type
##          : $family_id               => Family id
##          : $FILEHANDLE              => Filehandle to write to
##          : $file_info_href          => File_info hash {REF
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $parameter_href          => Parameter hash {REF}
##          : $program_info_path       => The program info path
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $stderr_path             => Stderr path of the block script
##          : $temp_directory          => Temporary directory

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
    my $stderr_path;

    ## Default(s)
    my $call_type;
    my $family_id;
    my $outaligner_dir;
    my $temp_directory;
    my $xargs_file_counter;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        call_type =>
          { default => q{BOTH}, store => \$call_type, strict_type => 1, },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
            strict_type => 1,
        },
        FILEHANDLE     => { required => 1, store => \$FILEHANDLE },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_path               => { store => \$file_path, strict_type => 1, },
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
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        program_info_path =>
          { store => \$program_info_path, strict_type => 1, },
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        stderr_path => {
            defined     => 1,
            required    => 1,
            store       => \$stderr_path,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
        xargs_file_counter => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
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

    Readonly my $VCFPARSER_OUTFILE_COUNT =>
      $active_parameter_href->{vcfparser_outfile_count} - 1;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set program mode
    my $program_mode = $active_parameter_href->{$program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$program_name}{chain};

    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    my $config_file_path =
      catfile( $active_parameter_href->{snpeff_path}, q{snpEff.config} );

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

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir );
    my $outfamily_directory = $infamily_directory;

    # Used downstream
    $parameter_href->{$program_name}{indirectory} = $outfamily_directory;

    ## Tags
    my $infile_tag = $file_info_href->{$family_id}{vcfparser}{file_tag};
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
            file_suffix    => $parameter_href->{$program_name}{outfile_suffix},
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

    my $vcfparser_analysis_type = $EMPTY_STR;

    # Set default
    my $vcfparser_contigs_ref = \@{ $file_info_href->{contigs_size_ordered} };

    ## Determined by vcfparser output
  VCFPARSER_OUTFILE:
    for my $vcfparser_outfile_counter ( 0 .. $VCFPARSER_OUTFILE_COUNT ) {

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

        my $java_snpeff_jar =
          catfile( $active_parameter_href->{snpeff_path}, q{snpEff.jar} );

        # Annotate using snpeff
        if ( $active_parameter_href->{snpeff_ann} == 1 ) {
            ## Create file commands for xargs
            ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
                {
                    core_number   => $core_number,
                    first_command => q{java},
                    FILEHANDLE    => $FILEHANDLE,
                    file_path     => $file_path,
                    java_jar      => $java_snpeff_jar,
                    java_use_large_pages =>
                      $active_parameter_href->{java_use_large_pages},
                    memory_allocation  => q{Xmx4g -XX:-UseConcMarkSweepGC},
                    program_info_path  => $program_info_path,
                    temp_directory     => $temp_directory,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                    xargs_file_counter => $xargs_file_counter,
                }
            );

          CONTIG:
            foreach my $contig ( @{$vcfparser_contigs_ref} ) {

                my $snpeff_file_path_prefix =
                    $file_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $vcfparser_analysis_type;

                snpeff_ann(
                    {
                        config_file_path => $config_file_path,
                        FILEHANDLE       => $XARGSFILEHANDLE,
                        genome_build_version =>
                          $active_parameter_href->{snpeff_genome_build_version},
                        infile_path => $snpeff_file_path_prefix
                          . $infile_suffix,
                        stderrfile_path => $xargs_file_path_prefix
                          . $UNDERSCORE
                          . $contig
                          . $DOT
                          . q{stderr.txt},
                        stdoutfile_path => $snpeff_file_path_prefix
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

      ANNOTATION:
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
                    defined
                    $active_parameter_href->{snpsift_annotation_outinfo_key}
                    {$annotation_file} )
                {

                    $name_prefix =
                      $active_parameter_href->{snpsift_annotation_outinfo_key}
                      {$annotation_file};
                }

                # Database
                $info_key = $annotation_info_key;
            }

          CONTIG:
            foreach my $contig ( @{$vcfparser_contigs_ref} ) {

                ## Get contig specific parameters
                my $infile_path;

                # First file per contig
                if ( not $annotation_file_counter ) {

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
                        config_file_path => $config_file_path,
                        database_path    => $annotation_file,
                        FILEHANDLE       => $XARGSFILEHANDLE,
                        infile_path      => $infile_path,
                        info             => $info_key,
                        name_prefix      => $name_prefix,
                        stderrfile_path  => $xargs_file_path_prefix
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

          CONTIG:
            foreach my $contig ( @{$vcfparser_contigs_ref} ) {

                snpsift_dbnsfp(
                    {
                        annotate_fields_ref => \@{
                            $active_parameter_href->{snpsift_dbnsfp_annotations}
                        },
                        config_file_path => $config_file_path,
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

    if ( $program_mode == 1 ) {

        ## Collect QC metadata info for later use
        my $qc_snpeff_outfile =
            $outfile_prefix
          . $UNDERSCORE
          . $file_info_href->{contigs_size_ordered}[0]
          . $vcfparser_analysis_type
          . $outfile_suffix;
        add_program_outfile_to_sample_info(
            {
                path => catfile( $outfamily_directory, $qc_snpeff_outfile ),
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );
    }

    # Track the number of created xargs scripts per module for Block algorithm
    return $xargs_file_counter;
}
1;
