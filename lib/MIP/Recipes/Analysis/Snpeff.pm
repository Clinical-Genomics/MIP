package MIP::Recipes::Analysis::Snpeff;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use POSIX;
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ %ANALYSIS $ASTERISK $DOT $EMPTY_STR $NEWLINE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_snpeff };

}

## Constants
Readonly my $JAVA_GUEST_OS_MEMORY => $ANALYSIS{JAVA_GUEST_OS_MEMORY};

sub analysis_snpeff {

## Function : Snpeff annotates variants from different sources.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $file_path;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
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
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
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
        profile_base_command => {
            default     => q{sbatch},
            store       => \$profile_base_command,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
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

    use MIP::Cluster qw{ get_parallel_processes };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Variantcalling::Snpeff qw{ snpeff_ann };
    use MIP::Program::Variantcalling::Snpsift qw{ snpsift_annotate snpsift_dbnsfp };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger( uc q{mip_analyse} );

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $case_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my %infile_path        = %{ $io{in}{file_path_href} };
    my $job_id_chain       = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode = $active_parameter_href->{$recipe_name};
    my %snpsift_annotation_outinfo_key =
      %{ $active_parameter_href->{snpsift_annotation_outinfo_key} };
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number = $recipe_resource{core_number};

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $case_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => [ ( keys %infile_path ) ],
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
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

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $case_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            memory_allocation               => $recipe_resource{memory},
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
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
    my $java_snpeff_jar = catfile( $active_parameter_href->{snpeff_path}, q{snpEff.jar} );

    # Annotate using snpeff
    if ( $active_parameter_href->{snpeff_ann} ) {

        ## Snpeff requries more memory than snpsift
        Readonly my $JAVA_MEMORY_ALLOCATION => 4;
        my $process_memory_allocation = $JAVA_MEMORY_ALLOCATION + $JAVA_GUEST_OS_MEMORY;

        # Constrain parallelization to match available memory
        my $parallel_processes = get_parallel_processes(
            {
                process_memory_allocation => $process_memory_allocation,
                recipe_memory_allocation  => $recipe_resource{memory},
                core_number               => $core_number,
            }
        );

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                core_number          => $parallel_processes,
                first_command        => q{java},
                FILEHANDLE           => $FILEHANDLE,
                file_path            => $recipe_file_path,
                java_jar             => $java_snpeff_jar,
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                memory_allocation    => q{Xmx}
                  . $JAVA_MEMORY_ALLOCATION
                  . q{g -XX:-UseConcMarkSweepGC},
                recipe_info_path   => $recipe_info_path,
                temp_directory     => $temp_directory,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                xargs_file_counter => $xargs_file_counter,
            }
        );

      CONTIG:
        foreach my $contig ( keys %infile_path ) {

            my $ann_outfile_path = $outfile_path{$contig} . $DOT . $xargs_file_counter;
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
                file_path     => $recipe_file_path,
                first_command => q{java},
                java_jar =>
                  catfile( $active_parameter_href->{snpeff_path}, q{SnpSift.jar} ),
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                memory_allocation    => q{Xmx2g -XX:-UseConcMarkSweepGC},
                recipe_info_path     => $recipe_info_path,
                temp_directory       => $temp_directory,
                XARGSFILEHANDLE      => $XARGSFILEHANDLE,
                xargs_file_counter   => $xargs_file_counter,
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

                $name_prefix = $snpsift_annotation_outinfo_key{$annotation_file};
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
                  $outfile_path{$contig} . $DOT . $annotation_infile_number;
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

    if ( exists $active_parameter_href->{snpsift_dbnsfp_annotations}
        and @{ $active_parameter_href->{snpsift_dbnsfp_annotations} } )
    {

        ## SnpSiftDbNSFP Annotation
        say {$FILEHANDLE} q{## SnpSift DnNSFP Annotation};

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                core_number       => $core_number,
                FILEHANDLE        => $FILEHANDLE,
                file_path         => $recipe_file_path,
                first_command     => q{java},
                memory_allocation => q{Xmx2g -XX:-UseConcMarkSweepGC},
                java_jar =>
                  catfile( $active_parameter_href->{snpeff_path}, q{SnpSift.jar} ),
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                recipe_info_path     => $recipe_info_path,
                temp_directory       => $temp_directory,
                XARGSFILEHANDLE      => $XARGSFILEHANDLE,
                xargs_file_counter   => $xargs_file_counter,
            }
        );

        ## Decrement to make last output file as input to dbnsfp
        my $annotation_infile_number = $xargs_file_counter - 1;

      CONTIG:
        foreach my $contig ( keys %infile_path ) {

            my $dbnsfp_outfile_path = $outfile_path{$contig} . $DOT . $xargs_file_counter;
            snpsift_dbnsfp(
                {
                    annotate_fields_ref =>
                      \@{ $active_parameter_href->{snpsift_dbnsfp_annotations} },
                    config_file_path => $snpeff_config_file_path,
                    database_path    => $active_parameter_href->{snpsift_dbnsfp_file},
                    FILEHANDLE       => $XARGSFILEHANDLE,
                    infile_path      => $outfile_path{$contig}
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

    ## Move contig files into final outfile
    say {$FILEHANDLE} q{## Move contig files into final outfile};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Decrement to make last output file as input to renaming
    my $annotation_infile_number = $xargs_file_counter - 1;

  CONTIG:
    foreach my $contig ( keys %infile_path ) {

        gnu_mv(
            {
                FILEHANDLE  => $XARGSFILEHANDLE,
                infile_path => $outfile_path{$contig} . $DOT . $annotation_infile_number,
                stderrfile_path => $xargs_file_path_prefix
                  . $DOT
                  . $contig
                  . $DOT
                  . q{stderr.txt},
                outfile_path => $outfile_path{$contig}
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_paths[0],
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );
        submit_recipe(
            {
                base_command            => $profile_base_command,
                case_id                 => $case_id,
                dependency_method       => q{sample_to_case},
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                job_id_chain            => $job_id_chain,
                recipe_file_path        => $recipe_file_path,
                sample_ids_ref          => \@{ $active_parameter_href->{sample_ids} },
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

1;
