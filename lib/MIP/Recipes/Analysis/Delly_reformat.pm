package MIP::Recipes::Analysis::Delly_reformat;

use 5.026;
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
    our $VERSION = 1.05;

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

## Function : Merge, regenotype, and filter using Delly version 0.7.8
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $reference_dir           => MIP reference directory
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
    my $reference_dir;
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
        reference_dir_ref => {
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
        xargs_file_counter => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Delete::List qw{ delete_contig_elements };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_parameters get_recipe_attributes };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_merge bcftools_index bcftools_view };
    use MIP::Program::Variantcalling::Delly qw{ delly_call delly_merge };
    use MIP::Program::Variantcalling::Picardtools qw{ picardtools_sortvcf };
    use MIP::Processmanagement::Processes qw{ print_wait submit_recipe };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::QC::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constants
    Readonly my $SV_MAX_SIZE => 100_000_000;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger( uc q{mip_analyse} );

    ## Unpack parameters
    my $job_id_chain = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode        = $active_parameter_href->{$recipe_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set and get the io files per chain, id and stream
    my %io = parse_io_outfiles(
        {
            chain_id               => $job_id_chain,
            id                     => $case_id,
            file_info_href         => $file_info_href,
            file_name_prefixes_ref => [$case_id],
            outdata_dir            => $active_parameter_href->{outdata_dir},
            parameter_href         => $parameter_href,
            recipe_name            => $recipe_name,
            temp_directory         => $temp_directory,
        }
    );

    my $outdir_path_prefix       = $io{out}{dir_path_prefix};
    my $outfile_path_prefix      = $io{out}{file_path_prefix};
    my $outfile_suffix           = $io{out}{file_suffix};
    my $outfile_path             = $outfile_path_prefix . $outfile_suffix;
    my $temp_outfile_path_prefix = $io{temp}{file_path_prefix};
    my $temp_outfile_suffix      = $io{temp}{file_suffix};
    my $temp_outfile_path        = $temp_outfile_path_prefix . $temp_outfile_suffix;

    ## Filehandles
    # Create anonymous filehandles
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
            process_time                    => $time,
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    ## Collect infiles for dependence recipes streams for all sample_ids
    my %recipe_tag_keys = (
        gatk_baserecalibration => q{out},
        $recipe_name           => q{in},
    );

    my %delly_sample_file_info;
    my $process_batches_count = 1;
    while ( my ( $sample_id_index, $sample_id ) =
        each @{ $active_parameter_href->{sample_ids} } )
    {

      PROGRAM_TAG:
        while ( my ( $recipe_tag, $stream ) = each %recipe_tag_keys ) {

            ## Get the io infiles per chain and id
            my %sample_io = get_io_files(
                {
                    id             => $sample_id,
                    file_info_href => $file_info_href,
                    parameter_href => $parameter_href,
                    recipe_name    => $recipe_tag,
                    stream         => $stream,
                    temp_directory => $temp_directory,
                }
            );
            my $infile_path_prefix = $sample_io{$stream}{file_path_prefix};
            my $infile_suffix      = $sample_io{$stream}{file_suffix};
            my $infile_path =
              $infile_path_prefix . substr( $infile_suffix, 0, 2 ) . $ASTERISK;
            my $temp_infile_path_prefix = $sample_io{temp}{file_path_prefix};
            my $temp_infile_path        = $temp_infile_path_prefix . $infile_suffix;

            $delly_sample_file_info{$sample_id}{in}{$infile_suffix} =
              $temp_infile_path;

            $process_batches_count = print_wait(
                {
                    FILEHANDLE            => $FILEHANDLE,
                    max_process_number    => $core_number,
                    process_batches_count => $process_batches_count,
                    process_counter       => $sample_id_index,
                }
            );

            ## Copy file(s) to temporary directory
            say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
            migrate_file(
                {
                    FILEHANDLE   => $FILEHANDLE,
                    infile_path  => $infile_path,
                    outfile_path => $temp_directory,
                }
            );
        }
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Delly call bcf sample infiles
    my @delly_merge_temp_infile_paths =
      map { $delly_sample_file_info{$_}{in}{q{.bcf}} }
      @{ $active_parameter_href->{sample_ids} };

    ## We have something to merge
    if ( scalar @{ $active_parameter_href->{sample_ids} } > 1 ) {

        ### Delly merge
        say {$FILEHANDLE} q{## delly merge} . $NEWLINE;

        say {$FILEHANDLE} q{## Fix locale bug using old centosOS and Boost library};
        say {$FILEHANDLE} q?LC_ALL="C"; export LC_ALL ?, $NEWLINE . $NEWLINE;

        ## Get parameters
        my $xargs_file_path_prefix;

        delly_merge(
            {
                FILEHANDLE       => $FILEHANDLE,
                infile_paths_ref => \@delly_merge_temp_infile_paths,
                min_size         => 0,
                max_size         => $SV_MAX_SIZE,
                outfile_path     => $temp_outfile_path_prefix
                  . $UNDERSCORE
                  . q{merged}
                  . $DOT . q{bcf},
                stderrfile_path => $recipe_file_path
                  . $UNDERSCORE
                  . q{merged}
                  . $DOT
                  . q{stderr.txt},
                stdoutfile_path => $recipe_file_path
                  . $UNDERSCORE
                  . q{merged}
                  . $DOT
                  . q{stdout.txt},
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Delly call regenotype
        say {$FILEHANDLE} q{## delly call regenotype};

        ## Store outfiles
        my @delly_genotype_temp_outfile_paths;

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

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Assemble file path
            my $alignment_sample_file_path =
              $delly_sample_file_info{$sample_id}{in}{q{.bam}};
            my $bcf_sample_outfile_path =
                $temp_outfile_path_prefix
              . $UNDERSCORE
              . q{merged}
              . $UNDERSCORE . q{geno}
              . $UNDERSCORE
              . $sample_id
              . $DOT . q{bcf};
            push @delly_genotype_temp_outfile_paths, $bcf_sample_outfile_path;
            delly_call(
                {
                    exclude_file_path => $active_parameter_href->{delly_exclude_file},
                    FILEHANDLE        => $XARGSFILEHANDLE,
                    genotypefile_path => $temp_outfile_path_prefix
                      . $UNDERSCORE
                      . q{merged}
                      . $DOT . q{bcf},
                    infile_path        => $alignment_sample_file_path,
                    outfile_path       => $bcf_sample_outfile_path,
                    referencefile_path => $referencefile_path,
                    stderrfile_path    => $xargs_file_path_prefix
                      . $UNDERSCORE
                      . $sample_id
                      . $DOT
                      . q{stderr.txt},
                    stdoutfile_path => $xargs_file_path_prefix
                      . $UNDERSCORE
                      . $sample_id
                      . $DOT
                      . q{stdout.txt},
                }
            );
            say {$XARGSFILEHANDLE} $NEWLINE;
        }

        close $XARGSFILEHANDLE
          or $log->logcroak(q{Could not close XARGSFILEHANDLE});

        ### Merge calls
        say {$FILEHANDLE} q{## bcftools merge};

        bcftools_merge(
            {
                FILEHANDLE       => $FILEHANDLE,
                infile_paths_ref => \@delly_genotype_temp_outfile_paths,
                outfile_path     => $temp_outfile_path_prefix
                  . $UNDERSCORE
                  . q{to_sort}
                  . $outfile_suffix,
                output_type     => q{v},
                stderrfile_path => $xargs_file_path_prefix . $DOT . q{stderr.txt},
                stdoutfile_path => $xargs_file_path_prefix . $DOT . q{stdout.txt},
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }
    else {

        # Only one sample
        say {$FILEHANDLE} q{## Only one sample - skip merging and regenotyping};
        say {$FILEHANDLE}
q{## Reformat bcf infile to match outfile from regenotyping with multiple samples};

        bcftools_view(
            {
                FILEHANDLE   => $FILEHANDLE,
                output_type  => q{v},
                infile_path  => $delly_merge_temp_infile_paths[0],
                outfile_path => $temp_outfile_path_prefix
                  . $UNDERSCORE
                  . q{to_sort}
                  . $outfile_suffix,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Writes sbatch code to supplied filehandle to sort variants in vcf format
    say {$FILEHANDLE} q{## Picard SortVcf};
    picardtools_sortvcf(
        {
            FILEHANDLE => $FILEHANDLE,
            infile_paths_ref =>
              [ $temp_outfile_path_prefix . $UNDERSCORE . q{to_sort} . $outfile_suffix ],
            java_jar =>
              catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
            java_use_large_pages => $active_parameter_href->{java_use_large_pages},
            memory_allocation    => q{Xmx2g},
            outfile_path         => $temp_outfile_path_prefix . $DOT . q{vcf},
            referencefile_path   => $referencefile_path,
            sequence_dictionary  => catfile(
                $reference_dir,
                $file_info_href->{human_genome_reference_name_prefix} . $DOT . q{dict}
            ),
            temp_directory => $temp_directory,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} $NEWLINE . q{## Copy file from temporary directory};
    migrate_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $temp_outfile_path,
            outfile_path => $outdir_path_prefix,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $recipe_mode == 1 ) {

        set_recipe_outfile_in_sample_info(
            {
                recipe_name      => q{delly},
                path             => $outfile_path,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command            => $profile_base_command,
                dependency_method       => q{sample_to_case},
                case_id                 => $case_id,
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
