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
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_delly_reformat analysis_delly_reformat_old };

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
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
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
    my $parameter_href;
    my $program_name;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
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
    use MIP::Get::Parameter qw{ get_module_parameters get_program_attributes };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_merge bcftools_index bcftools_view };
    use MIP::Program::Variantcalling::Delly qw{ delly_call delly_merge };
    use MIP::Program::Variantcalling::Picardtools qw{ picardtools_sortvcf };
    use MIP::Processmanagement::Processes qw{ print_wait };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constants
    Readonly my $SV_MAX_SIZE => 100_000_000;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    my $job_id_chain = get_program_attributes(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
            attribute      => q{chain},
        }
    );
    my $program_mode       = $active_parameter_href->{$program_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Set and get the io files per chain, id and stream
    my %io = parse_io_outfiles(
        {
            chain_id               => $job_id_chain,
            id                     => $family_id,
            file_info_href         => $file_info_href,
            file_name_prefixes_ref => [$family_id],
            outdata_dir            => $active_parameter_href->{outdata_dir},
            parameter_href         => $parameter_href,
            program_name           => $program_name,
            temp_directory         => $temp_directory,
        }
    );

    my $outdir_path_prefix       = $io{out}{dir_path_prefix};
    my $outfile_path_prefix      = $io{out}{file_path_prefix};
    my $outfile_suffix           = $io{out}{file_suffix};
    my $outfile_path             = $outfile_path_prefix . $outfile_suffix;
    my $temp_outfile_path_prefix = $io{temp}{file_path_prefix};
    my $temp_outfile_suffix      = $io{temp}{file_suffix};
    my $temp_outfile_path = $temp_outfile_path_prefix . $temp_outfile_suffix;

    ## Filehandles
    # Create anonymous filehandles
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

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
            program_directory               => $program_name,
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    ## Collect infiles for dependence programs streams for all sample_ids
    my %program_tag_keys = (
        gatk_baserecalibration => q{out},
        $program_name          => q{in},
    );

    my %delly_sample_file_info;
    my $process_batches_count = 1;
    while ( my ( $sample_id_index, $sample_id ) =
        each @{ $active_parameter_href->{sample_ids} } )
    {

      PROGRAM_TAG:
        while ( my ( $program_tag, $stream ) = each %program_tag_keys ) {

            ## Get the io infiles per chain and id
            my %sample_io = get_io_files(
                {
                    id             => $sample_id,
                    file_info_href => $file_info_href,
                    parameter_href => $parameter_href,
                    program_name   => $program_tag,
                    stream         => $stream,
                    temp_directory => $temp_directory,
                }
            );
            my $infile_path_prefix = $sample_io{$stream}{file_path_prefix};
            my $infile_suffix      = $sample_io{$stream}{file_suffix};
            my $infile_path =
              $infile_path_prefix . substr( $infile_suffix, 0, 2 ) . $ASTERISK;
            my $temp_infile_path_prefix = $sample_io{temp}{file_path_prefix};
            my $temp_infile_path = $temp_infile_path_prefix . $infile_suffix;

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

        say {$FILEHANDLE}
          q{## Fix locale bug using old centosOS and Boost library};
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
                stderrfile_path => $file_path
                  . $UNDERSCORE
                  . q{merged}
                  . $DOT
                  . q{stderr.txt},
                stdoutfile_path => $file_path
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
                file_path          => $file_path,
                program_info_path  => $program_info_path,
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
                    exclude_file_path =>
                      $active_parameter_href->{delly_exclude_file},
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
                stderrfile_path => $xargs_file_path_prefix
                  . $DOT
                  . q{stderr.txt},
                stdoutfile_path => $xargs_file_path_prefix
                  . $DOT
                  . q{stdout.txt},
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
            FILEHANDLE       => $FILEHANDLE,
            infile_paths_ref => [
                    $temp_outfile_path_prefix
                  . $UNDERSCORE
                  . q{to_sort}
                  . $outfile_suffix
            ],
            java_jar => catfile(
                $active_parameter_href->{picardtools_path},
                q{picard.jar}
            ),
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            memory_allocation   => q{Xmx2g},
            outfile_path        => $temp_outfile_path_prefix . $DOT . q{vcf},
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
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $temp_outfile_path,
            outfile_path => $outdir_path_prefix,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $program_mode == 1 ) {

        add_program_outfile_to_sample_info(
            {
                program_name     => q{delly},
                path             => $outfile_path,
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

sub analysis_delly_reformat_old {

## Function : Merge, regenotype, and filter using Delly
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => Variant call type
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
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
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        call_type =>
          { default => q{SV}, store => \$call_type, strict_type => 1, },
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
    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Get::Parameter qw{ get_module_parameters };
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

    ## Set program name
    my $program_mode = $active_parameter_href->{$program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$program_name}{chain};
    my $program_outdirectory_name =
      $parameter_href->{$program_name}{outdir_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Filehandles
    # Create anonymous filehandles
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    my $program_directory =
      catfile( $outaligner_dir, $program_outdirectory_name );

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
            program_directory               => $program_directory,
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ## Assign directories
    my $outfamily_directory = catfile(
        $active_parameter_href->{outdata_dir},
        $active_parameter_href->{family_id},
        $active_parameter_href->{outaligner_dir},
        $program_outdirectory_name,
    );

    # Used downstream
    $parameter_href->{$program_name}{indirectory} = $outfamily_directory;

    ## Tags
    my $outfile_tag =
      $file_info_href->{$family_id}{$program_name}{file_tag};

    ### Assign suffix
    ## Files
    my $outfile_prefix = $family_id . $outfile_tag . $UNDERSCORE . $call_type;
    my $file_suffix    = $parameter_href->{$program_name}{outfile_suffix};

    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            file_suffix    => $file_suffix,
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

    ## Paths
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

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
    my @program_tag_keys = (qw{ gatk_baserecalibration delly_call });

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Assign directories
        my $insample_directory_bam =
          catdir( $active_parameter_href->{outdata_dir},
            $sample_id, $outaligner_dir );

        my $insample_directory_bcf =
          catdir( $active_parameter_href->{outdata_dir},
            $sample_id, $outaligner_dir,
            $parameter_href->{delly_call}{outdir_name} );

      INFILE_TAG:
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
            if ( $infile_tag_key eq q{delly_call} ) {

                ## Assign suffix
                $suffix{$infile_tag_key} = get_file_suffix(
                    {
                        parameter_href => $parameter_href,
                        program_name   => $infile_tag_key,
                        suffix_key     => q{outfile_suffix},
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
                                contigs_ref => \@contigs,
                                core_number => $core_number,
                                FILEHANDLE  => $FILEHANDLE,
                                file_ending => $file_ending,
                                file_path   => $file_path,
                                indirectory => $insample_directory_bcf,
                                infile => $merged_infile_prefix . $infile_tag,
                                program_info_path  => $program_info_path,
                                temp_directory     => $temp_directory,
                                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                                xargs_file_counter => $xargs_file_counter,
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
            else {

                #BAMs
                ## Assign suffix
                $suffix{$infile_tag_key} = get_file_suffix(
                    {
                        jobid_chain =>
                          $parameter_href->{$infile_tag_key}{chain},
                        parameter_href => $parameter_href,
                        suffix_key     => q{alignment_file_suffix},
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
                        contigs_ref => \@contigs,
                        core_number => ( $core_number - 1 )
                        ,    # Compensate for cp of entire BAM TRA, see above
                        FILEHANDLE  => $FILEHANDLE,
                        file_ending => $file_ending,
                        file_path   => $file_path,
                        indirectory => $insample_directory_bam,
                        infile      => $merged_infile_prefix . $infile_tag,
                        program_info_path  => $program_info_path,
                        temp_directory     => $temp_directory,
                        XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                        xargs_file_counter => $xargs_file_counter,
                    }
                );
            }
            say {$FILEHANDLE} q{wait}, $NEWLINE;
        }
    }

    if ( scalar @{ $active_parameter_href->{sample_ids} } > 1 ) {

        ### Delly merge
        say {$FILEHANDLE} q{## delly merge} . $NEWLINE;

        say {$FILEHANDLE}
          q{## Fix locale bug using old centosOS and Boost library};
        say {$FILEHANDLE} q?LC_ALL="C"; export LC_ALL ?, $NEWLINE . $NEWLINE;

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

      SV_TYPE:
        foreach my $sv_type ( @{ $active_parameter_href->{delly_types} } ) {

            my $sv_max_size = 100_000_000;

            if ( $sv_type ne q{TRA} ) {

              CONTIG:
                foreach my $contig (@contigs) {

                    ## Assemble file paths by adding file ending
                    my @file_paths = map {
                            $infile_path_prefix{$_}{delly_call}
                          . $UNDERSCORE
                          . $contig
                          . $UNDERSCORE
                          . $sv_type
                          . $suffix{delly_call}
                    } @{ $active_parameter_href->{sample_ids} };

                    delly_merge(
                        {
                            FILEHANDLE       => $XARGSFILEHANDLE,
                            infile_paths_ref => \@file_paths,
                            min_size         => 0,
                            max_size         => $sv_max_size,
                            outfile_path     => $outfile_path_prefix
                              . $UNDERSCORE
                              . $contig
                              . $UNDERSCORE
                              . $sv_type
                              . $DOT . q{bcf},
                            stderrfile_path => $xargs_file_path_prefix
                              . $DOT
                              . $contig
                              . $DOT
                              . $sv_type
                              . $DOT
                              . q{stderr.txt},
                            stdoutfile_path => $xargs_file_path_prefix
                              . $DOT
                              . $contig
                              . $DOT
                              . $sv_type
                              . $DOT
                              . q{stdout.txt},
                            sv_type => $sv_type,
                        }
                    );
                    say {$XARGSFILEHANDLE} $NEWLINE;
                }
            }
            else {

                ## Assemble file paths by adding file ending
                my @file_paths = map {
                        $infile_path_prefix{$_}{delly_call}
                      . $UNDERSCORE
                      . $sv_type
                      . $suffix{delly_call}
                } @{ $active_parameter_href->{sample_ids} };

                delly_merge(
                    {
                        FILEHANDLE       => $XARGSFILEHANDLE,
                        infile_paths_ref => \@file_paths,
                        min_size         => 0,
                        max_size         => $sv_max_size,
                        outfile_path     => $outfile_path_prefix
                          . $UNDERSCORE
                          . $sv_type
                          . $DOT . q{bcf},
                        stderrfile_path => $xargs_file_path_prefix
                          . $DOT
                          . $sv_type
                          . $DOT
                          . q{stderr.txt},
                        stdoutfile_path => $xargs_file_path_prefix
                          . $DOT
                          . $sv_type
                          . $DOT
                          . q{stdout.txt},
                        sv_type => $sv_type,
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
                    core_number        => $core_number,
                    FILEHANDLE         => $FILEHANDLE,
                    file_path          => $file_path,
                    program_info_path  => $program_info_path,
                    XARGSFILEHANDLE    => $XARGSFILEHANDLE,
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
                          {gatk_baserecalibration}
                          . $UNDERSCORE
                          . $contig
                          . $suffix{gatk_baserecalibration};

                        delly_call(
                            {
                                exclude_file_path =>
                                  $active_parameter_href->{delly_exclude_file},
                                FILEHANDLE        => $XARGSFILEHANDLE,
                                genotypefile_path => $outfile_path_prefix
                                  . $UNDERSCORE
                                  . $contig
                                  . $UNDERSCORE
                                  . $sv_type
                                  . $suffix{delly_call},
                                infile_path  => $alignment_sample_file_path,
                                outfile_path => $file_path_prefix{$sample_id}
                                  . $UNDERSCORE
                                  . $contig
                                  . $UNDERSCORE
                                  . $sv_type
                                  . $UNDERSCORE . q{geno}
                                  . $suffix{delly_call},
                                referencefile_path => $referencefile_path,
                                stderrfile_path    => $xargs_file_path_prefix
                                  . $DOT
                                  . $contig
                                  . $DOT
                                  . $sv_type
                                  . $DOT
                                  . q{stderr.txt},
                                stdoutfile_path => $xargs_file_path_prefix
                                  . $DOT
                                  . $contig
                                  . $DOT
                                  . $sv_type
                                  . $DOT
                                  . q{stdout.txt},
                                sv_type => $sv_type,
                            }
                        );
                        say {$XARGSFILEHANDLE} $NEWLINE;
                    }
                }
                else {

                    ## Assemble file path
                    my $alignment_sample_file_path =
                        $infile_path_prefix{$sample_id}{gatk_baserecalibration}
                      . $suffix{gatk_baserecalibration};

                    delly_call(
                        {
                            exclude_file_path =>
                              $active_parameter_href->{delly_exclude_file},
                            FILEHANDLE        => $XARGSFILEHANDLE,
                            genotypefile_path => $outfile_path_prefix
                              . $UNDERSCORE
                              . $sv_type
                              . $suffix{delly_call},
                            infile_path  => $alignment_sample_file_path,
                            outfile_path => $file_path_prefix{$sample_id}
                              . $UNDERSCORE
                              . $sv_type
                              . $UNDERSCORE . q{geno}
                              . $suffix{delly_call},
                            referencefile_path => $referencefile_path,
                            stderrfile_path    => $xargs_file_path_prefix
                              . $DOT
                              . $sv_type
                              . $DOT
                              . q{stderr.txt},
                            stdoutfile_path => $xargs_file_path_prefix
                              . $DOT
                              . $sv_type
                              . $DOT
                              . q{stdout.txt},
                            sv_type => $sv_type,
                        }
                    );
                    say {$XARGSFILEHANDLE} $NEWLINE;
                }
            }
        }

        ### Merge calls
        say {$FILEHANDLE} q{## bcftools merge};

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
                          . $suffix{delly_call}
                    } @{ $active_parameter_href->{sample_ids} };

                    bcftools_merge(
                        {
                            FILEHANDLE       => $XARGSFILEHANDLE,
                            infile_paths_ref => \@file_paths,
                            outfile_path     => $outfile_path_prefix
                              . $UNDERSCORE
                              . $contig
                              . $UNDERSCORE
                              . $sv_type
                              . $suffix{delly_call},
                            output_type     => q{b},
                            stderrfile_path => $xargs_file_path_prefix
                              . $DOT
                              . $contig
                              . $DOT
                              . $sv_type
                              . $DOT
                              . q{stderr.txt},
                            stdoutfile_path => $xargs_file_path_prefix
                              . $DOT
                              . $contig
                              . $DOT
                              . $sv_type
                              . $DOT
                              . q{stdout.txt},
                        }
                    );
                    print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

                    bcftools_index(
                        {
                            FILEHANDLE  => $XARGSFILEHANDLE,
                            infile_path => $outfile_path_prefix
                              . $UNDERSCORE
                              . $contig
                              . $UNDERSCORE
                              . $sv_type
                              . $suffix{delly_call},
                            output_type     => q{csi},
                            stderrfile_path => $xargs_file_path_prefix
                              . $DOT
                              . $contig
                              . $DOT
                              . $sv_type
                              . $UNDERSCORE
                              . q{index.stderr.txt},
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
                      . $suffix{delly_call}
                } @{ $active_parameter_href->{sample_ids} };

                bcftools_merge(
                    {
                        FILEHANDLE       => $XARGSFILEHANDLE,
                        infile_paths_ref => \@file_paths,
                        outfile_path     => $outfile_path_prefix
                          . $UNDERSCORE
                          . $sv_type
                          . $suffix{delly_call},
                        output_type     => q{b},
                        stderrfile_path => $xargs_file_path_prefix
                          . $DOT
                          . $sv_type
                          . $DOT
                          . q{stderr.txt},
                        stdoutfile_path => $xargs_file_path_prefix
                          . $DOT
                          . $sv_type
                          . $DOT
                          . q{stdout.txt},
                    }
                );
                print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

                bcftools_index(
                    {
                        FILEHANDLE  => $XARGSFILEHANDLE,
                        infile_path => $outfile_path_prefix
                          . $UNDERSCORE
                          . $sv_type
                          . $suffix{delly_call},
                        stderrfile_path => $xargs_file_path_prefix
                          . $DOT
                          . $sv_type
                          . $UNDERSCORE
                          . q{index.stderr.txt},
                        output_type => q{csi},
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
                      . $suffix{delly_call}
                } @contigs;
            }
            else {

                push @file_paths,
                    $outfile_path_prefix
                  . $UNDERSCORE
                  . $sv_type
                  . $suffix{delly_call};
            }

            bcftools_concat(
                {
                    allow_overlaps   => 1,
                    FILEHANDLE       => $FILEHANDLE,
                    infile_paths_ref => \@file_paths,
                    outfile_path     => $outfile_path_prefix
                      . $UNDERSCORE
                      . $sv_type
                      . $UNDERSCORE
                      . q{concat}
                      . $suffix{delly_call},
                    output_type     => q{b},
                    rm_dups         => q{all},
                    stderrfile_path => $program_info_path
                      . $UNDERSCORE
                      . $sv_type
                      . $UNDERSCORE
                      . q{concat.stderr.txt},
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            bcftools_index(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => $outfile_path_prefix
                      . $UNDERSCORE
                      . $sv_type
                      . $UNDERSCORE
                      . q{concat}
                      . $suffix{delly_call},
                    output_type     => q{csi},
                    stderrfile_path => $xargs_file_path_prefix
                      . $DOT
                      . $sv_type
                      . $UNDERSCORE
                      . q{index.stderr.txt},
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }
        ### Filter calls
        say {$FILEHANDLE} q{## Delly filter};

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
      SV_TYPE:
        foreach my $sv_type ( @{ $active_parameter_href->{delly_types} } ) {

            if ( $sv_type ne q{TRA} ) {

                delly_filter(
                    {
                        FILEHANDLE  => $XARGSFILEHANDLE,
                        filter_mode => q{germline},
                        infile_path => $outfile_path_prefix
                          . $UNDERSCORE
                          . $sv_type
                          . $UNDERSCORE
                          . q{concat}
                          . $suffix{delly_call},
                        outfile_path => $outfile_path_prefix
                          . $UNDERSCORE
                          . $sv_type
                          . $UNDERSCORE
                          . q{filtered}
                          . $suffix{delly_call},
                        stderrfile_path => $xargs_file_path_prefix
                          . $DOT
                          . $sv_type
                          . $DOT
                          . q{stderr.txt},
                        stdoutfile_path => $xargs_file_path_prefix
                          . $DOT
                          . $sv_type
                          . $DOT
                          . q{stdout.txt},
                        sv_type => $sv_type,

                    }
                );
                say {$XARGSFILEHANDLE} $NEWLINE;
            }
            else {

                delly_filter(
                    {
                        FILEHANDLE  => $XARGSFILEHANDLE,
                        filter_mode => q{germline},
                        infile_path => $outfile_path_prefix
                          . $UNDERSCORE
                          . $sv_type
                          . $suffix{delly_call},
                        outfile_path => $outfile_path_prefix
                          . $UNDERSCORE
                          . $sv_type
                          . $UNDERSCORE
                          . q{filtered}
                          . $suffix{delly_call},
                        stderrfile_path => $xargs_file_path_prefix
                          . $DOT
                          . $sv_type
                          . $DOT
                          . q{stderr.txt},
                        stdoutfile_path => $xargs_file_path_prefix
                          . $DOT
                          . $sv_type
                          . $DOT
                          . q{stdout.txt},
                        sv_type => $sv_type,
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
              . $suffix{delly_call}
        } @{ $active_parameter_href->{delly_types} };
    }
    else {    #Only one sample

        say {$FILEHANDLE} q{## Only one sample - skip merging and regenotyping};
        say {$FILEHANDLE}
          q{## bcftools concat - merge all SV types and contigs};

      SV_TYPE:
        foreach my $sv_type ( @{ $active_parameter_href->{delly_types} } ) {

            if ( $sv_type ne q{TRA} ) {

              CONTIG:
                foreach my $contig (@contigs) {

                    ## Assemble file paths by adding file ending
                    push @file_paths, map {
                            $infile_path_prefix{$_}{delly_call}
                          . $UNDERSCORE
                          . $contig
                          . $UNDERSCORE
                          . $sv_type
                          . $suffix{delly_call}
                    } @{ $active_parameter_href->{sample_ids} };
                }
            }
            else {

                push @file_paths, map {
                        $infile_path_prefix{$_}{delly_call}
                      . $UNDERSCORE
                      . $sv_type
                      . $suffix{delly_call}
                } @{ $active_parameter_href->{sample_ids} };
            }
        }
    }
    bcftools_concat(
        {
            allow_overlaps   => 1,
            FILEHANDLE       => $FILEHANDLE,
            infile_paths_ref => \@file_paths,
            outfile_path     => $outfile_path_prefix
              . $UNDERSCORE
              . q{concat}
              . $outfile_suffix,
            output_type     => q{v},
            rm_dups         => q{all},
            stderrfile_path => $program_info_path
              . $UNDERSCORE
              . q{concat.stderr.txt},
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
            memory_allocation   => q{Xmx2g},
            outfile_path        => $outfile_path_prefix . $outfile_suffix,
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
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $outfile_path_prefix . $outfile_suffix,
            outfile_path => $outfamily_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    if ( $program_mode == 1 ) {

        add_program_outfile_to_sample_info(
            {
                program_name => q{delly},
                path         => catfile(
                    $outfamily_directory, $outfile_prefix . $outfile_suffix
                ),
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

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});
    close $XARGSFILEHANDLE
      or $log->logcroak(q{Could not close $XARGSFILEHANDLE});
    return;
}

1;
