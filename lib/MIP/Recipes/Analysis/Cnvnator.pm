package MIP::Recipes::Analysis::Cnvnator;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw{ first_value };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $AMPERSAND $ASTERISK $DOT $EMPTY_STR $LOG_NAME $NEWLINE $SEMICOLON $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_cnvnator };

}

sub analysis_cnvnator {

## Function : Call structural variants using cnvnator
## Returns  :
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => The file_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $reference_dir           => MIP reference directory
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $sample_id               => Sample id
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_id;
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
        reference_dir => {
            default     => $arg_href->{active_parameter_href}{reference_dir},
            store       => \$reference_dir,
            strict_type => 1,
        },
        sample_id => {
            default     => 1,
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
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
            allow       => qr/ ^\d+$ /sxm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Contigs qw{ delete_contig_elements };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Samtools qw{ samtools_create_chromosome_files };
    use MIP::Program::Bcftools
      qw{ bcftools_annotate bcftools_concat bcftools_create_reheader_samples_file bcftools_rename_vcf_samples };
    use MIP::Program::Cnvnator
      qw{ cnvnator_read_extraction cnvnator_histogram cnvnator_statistics cnvnator_partition cnvnator_calling cnvnator_convert_to_vcf };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my %infile_path        = %{ $io{in}{file_path_href} };

    my $human_genome_reference = $active_parameter_href->{human_genome_reference};
    my $job_id_chain           = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number = $recipe_resource{core_number};
    my $xargs_file_path_prefix;

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => [$infile_name_prefix],
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );

    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my $outfile_path        = $io{out}{file_path};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle      = IO::Handle->new();
    my $xargsfilehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $sample_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe_resource{memory},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            process_time          => $recipe_resource{time},
            temp_directory        => $temp_directory,
        }
    );

    ### SHELL:

    ### Update contigs
    ## Removes an element from array and return new array while leaving orginal contigs_ref untouched
    # Skip contig MT as it is not applicable in this analysis
    my @size_ordered_contigs = delete_contig_elements(
        {
            contigs_ref        => \@{ $file_info_href->{contigs_size_ordered} },
            remove_contigs_ref => [qw{ MT M }],
        }
    );
    my @contigs = delete_contig_elements(
        {
            contigs_ref        => \@{ $file_info_href->{contigs} },
            remove_contigs_ref => [qw{ MT M }],
        }
    );

    ## Add contigs to vcfheader
    _add_contigs_to_vcfheader(
        {
            filehandle             => $filehandle,
            human_genome_reference => $human_genome_reference,
            temp_directory         => $outdir_path_prefix,
        }
    );

    ## Call to Samtools to create chromosome sequence files (.fa) used by Cnvnator
    say {$filehandle} q{## Create by cnvnator required 'chr.fa' files};
    samtools_create_chromosome_files(
        {
            filehandle         => $filehandle,
            infile_path        => $human_genome_reference,
            max_process_number => $core_number,
            regions_ref        => \@size_ordered_contigs,
            suffix             => $DOT . q{fa},
            temp_directory     => $outdir_path_prefix,
        }
    );
    say {$filehandle} q{wait}, $NEWLINE;

    ## cnvnator
    say {$filehandle} q{## cnvnator};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            filehandle         => $filehandle,
            xargsfilehandle    => $xargsfilehandle,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig (@size_ordered_contigs) {

        my $stdbasefile_path_prefix = $xargs_file_path_prefix . $DOT . $contig;

        ## Assemble parameter
        # Output ROOT file
        my $root_file                  = $infile_path{$contig} . $DOT . q{root};
        my $outfile_path_prefix_contig = $outfile_path_prefix . $UNDERSCORE . $contig;

        cnvnator_read_extraction(
            {
                infile_paths_ref => [ $infile_path{$contig} ],
                outfile_path     => $root_file,
                regions_ref      => [$contig],
                stdoutfile_path  => $stdbasefile_path_prefix . $DOT . q{stdout.txt},
                stderrfile_path  => $stdbasefile_path_prefix . $DOT . q{stderr.txt},
                filehandle       => $xargsfilehandle,
            }
        );
        print {$xargsfilehandle} $SEMICOLON . $SPACE;

        cnvnator_histogram(
            {
                infile_path             => $root_file,
                regions_ref             => [$contig],
                cnv_bin_size            => $active_parameter_href->{cnv_bin_size},
                referencedirectory_path => $outdir_path_prefix,
                filehandle              => $xargsfilehandle,
                stdoutfile_path         => $stdbasefile_path_prefix
                  . $UNDERSCORE
                  . q{histogram.stdout.txt},
                stderrfile_path => $stdbasefile_path_prefix
                  . $UNDERSCORE
                  . q{histogram.stderr.txt},
            }
        );
        print {$xargsfilehandle} $SEMICOLON . $SPACE;

        cnvnator_statistics(
            {
                infile_path     => $root_file,
                regions_ref     => [$contig],
                cnv_bin_size    => $active_parameter_href->{cnv_bin_size},
                filehandle      => $xargsfilehandle,
                stdoutfile_path => $stdbasefile_path_prefix
                  . $UNDERSCORE
                  . q{statistics.stdout.txt},
                stderrfile_path => $stdbasefile_path_prefix
                  . $UNDERSCORE
                  . q{statistics.stderr.txt},
            }
        );
        print {$xargsfilehandle} $SEMICOLON . $SPACE;

        cnvnator_partition(
            {
                infile_path     => $root_file,
                regions_ref     => [$contig],
                cnv_bin_size    => $active_parameter_href->{cnv_bin_size},
                filehandle      => $xargsfilehandle,
                stdoutfile_path => $stdbasefile_path_prefix
                  . $UNDERSCORE
                  . q{partition.stdout.txt},
                stderrfile_path => $stdbasefile_path_prefix
                  . $UNDERSCORE
                  . q{partition.stderr.txt},
            }
        );
        print {$xargsfilehandle} $SEMICOLON . $SPACE;

        cnvnator_calling(
            {
                infile_path     => $root_file,
                stdoutfile_path => $outfile_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $DOT
                  . q{cnvnator},
                regions_ref     => [$contig],
                cnv_bin_size    => $active_parameter_href->{cnv_bin_size},
                filehandle      => $xargsfilehandle,
                stderrfile_path => $stdbasefile_path_prefix
                  . $UNDERSCORE
                  . q{calling.stderr.txt},
            }
        );
        print {$xargsfilehandle} $SEMICOLON . $SPACE;

        cnvnator_convert_to_vcf(
            {
                infile_path     => $outfile_path_prefix_contig . $DOT . q{cnvnator},
                stdoutfile_path => $outfile_path_prefix_contig . $outfile_suffix,
                stderrfile_path => $stdbasefile_path_prefix
                  . $UNDERSCORE
                  . q{convert_to_vcf.stderr.txt},
                referencedirectory_path => $outdir_path_prefix,
                filehandle              => $xargsfilehandle,
            }
        );
        say {$xargsfilehandle} $NEWLINE;
    }

    ## Write sbatch code to supplied filehandle to concatenate variants in vcf format.
    say {$filehandle} q{## Format the VCF};

    ## Store infiles for bcftools concat
    my @concat_infile_paths;

    say {$filehandle} q{## Adding sample id instead of file prefix};
    my $rename_sample = first_value { $sample_id eq $_ }
    @{ $active_parameter_href->{sample_ids} };

    bcftools_create_reheader_samples_file(
        {
            filehandle     => $filehandle,
            sample_ids_ref => [$rename_sample],
            temp_directory => $outdir_path_prefix,
        }
    );

  CONTIG:
    foreach my $contig (@contigs) {

        ## Name intermediary files
        my $cnvnator_outfile_path =
          $outfile_path_prefix . $UNDERSCORE . $contig . $outfile_suffix;
        my $fixed_vcffile_path_prefix =
          $outfile_path_prefix . $UNDERSCORE . $contig . $UNDERSCORE . q{fixed};
        my $fixed_header_vcffile_path =
            $outfile_path_prefix
          . $UNDERSCORE
          . $contig
          . $UNDERSCORE
          . q{annot}
          . $outfile_suffix;

        ## Save infiles for bcftools annotate
        push @concat_infile_paths, $fixed_header_vcffile_path;

        bcftools_rename_vcf_samples(
            {
                create_sample_file  => 0,
                filehandle          => $filehandle,
                index               => 0,
                infile              => $cnvnator_outfile_path,
                outfile_path_prefix => $fixed_vcffile_path_prefix,
                output_type         => q{v},
                temp_directory      => $outdir_path_prefix,
                sample_ids_ref      => [$rename_sample],
            }
        );

        ## Add contigs to header
        bcftools_annotate(
            {
                filehandle      => $filehandle,
                headerfile_path => catfile( $outdir_path_prefix, q{contig_header.txt} ),
                infile_path     => $fixed_vcffile_path_prefix . $outfile_suffix,
                outfile_path    => $fixed_header_vcffile_path,
                output_type     => q{v},
            }
        );
        say {$filehandle} $NEWLINE;
    }

    say {$filehandle} q{## Concatenate VCFs};
    bcftools_concat(
        {
            filehandle       => $filehandle,
            infile_paths_ref => \@concat_infile_paths,
            outfile_path     => $outfile_path,
            output_type      => q{v},
            rm_dups          => 0,
        }
    );
    say {$filehandle} $NEWLINE;

    close $filehandle;

    if ( $recipe_mode == 1 ) {

        set_recipe_outfile_in_sample_info(
            {
                sample_info_href => $sample_info_href,
                recipe_name      => q{cnvnator},
                path             => $outfile_path,
            }
        );

        submit_recipe(
            {
                base_command      => $profile_base_command,
                case_id           => $case_id,
                dependency_method => q{sample_to_sample},
                job_id_chain      => $job_id_chain,
                job_id_href       => $job_id_href,
                log               => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path     => $recipe_file_path,
                job_reservation_name => $active_parameter_href->{job_reservation_name},
                sample_id            => $sample_id,
                submission_profile   => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub _add_contigs_to_vcfheader {

## Function : Change the sample ID in the VCF header
##          : $filehandle             => Filehandle to write to
##          : $human_genome_reference => Human genome reference
##          : $temp_directory         => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $human_genome_reference;
    my $temp_directory;

    my $tmpl = {
        filehandle => {
            defined  => 1,
            required => 1,
            store    => \$filehandle,
        },
        human_genome_reference => {
            defined  => 1,
            required => 1,
            store    => \$human_genome_reference,
        },
        temp_directory => {
            defined  => 1,
            required => 1,
            store    => \$temp_directory,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Executable qw{ get_executable_base_command };

    my $regexp;

    my @commands = ( get_executable_base_command( { base_command => q{perl}, } ), );

    # Execute perl
    print {$filehandle} join $SPACE, @commands;
    $regexp = q? -nae '?;

    ## Write contig ID from header
    $regexp .= q?{print "##contig=<ID=".$F[0].",length=".$F[1].">", "\n"}'?;

    print {$filehandle} $regexp . $SPACE;

    # Reference fai file
    print {$filehandle} $human_genome_reference . $DOT . q{fai} . $SPACE;

    say {$filehandle} q{>} . $SPACE . catfile( $temp_directory, q{contig_header.txt} ),
      $NEWLINE;

    return;
}

1;
