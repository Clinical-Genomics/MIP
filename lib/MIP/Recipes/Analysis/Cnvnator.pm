package MIP::Recipes::Analysis::Cnvnator;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw{ first_value };
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.08;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_cnvnator };

}

## Constants
Readonly my $AMPERSAND  => q{&};
Readonly my $ASTERISK   => q{*};
Readonly my $DOT        => q{.};
Readonly my $EMPTY_STR  => q{};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SEMICOLON  => q{;};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_cnvnator {

## Function : Call structural variants using cnvnator
## Returns  :
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => The file_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $reference_dir           => MIP reference directory
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $sample_id               => Sample id
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
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
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

    use MIP::Delete::List qw{ delete_contig_elements };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_module_parameters get_program_attributes };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Alignment::Samtools qw{ samtools_create_chromosome_files };
    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_annotate bcftools_concat bcftools_create_reheader_samples_file bcftools_rename_vcf_samples };
    use MIP::Program::Variantcalling::Cnvnator
      qw{ cnvnator_read_extraction cnvnator_histogram cnvnator_statistics cnvnator_partition cnvnator_calling cnvnator_convert_to_vcf };
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
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            program_name   => $program_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my $indir_path_prefix  = $io{in}{dir_path_prefix};
    my $infile_suffix      = $io{in}{file_suffix};
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my %temp_infile_path   = %{ $io{temp}{file_path_href} };

    my $human_genome_reference = $active_parameter_href->{human_genome_reference};
    my $job_id_chain           = get_program_attributes(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
            attribute      => q{chain},
        }
    );
    my $program_mode = $active_parameter_href->{$program_name};
    my ( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
    );
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
                program_name           => $program_name,
                temp_directory         => $temp_directory,
            }
        )
    );

    my $outdir_path_prefix       = $io{out}{dir_path_prefix};
    my $outfile_path_prefix      = $io{out}{file_path_prefix};
    my $outfile_suffix           = $io{out}{file_suffix};
    my $outfile_path             = $outfile_path_prefix . $outfile_suffix;
    my $temp_outfile_path_prefix = $io{temp}{file_path_prefix};
    my $temp_outfile_path        = $temp_outfile_path_prefix . $outfile_suffix;

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $recipe_file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $sample_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            program_directory               => $program_name,
            program_name                    => $program_name,
            process_time                    => $time,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    ### Update contigs
    ## Removes an element from array and return new array while leaving orginal elements_ref untouched
    # Skip contig MT as it is not applicable in this analysis
    my @size_ordered_contigs = delete_contig_elements(
        {
            elements_ref       => \@{ $file_info_href->{contigs_size_ordered} },
            remove_contigs_ref => [qw{ MT M }],
        }
    );
    my @contigs = delete_contig_elements(
        {
            elements_ref       => \@{ $file_info_href->{contigs} },
            remove_contigs_ref => [qw{ MT M }],
        }
    );

    ## Add contigs to vcfheader
    _add_contigs_to_vcfheader(
        {
            FILEHANDLE             => $FILEHANDLE,
            human_genome_reference => $human_genome_reference,
            temp_directory         => $temp_directory,
        }
    );

    ## Call to Samtools to create chromosome sequence files (.fa) used by Cnvnator
    say {$FILEHANDLE} q{## Create by cnvnator required 'chr.fa' files};
    samtools_create_chromosome_files(
        {
            FILEHANDLE         => $FILEHANDLE,
            infile_path        => $human_genome_reference,
            max_process_number => $core_number,
            regions_ref        => \@size_ordered_contigs,
            suffix             => $DOT . q{fa},
            temp_directory     => $temp_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    ($xargs_file_counter) = xargs_migrate_contig_files(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            contigs_ref        => \@size_ordered_contigs,
            file_path          => $recipe_file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
            infile             => $infile_name_prefix,
            indirectory        => $indir_path_prefix,
            file_ending        => substr( $infile_suffix, 0, 2 ) . $ASTERISK,
            temp_directory     => $temp_directory,
        }
    );

    ## cnvnator
    say {$FILEHANDLE} q{## cnvnator};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            file_path          => $recipe_file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig (@size_ordered_contigs) {

        my $stdbasefile_path_prefix = $xargs_file_path_prefix . $DOT . $contig;

        ## Assemble parameter
        # Output ROOT file
        my $root_file = $temp_infile_path{$contig} . $DOT . q{root};
        my $temp_outfile_path_prefix_contig =
          $temp_outfile_path_prefix . $UNDERSCORE . $contig;

        cnvnator_read_extraction(
            {
                infile_paths_ref => [ $temp_infile_path{$contig} ],
                outfile_path     => $root_file,
                regions_ref      => [$contig],
                unique           => 1,
                stdoutfile_path  => $stdbasefile_path_prefix . $DOT . q{stdout.txt},
                stderrfile_path  => $stdbasefile_path_prefix . $DOT . q{stderr.txt},
                FILEHANDLE       => $XARGSFILEHANDLE,
            }
        );
        print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

        cnvnator_histogram(
            {
                infile_path             => $root_file,
                regions_ref             => [$contig],
                cnv_bin_size            => $active_parameter_href->{cnv_bin_size},
                referencedirectory_path => $temp_directory,
                FILEHANDLE              => $XARGSFILEHANDLE,
                stdoutfile_path         => $stdbasefile_path_prefix
                  . $UNDERSCORE
                  . q{histogram.stdout.txt},
                stderrfile_path => $stdbasefile_path_prefix
                  . $UNDERSCORE
                  . q{histogram.stderr.txt},
            }
        );
        print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

        cnvnator_statistics(
            {
                infile_path     => $root_file,
                regions_ref     => [$contig],
                cnv_bin_size    => $active_parameter_href->{cnv_bin_size},
                FILEHANDLE      => $XARGSFILEHANDLE,
                stdoutfile_path => $stdbasefile_path_prefix
                  . $UNDERSCORE
                  . q{statistics.stdout.txt},
                stderrfile_path => $stdbasefile_path_prefix
                  . $UNDERSCORE
                  . q{statistics.stderr.txt},
            }
        );
        print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

        cnvnator_partition(
            {
                infile_path     => $root_file,
                regions_ref     => [$contig],
                cnv_bin_size    => $active_parameter_href->{cnv_bin_size},
                FILEHANDLE      => $XARGSFILEHANDLE,
                stdoutfile_path => $stdbasefile_path_prefix
                  . $UNDERSCORE
                  . q{partition.stdout.txt},
                stderrfile_path => $stdbasefile_path_prefix
                  . $UNDERSCORE
                  . q{partition.stderr.txt},
            }
        );
        print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

        cnvnator_calling(
            {
                infile_path     => $root_file,
                stdoutfile_path => $temp_outfile_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $DOT
                  . q{cnvnator},
                regions_ref     => [$contig],
                cnv_bin_size    => $active_parameter_href->{cnv_bin_size},
                FILEHANDLE      => $XARGSFILEHANDLE,
                stderrfile_path => $stdbasefile_path_prefix
                  . $UNDERSCORE
                  . q{calling.stderr.txt},
            }
        );
        print {$XARGSFILEHANDLE} $SEMICOLON . $SPACE;

        cnvnator_convert_to_vcf(
            {
                infile_path     => $temp_outfile_path_prefix_contig . $DOT . q{cnvnator},
                stdoutfile_path => $temp_outfile_path_prefix_contig . $outfile_suffix,
                stderrfile_path => $stdbasefile_path_prefix
                  . $UNDERSCORE
                  . q{convert_to_vcf.stderr.txt},
                referencedirectory_path => $temp_directory,
                FILEHANDLE              => $XARGSFILEHANDLE,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## Write sbatch code to supplied filehandle to concatenate variants in vcf format.
    say {$FILEHANDLE} q{## Format the VCF};

    ## Store infiles for bcftools concat
    my @concat_infile_paths;

    say {$FILEHANDLE} q{## Adding sample id instead of file prefix};
    my $rename_sample = first_value { $sample_id eq $_ }
    @{ $active_parameter_href->{sample_ids} };

    bcftools_create_reheader_samples_file(
        {
            FILEHANDLE     => $FILEHANDLE,
            sample_ids_ref => [$rename_sample],
            temp_directory => $temp_directory,
        }
    );

  CONTIG:
    foreach my $contig (@contigs) {

        ## Name intermediary files
        my $cnvnator_outfile_path =
          $temp_outfile_path_prefix . $UNDERSCORE . $contig . $outfile_suffix;
        my $fixed_vcffile_path_prefix =
          $temp_outfile_path_prefix . $UNDERSCORE . $contig . $UNDERSCORE . q{fixed};
        my $fixed_header_vcffile_path =
            $temp_outfile_path_prefix
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
                FILEHANDLE          => $FILEHANDLE,
                index               => 0,
                infile              => $cnvnator_outfile_path,
                outfile_path_prefix => $fixed_vcffile_path_prefix,
                output_type         => q{v},
                temp_directory      => $temp_directory,
                sample_ids_ref      => [$rename_sample],
            }
        );

        ## Add contigs to header
        bcftools_annotate(
            {
                FILEHANDLE      => $FILEHANDLE,
                headerfile_path => catfile( $temp_directory, q{contig_header.txt} ),
                infile_path     => $fixed_vcffile_path_prefix . $outfile_suffix,
                outfile_path    => $fixed_header_vcffile_path,
                output_type     => q{v},
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    say {$FILEHANDLE} q{## Concatenate VCFs};
    bcftools_concat(
        {
            FILEHANDLE       => $FILEHANDLE,
            infile_paths_ref => \@concat_infile_paths,
            outfile_path     => $outfile_path,
            output_type      => q{v},
            rm_dups          => 0,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => q{cnvnator},
                path             => $outfile_path,
            }
        );

        submit_recipe(
            {
                dependency_method       => q{sample_to_sample},
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                job_id_chain            => $job_id_chain,
                recipe_file_path        => $recipe_file_path,
                sample_id               => $sample_id,
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return;
}

sub _add_contigs_to_vcfheader {

## Function : Change the sample ID in the VCF header
##          : $FILEHANDLE             => Filehandle to write to
##          : $human_genome_reference => Human genome reference
##          : $temp_directory         => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $human_genome_reference;
    my $temp_directory;

    my $tmpl = {
        FILEHANDLE => {
            defined  => 1,
            required => 1,
            store    => \$FILEHANDLE,
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

    my $regexp;

    ## Execute perl
    $regexp = q?perl -nae '?;

    ## Write contig ID from header
    $regexp .= q?{print "##contig=<ID=".$F[0].",length=".$F[1].">", "\n"}'?;

    print {$FILEHANDLE} $regexp . $SPACE;

    # Reference fai file
    print {$FILEHANDLE} $human_genome_reference . $DOT . q{fai} . $SPACE;

    say {$FILEHANDLE} q{>} . $SPACE . catfile( $temp_directory, q{contig_header.txt} ),
      $NEWLINE;

    return;
}

1;
