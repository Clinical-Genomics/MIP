package MIP::Recipes::Analysis::Expansionhunter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw{ each_array };
use MIP::Constants qw{ $AMPERSAND $ASTERISK $DOT $LOG_NAME $NEWLINE $UNDERSCORE };
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_expansionhunter };

}

sub analysis_expansionhunter {

## Function : Call expansions of Short Tandem Repeats (STR) using Expansion Hunter
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
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
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

    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ print_wait submit_recipe };
    use MIP::Program::Bcftools
      qw{ bcftools_index bcftools_norm bcftools_rename_vcf_samples bcftools_view };
    use MIP::Program::Expansionhunter qw{ expansionhunter };
    use MIP::Program::Htslib qw{ htslib_bgzip };
    use MIP::Program::Samtools qw{ samtools_index samtools_sort };
    use MIP::Program::Stranger qw{ stranger };
    use MIP::Program::Svdb qw{ svdb_merge };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info
      qw{ get_pedigree_sample_id_attributes set_file_path_to_store set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    Readonly my $SCALING_FACTOR => 0.75;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my $human_genome_reference = $arg_href->{active_parameter_href}{human_genome_reference};
    my $variant_catalog_file_path =
      $active_parameter_href->{expansionhunter_variant_catalog_file_path};
    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set and get the io files per chain, id and stream
    my %io = parse_io_outfiles(
        {
            chain_id               => $recipe{job_id_chain},
            id                     => $case_id,
            file_info_href         => $file_info_href,
            file_name_prefixes_ref => [$case_id],
            outdata_dir            => $active_parameter_href->{outdata_dir},
            parameter_href         => $parameter_href,
            recipe_name            => $recipe_name,
        }
    );

    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_constant_suffix};
    my $outfile_path        = $outfile_path_prefix . $outfile_suffix;

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe{core_number},
            directory_id          => $case_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe{memory},
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            temp_directory        => $temp_directory,
        }
    );

    ### SHELL:

    my %exphun_sample_file;
    my @decompose_outfile_paths;

    ## Collect infiles for all sample_ids to enable migration to temporary directory
  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Get the io infiles per chain and id
        my %sample_io = get_io_files(
            {
                id             => $sample_id,
                file_info_href => $file_info_href,
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
                stream         => q{in},
                temp_directory => $temp_directory,
            }
        );
        my $infile_path_prefix = $sample_io{in}{file_path_prefix};
        my $infile_suffix      = $sample_io{in}{file_suffix};
        my $infile_path        = $infile_path_prefix . $infile_suffix;

        say {$filehandle} q{## Run ExpansionHunter};

        # Get parameter
        my $sample_id_sex = get_pedigree_sample_id_attributes(
            {
                attribute        => q{sex},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        ## Special case for unknown sex
        $sample_id_sex = _update_sample_id_sex( { sample_id_sex => $sample_id_sex, } );

        my $sample_outfile_path_prefix = $outfile_path_prefix . $UNDERSCORE . $sample_id;
        expansionhunter(
            {
                filehandle                => $filehandle,
                infile_path               => $infile_path,
                outfile_path_prefix       => $sample_outfile_path_prefix,
                reference_genome_path     => $human_genome_reference,
                sex                       => $sample_id_sex,
                threads                   => $recipe{core_number},
                variant_catalog_file_path => $variant_catalog_file_path,
            }
        );
        say {$filehandle} $NEWLINE;

        my $exphun_sample_vcf = $sample_outfile_path_prefix . $outfile_suffix;
        my $exphun_sample_bam = $sample_outfile_path_prefix . q{.bam};
        $exphun_sample_file{$sample_id}{out_vcf} = $exphun_sample_vcf;
        $exphun_sample_file{$sample_id}{out_bam} = $exphun_sample_bam;

        my $thread_memory = int( $recipe{memory} / $recipe{core_number} );
        $thread_memory = int( $thread_memory * $SCALING_FACTOR );
        samtools_sort(
            {
                filehandle            => $filehandle,
                infile_path           => $sample_outfile_path_prefix . q{_realigned.bam},
                max_memory_per_thread => $thread_memory . q{G},
                outfile_path          => $exphun_sample_bam,
                output_format         => q{bam},
                temp_file_path_prefix => $temp_directory,
                thread_number         => $recipe{core_number},
            }
        );
        say {$filehandle} $NEWLINE;

        samtools_index(
            {
                bai_format  => 1,
                filehandle  => $filehandle,
                infile_path => $exphun_sample_bam,
            }
        );
        say {$filehandle} $NEWLINE;

        htslib_bgzip(
            {
                filehandle      => $filehandle,
                infile_path     => $exphun_sample_vcf,
                stdoutfile_path => $exphun_sample_vcf . q{.gz},
                threads         => $recipe{core_number},
                write_to_stdout => 1,
            }
        );
        say {$filehandle} $NEWLINE;

        bcftools_index(
            {
                filehandle  => $filehandle,
                infile_path => $exphun_sample_vcf . q{.gz},
                output_type => q{tbi},
            }
        );
        say {$filehandle} $NEWLINE;

        say {$filehandle} q{## Split multiallelic variants};
        my $norm_outfile_path =
            $outfile_path_prefix
          . $UNDERSCORE
          . q{decompose}
          . $UNDERSCORE
          . $sample_id
          . $outfile_suffix;
        bcftools_norm(
            {
                filehandle   => $filehandle,
                infile_path  => $exphun_sample_vcf . q{.gz},
                multiallelic => q{-},
                outfile_path => $norm_outfile_path,
                threads      => $recipe{core_number},
            }
        );
        say {$filehandle} $NEWLINE;

        push @decompose_outfile_paths, $norm_outfile_path;
    }

    say {$filehandle} q{## Merge calls};
    my $svdb_outfile_path =
      $outfile_path_prefix . $UNDERSCORE . q{decompose_svdbmerge} . $outfile_suffix;
    svdb_merge(
        {
            filehandle       => $filehandle,
            infile_paths_ref => \@decompose_outfile_paths,
            notag            => 1,
            stdoutfile_path  => $svdb_outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Stranger annotation};
    my $stranger_outfile_path =
      $outfile_path_prefix . $UNDERSCORE . q{decompose_svdbmerge_ann} . $outfile_suffix;
    stranger(
        {
            case_id           => $case_id,
            filehandle        => $filehandle,
            infile_path       => $svdb_outfile_path,
            repeats_file_path => $variant_catalog_file_path,
            stdoutfile_path   => $stranger_outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Adding sample id instead of file prefix};
    bcftools_rename_vcf_samples(
        {
            filehandle          => $filehandle,
            index               => 1,
            index_type          => q{csi},
            infile              => $stranger_outfile_path,
            outfile_path_prefix => $outfile_path_prefix,
            output_type         => q{z},
            temp_directory      => $temp_directory,
            sample_ids_ref      => \@{ $active_parameter_href->{sample_ids} },
        }
    );

    close $filehandle;

    if ( $recipe{mode} == 1 ) {

        set_recipe_outfile_in_sample_info(
            {
                sample_info_href => $sample_info_href,
                recipe_name      => $recipe_name,
                path             => $outfile_path . $DOT . q{gz},
            }
        );

        set_file_path_to_store(
            {
                format           => q{vcf},
                id               => $case_id,
                path             => $outfile_path . $DOT . q{gz},
                path_index       => $outfile_path . $DOT . q{gz} . $DOT . q{csi},
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
                tag              => q{sv} . $UNDERSCORE . q{str},
            }
        );

        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            set_file_path_to_store(
                {
                    format           => q{vcf},
                    id               => $sample_id,
                    path             => $exphun_sample_file{$sample_id}{out_vcf},
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                    tag              => q{str_variants},
                }
            );
            set_file_path_to_store(
                {
                    format           => q{bam},
                    id               => $sample_id,
                    path             => $exphun_sample_file{$sample_id}{out_bam},
                    path_index       => $exphun_sample_file{$sample_id}{out_bam} . q{.bai},
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                    tag              => q{str_alignment},
                }
            );
        }

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_case},
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
                log                               => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_ids_ref     => \@{ $active_parameter_href->{sample_ids} },
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub _update_sample_id_sex {

    ## Function : If sample sex is unknown return undef otherwise return $sample_id_sex
    ## Returns  : $sample_id_sex or undef
    ## Arguments: $sample_id_sex => Sex of sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_id_sex;

    my $tmpl = {
        sample_id_sex => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id_sex,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( $sample_id_sex eq q{unknown} );

    return $sample_id_sex;
}

1;
