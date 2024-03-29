package MIP::Recipes::Analysis::Retroseq;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $DOT $LOG_NAME $NEWLINE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_retroseq };

}

sub analysis_retroseq {

## Function : Discover and call mobile elememnts using RetroSeq.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Recipe name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

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
        sample_id => {
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bcftools qw{ bcftools_sort };
    use MIP::Program::Gnu::Coreutils qw{ gnu_echo };
    use MIP::Program::Htslib qw{ htslib_tabix };
    use MIP::Program::Picardtools qw{ picardtools_updatevcfsequencedictionary };
    use MIP::Program::Retroseq qw{ retroseq_call retroseq_discover };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

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
            temp_directory => $temp_directory,
        }
    );

    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path_prefix = $io{in}{file_path_prefix};
    my $infile_suffix      = $io{in}{file_suffix};
    my $infile_path        = $infile_path_prefix . $infile_suffix;

    my $analysis_type = $active_parameter_href->{analysis_type}{$sample_id};

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $recipe{job_id_chain},
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => [$infile_name_prefix],
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
                temp_directory         => $temp_directory,
            }
        )
    );

    my $outfile_name_prefix = $io{out}{file_name_prefix};
    my $outdir_path         = $io{out}{dir_path};
    my $outfile_path        = $io{out}{file_path};
    my $outfile_path_prefix = $io{out}{file_path_prefix};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe{core_number},
            directory_id          => $sample_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe{memory},
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
        }
    );

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    say {$filehandle} q{## Prepare mobile element references};

    my @mobile_element_references;
    while ( my ( $me_reference_path, $me_type ) =
        each %{ $active_parameter_href->{mobile_element_reference} } )
    {

        push @mobile_element_references, $me_type . q{\t} . $me_reference_path . q{\n};
    }

    my $mobile_element_reference_file = catfile( $outdir_path, q{mobile_element_references.tsv} );
    gnu_echo(
        {
            enable_interpretation => 1,
            filehandle            => $filehandle,
            no_trailing_newline   => 1,
            outfile_path          => $mobile_element_reference_file,
            strings_ref           => \@mobile_element_references,
        }
    );
    say ${filehandle} $NEWLINE;

    say {$filehandle} q{## Discover mobile elements};
    my $discover_outfile_path = $outfile_path_prefix . q{.bed};
    retroseq_discover(
        {
            filehandle              => $filehandle,
            infile_path             => $infile_path,
            mobile_element_tsv_path => $mobile_element_reference_file,
            outfile_path            => $discover_outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Call mobile elements};
    my $call_outfile_path = $outfile_path_prefix . q{_call.vcf};
    retroseq_call(
        {
            filehandle           => $filehandle,
            infile_path          => $infile_path,
            retroseq_bed_path    => $discover_outfile_path,
            outfile_path         => $call_outfile_path,
            reference_fasta_path => $active_parameter_href->{human_genome_reference},
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Cleaning up the vcf};
    my $picardtools_outfile_path = $outfile_path_prefix . q{_header_fix.vcf};
    picardtools_updatevcfsequencedictionary(
        {
            filehandle   => $filehandle,
            infile_path  => $call_outfile_path,
            java_jar     => catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
            outfile_path => $picardtools_outfile_path,
            sequence_dictionary => $infile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    my $sort_memory = $recipe{memory} - 2;
    bcftools_sort(
        {
            filehandle     => $filehandle,
            infile_path    => $picardtools_outfile_path,
            max_mem        => $sort_memory . q{G},
            outfile_path   => $outfile_path,
            output_type    => q{z},
            temp_directory => catdir( $temp_directory, q{bcftools_sort} ),
        }
    );
    say {$filehandle} $NEWLINE;

    htslib_tabix(
        {
            filehandle  => $filehandle,
            infile_path => $outfile_path,
            preset      => q{vcf},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Close filehandles
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => $outfile_name_prefix,
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_sample},
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
                log                               => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_id          => $sample_id,
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

1;
