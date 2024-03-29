package MIP::Recipes::Analysis::Mt_annotation;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile devnull };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_mt_annotation };

}

sub analysis_mt_annotation {

## Function : Annotate your mitochondrial variants here, add'l info field
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Recipe name
##          : $sample_info_href        => Info on samples and case hash {REF}

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
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::File::Path qw { remove_file_path_suffix };
    use MIP::Program::Bcftools qw{ bcftools_view };
    use MIP::Program::Gnu::Coreutils qw { gnu_cp };
    use MIP::Program::Gnu::Software::Gnu_sed qw { gnu_sed };
    use MIP::Program::HmtNote qw{ hmtnote_annotate };
    use MIP::Program::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $case_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my %infile_path        = %{ $io{in}{file_path_href} };

    my @contigs_size_ordered = @{ $file_info_href->{contigs_size_ordered} };
    my %recipe               = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $recipe{job_id_chain},
                id               => $case_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => \@contigs_size_ordered,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
            }
        )
    );

    my @outfile_paths = @{ $io{out}{file_paths} };
    my %outfile_path  = %{ $io{out}{file_path_href} };
    my %outfile_name  = %{ $io{out}{file_name_href} };
    my $outdir_path   = $io{out}{dir_path};

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
        }
    );

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    say {$filehandle} q{## Copying non-MT variant files and annotating mt vcfs};
    foreach my $contig ( keys %infile_path ) {

        if ( $contig =~ / MT|M /xsm ) {

            my $outfile_no_suffix = remove_file_path_suffix(
                {
                    file_path         => $outfile_name{$contig},
                    file_suffixes_ref => [qw{.gz}],
                }
            );
            my $temp_outfile_path =
              catfile( $active_parameter_href->{temp_directory}, $outfile_no_suffix );

            ## Hmtnote inserts whitespace in the INFO field which violates the VCF v4.2 specification
            say {$filehandle} q{## Annotate MT variants};
            hmtnote_annotate(
                {
                    filehandle   => $filehandle,
                    infile_path  => $infile_path{$contig},
                    offline      => $active_parameter_href->{mt_offline},
                    outfile_path => $temp_outfile_path,
                }
            );
            print {$filehandle} $NEWLINE;

            say {$filehandle} q{## Remove whitespace from vcf};
            my $outfile_path_no_suffix = catfile( $outdir_path, $outfile_no_suffix );
            bcftools_view(
                {
                    filehandle   => $filehandle,
                    header_only  => 1,
                    infile_path  => $temp_outfile_path,
                    outfile_path => $outfile_path_no_suffix,
                }
            );
            print {$filehandle} $NEWLINE;

            bcftools_view(
                {
                    filehandle  => $filehandle,
                    infile_path => $temp_outfile_path,
                    no_header   => 1,
                }
            );
            print {$filehandle} $PIPE . $SPACE;

            gnu_sed(
                {
                    filehandle             => $filehandle,
                    script                 => q{'s/ /_/g'},
                    stdoutfile_path_append => $outfile_path_no_suffix,
                }
            );
            print {$filehandle} $NEWLINE;

            say {$filehandle} q{## Compress and index};
            htslib_bgzip(
                {
                    filehandle  => $filehandle,
                    force       => 1,
                    infile_path => $outfile_path_no_suffix,
                }
            );
            print {$filehandle} $NEWLINE;

            htslib_tabix(
                {
                    filehandle  => $filehandle,
                    infile_path => $outfile_path{$contig},
                    preset      => q{vcf},

                }
            );
            say {$filehandle} $NEWLINE;
        }
        else {

            gnu_cp(
                {
                    filehandle   => $filehandle,
                    infile_path  => $infile_path{$contig},
                    outfile_path => $outfile_path{$contig},
                }
            );
            print {$filehandle} $NEWLINE;

            gnu_cp(
                {
                    filehandle   => $filehandle,
                    infile_path  => $infile_path{$contig} . q{.tbi},
                    outfile_path => $outfile_path{$contig} . q{.tbi},
                }
            );
            say {$filehandle} $NEWLINE;

        }
    }

    ## Close filehandle
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

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

1;
