package MIP::Recipes::Analysis::Me_annotate;

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
use MIP::Constants qw{ $EMPTY_STR $DOT $LOG_NAME $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_me_annotate };

}

sub analysis_me_annotate {

## Function : Annotate mobile elements.
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
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bcftools qw{ bcftools_annotate };
    use MIP::Program::Gnu::Coreutils qw( gnu_mv );
    use MIP::Program::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Program::Svdb qw{ svdb_query };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{ set_file_path_to_store set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Validate::Data qw{ %constraint };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

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
    my $infile_path        = $io{in}{file_path};
    my $infile_path_prefix = $io{in}{file_path_prefix};
    my $infile_suffix      = $io{in}{file_suffix};

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number = $recipe{core_number};
    my $memory      = $recipe{memory};

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $recipe{job_id_chain},
                id                     => $case_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => [$infile_name_prefix],
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );

    my $outfile_path        = $io{out}{file_path};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $case_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $memory,
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
        }
    );

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    my $svdb_infile_path = $infile_path;
    ## Default for when there are no svdb annotation files
    my $svdb_outfile_path    = $infile_path;
    my $alt_file_tag         = $UNDERSCORE . q{svdbq};
    my $query_outfile_suffix = $DOT . q{vcf};
    my $outfile_tracker      = 0;

  QUERIES:
    while ( my ( $query_db_file, $query_db_tag_info ) =
        each %{ $active_parameter_href->{me_annotate_query_files} } )
    {

        if ($outfile_tracker) {

            $svdb_infile_path =
                $outfile_path_prefix
              . $alt_file_tag
              . $query_outfile_suffix
              . $DOT
              . $outfile_tracker;
        }
        $outfile_tracker++;
        $svdb_outfile_path =
          $outfile_path_prefix . $alt_file_tag . $query_outfile_suffix . $DOT . $outfile_tracker;

        ## Get parameters
# Split query_db_tag to decide svdb input query fields
# FORMAT: filename|OUT_FREQUENCY_INFO_KEY|OUT_ALLELE_COUNT_INFO_KEY|IN_FREQUENCY_INFO_KEY|IN_ALLELE_COUNT_INFO_KEY
        my ( $query_db_tag, $out_frequency_tag_suffix, $out_allele_count_tag_suffix,
            $in_frequency_tag, $in_allele_count_tag )
          = split /[|]/sxm, $query_db_tag_info;
        my $out_frequency_tag    = $out_frequency_tag_suffix    ||= $EMPTY_STR;
        my $out_allele_count_tag = $out_allele_count_tag_suffix ||= $EMPTY_STR;

        svdb_query(
            {
                bnd_distance         => $active_parameter_href->{me_annotate_query_bnd_distance},
                dbfile_path          => $query_db_file,
                filehandle           => $filehandle,
                infile_path          => $svdb_infile_path,
                in_frequency_tag     => $in_frequency_tag,
                in_allele_count_tag  => $in_allele_count_tag,
                out_frequency_tag    => $query_db_tag . $out_frequency_tag,
                out_allele_count_tag => $query_db_tag . $out_allele_count_tag,
                stdoutfile_path      => $svdb_outfile_path,
                overlap              => $active_parameter_href->{me_annotate_query_overlap},
            }
        );
        say {$filehandle} $NEWLINE;
    }

    if ( not $constraint{is_gzipped}->($svdb_outfile_path) ) {

        htslib_bgzip(
            {
                filehandle      => $filehandle,
                infile_path     => $svdb_outfile_path,
                stdoutfile_path => $outfile_path,
                threads         => $core_number,

            }
        );
    }
    else {

        gnu_mv(
            {
                filehandle   => $filehandle,
                infile_path  => $svdb_outfile_path,
                outfile_path => $outfile_path,
            }
        );
    }
    say {$filehandle} $NEWLINE;

    htslib_tabix {
        (
            filehandle  => $filehandle,
            infile_path => $outfile_path,
        )
    };
    say {$filehandle} $NEWLINE;

    ## Close filehandle
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_case},
                log                               => $log,
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
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
