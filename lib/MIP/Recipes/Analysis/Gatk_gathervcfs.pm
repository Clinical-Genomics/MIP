package MIP::Recipes::Analysis::Gatk_gathervcfs;

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
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_gatk_gathervcfs };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_gatk_gathervcfs {

## Function : Gather VCFs produced after gatk_genotypegvcfs done per contig.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name            => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}

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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_parameters get_recipe_attributes };
    use MIP::Gnu::Coreutils qw(gnu_mv);
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Variantcalling::Bcftools qw{ bcftools_view_and_index_vcf };
    use MIP::Program::Variantcalling::Gatk qw{ gatk_gathervcfscloud gatk_selectvariants };
    use MIP::QC::Record
      qw{ add_processing_metafile_to_sample_info add_recipe_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
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

    my $job_id_chain = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $recipe_mode        = $active_parameter_href->{$recipe_name};
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $case_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => [$infile_name_prefix],
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
                temp_directory         => $temp_directory,
            }
        )
    );
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my $outfile_path        = $io{out}{file_path};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ($recipe_file_path) = setup_script(
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

    ## Gather vcf files
    my @gatk_infile_paths =
      map { $infile_path{$_} } @{ $file_info_href->{contigs} };

    ## GATK GatherVcfsCloud
    gatk_gathervcfscloud(
        {
            FILEHANDLE           => $FILEHANDLE,
            ignore_safety_checks => 0,
            infile_paths_ref     => \@gatk_infile_paths,
            memory_allocation    => q{Xmx4G},
            outfile_path         => $outfile_path,
            temp_directory       => $temp_directory,
            verbosity            => $active_parameter_href->{gatk_logging_level},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Produce a bcf compressed and index from vcf
    if ( $active_parameter_href->{gatk_gathervcfs_bcf_file} ) {

        # Exome analysis
        if ( $consensus_analysis_type eq q{wes} ) {

            say {$FILEHANDLE} q{### Remove extra reference samples};
            say {$FILEHANDLE} q{## GATK SelectVariants};
            gatk_selectvariants(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => $outfile_path,
                    java_use_large_pages =>
                      $active_parameter_href->{java_use_large_pages},
                    memory_allocation => q{Xmx2g},
                    outfile_path      => $outfile_path_prefix
                      . $UNDERSCORE
                      . q{incnonvariantloci}
                      . $outfile_suffix,
                    referencefile_path => $referencefile_path,
                    sample_names_ref   => \@{ $active_parameter_href->{sample_ids} },
                    temp_directory     => $temp_directory,
                    verbosity          => $active_parameter_href->{gatk_logging_level},
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            ## Move to original filename
            gnu_mv(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => $outfile_path_prefix
                      . $UNDERSCORE
                      . q{incnonvariantloci}
                      . $outfile_suffix,
                    outfile_path => $outfile_path,
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }

        ## Reformat variant calling file and index
        bcftools_view_and_index_vcf(
            {
                FILEHANDLE          => $FILEHANDLE,
                infile_path         => $outfile_path,
                outfile_path_prefix => $outfile_path_prefix,
                output_type         => q{b},
            }
        );
    }

    close $FILEHANDLE;

    if ( $recipe_mode == 1 ) {

        if ( $active_parameter_href->{gatk_gathervcfs_bcf_file} == 1 ) {

            my $gbcf_file_path = $outfile_path_prefix . $DOT . q{bcf};

            add_processing_metafile_to_sample_info(
                {
                    metafile_tag     => q{gbcf_file},
                    path             => $gbcf_file_path,
                    sample_info_href => $sample_info_href,
                }
            );
        }

        ## Collect QC metadata info for later use
        add_recipe_outfile_to_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
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
    return;
}

1;
