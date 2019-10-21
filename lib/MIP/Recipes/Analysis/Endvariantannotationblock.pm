package MIP::Recipes::Analysis::Endvariantannotationblock;

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

## MIPs lib/
use MIP::Constants qw{ $ASTERISK $DOT $EMPTY_STR $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_endvariantannotationblock };

}

sub analysis_endvariantannotationblock {

## Function : Concatenate ouput from variant annotation block
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $file_path               => File path
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
    my $file_path;
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
        file_path      => { store => \$file_path, strict_type => 1, },
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
        xargs_file_counter => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::Analysis qw{ get_vcf_parser_analysis_suffix };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Gnu::Software::Gnu_grep qw{ gnu_grep };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Utility::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Program::Variantcalling::Gatk qw{ gatk_concatenate_variants };
    use MIP::Sample_info qw{ set_most_complete_vcf set_recipe_metafile_in_sample_info };
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
            temp_directory => $temp_directory,
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path_prefix = $io{in}{file_path_prefix};

    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my @contigs                 = @{ $file_info_href->{contigs} };
    my $job_id_chain            = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    my @vcfparser_analysis_types = get_vcf_parser_analysis_suffix(
        {
            vcfparser_outfile_count =>
              $active_parameter_href->{sv_vcfparser_outfile_count},
        }
    );

    ## Set and get the io files per chain, id and stream
    my @set_outfile_name_prefixes =
      map { $infile_name_prefix . $_ } @vcfparser_analysis_types;
    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $case_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => \@vcfparser_analysis_types,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
                temp_directory   => $temp_directory,
            }
        )
    );

    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my @outfile_suffixes    = @{ $io{out}{file_suffixes} };
    my @outfile_paths       = @{ $io{out}{file_paths} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle      = IO::Handle->new();
    my $xargsfilehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $recipe_resource{core_number},
            filehandle                      => $filehandle,
            directory_id                    => $case_id,
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

  ANALYSIS_SUFFIXES:
    while ( my ( $analysis_suffix_index, $analysis_suffix ) = each @outfile_suffixes ) {

        my @concat_contigs = @contigs;
        my $infile_postfix = q{.vcf};
        my $metafile_tag   = q{research};

        ## Update contigs list using select file contigs
        if ( $analysis_suffix eq q{.selected.vcf} ) {

            @concat_contigs = @{ $file_info_href->{select_file_contigs} };
            $infile_postfix = $UNDERSCORE . q{selected.vcf};
            $metafile_tag   = q{clinical};
        }

        ## Writes sbatch code to supplied filehandle to concatenate variants in vcf format. Each array element is combined with the infile prefix and postfix.
        gatk_concatenate_variants(
            {
                active_parameter_href => $active_parameter_href,
                elements_ref          => \@concat_contigs,
                filehandle            => $filehandle,
                infile_prefix         => $infile_path_prefix,
                infile_postfix        => $infile_postfix,
                outfile_path_prefix   => $outfile_path_prefix,
                outfile_suffix        => $analysis_suffix,
            }
        );

        ## Remove variants in hgnc_id list from vcf
        if ( $active_parameter_href->{endvariantannotationblock_remove_genes_file} ) {

            my $grep_outfile_path =
              $outfile_path_prefix . $UNDERSCORE . q{filtered} . $analysis_suffix;
            ## Removes contig_names from contigs array if no male or other found
            gnu_grep(
                {
                    filehandle       => $filehandle,
                    filter_file_path => catfile(
                        $reference_dir,
                        $active_parameter_href->{sv_reformat_remove_genes_file}
                    ),
                    infile_path     => $outfile_paths[$analysis_suffix_index],
                    stdoutfile_path => $grep_outfile_path,
                    invert_match    => 1,
                }
            );
            say {$filehandle} $NEWLINE;

            ## Save filtered file
            $sample_info_href->{recipe}{$recipe_name}
              {reformat_remove_genes_file}{$metafile_tag}{path} = $grep_outfile_path;
        }

        if ( $active_parameter_href->{rankvariant_binary_file} ) {

            my $bgzip_outfile_path =
              $outfile_paths[$analysis_suffix_index] . $DOT . q{gz};
            ## Compress or decompress original file or stream to outfile (if supplied)
            htslib_bgzip(
                {
                    filehandle      => $filehandle,
                    infile_path     => $outfile_paths[$analysis_suffix_index],
                    stdoutfile_path => $bgzip_outfile_path,
                    write_to_stdout => 1,
                }
            );
            say {$filehandle} $NEWLINE;

            ## Index file using tabix
            htslib_tabix(
                {
                    filehandle  => $filehandle,
                    force       => 1,
                    infile_path => $bgzip_outfile_path,
                    preset      => substr $outfile_suffix,
                    1,
                }
            );
            say {$filehandle} $NEWLINE;
        }

        ## Adds the most complete vcf file to sample_info
        set_most_complete_vcf(
            {
                active_parameter_href     => $active_parameter_href,
                path                      => $outfile_paths[$analysis_suffix_index],
                recipe_name               => $recipe_name,
                sample_info_href          => $sample_info_href,
                vcfparser_outfile_counter => $analysis_suffix_index,
            }
        );

        if ( $recipe_mode == 1 ) {

            set_recipe_metafile_in_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    recipe_name      => $recipe_name,
                    metafile_tag     => $metafile_tag,
                    path             => $outfile_paths[$analysis_suffix_index],
                }
            );

            if ( $active_parameter_href->{rankvariant_binary_file} ) {

                $sample_info_href->{vcf_binary_file}{$metafile_tag}{path} =
                  $outfile_paths[$analysis_suffix_index] . $DOT . q{gz};
            }
        }
    }

    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe_mode == 1 ) {

        submit_recipe(
            {
                base_command            => $profile_base_command,
                dependency_method       => q{sample_to_case},
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_chain            => $job_id_chain,
                job_id_href             => $job_id_href,
                log                     => $log,
                recipe_file_path        => $recipe_file_path,
                sample_ids_ref          => \@{ $active_parameter_href->{sample_ids} },
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

1;
