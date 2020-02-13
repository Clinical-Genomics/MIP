package MIP::Recipes::Analysis::Gatk_cnnscorevariants;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw { uniq };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $ASTERISK $DASH $DOT $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_gatk_cnnscorevariants };

}

sub analysis_gatk_cnnscorevariants {

## Function : GATK CNNScoreVariants analysis recipe for single sample calling
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Pedigree qw{ create_fam_file gatk_pedigree_flag };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bcftools qw{ bcftools_norm };
    use MIP::Program::Gatk qw{ gatk_cnnscorevariants };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

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
    my $infile_path        = $io{in}{file_path};

    my $job_id_chain = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode        = $active_parameter_href->{$recipe_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my %recipe_resource    = get_recipe_resources(
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
    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my $outfile_path        = $io{out}{file_path};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $recipe_resource{core_number},
            directory_id                    => $case_id,
            filehandle                      => $filehandle,
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

    ## Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) =
      splitpath( $recipe_info_path . $DOT . q{stderr.txt} );

    ## Collect BAM infiles for dependence recipes streams for all sample_ids
    my %recipe_tag_keys = ( gatk_baserecalibration => q{out}, );

    ## Store sample id bam infiles
    my @bam_infiles_paths;
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
            my $infile_path_prefix_bam = $sample_io{$stream}{file_path_prefix};
            my $infile_suffix_bam      = $sample_io{$stream}{file_suffix};
            my $infile_path_bam        = $infile_path_prefix_bam . $infile_suffix_bam;
            push @bam_infiles_paths, $infile_path_bam;
        }
    }
    say {$filehandle} q{wait}, $NEWLINE;

    ## Create .fam file to be used in variant calling analyses
    my $fam_file_path = catfile( $outdir_path_prefix, $case_id . $DOT . q{fam} );
    create_fam_file(
        {
            active_parameter_href => $active_parameter_href,
            execution_mode        => q{system},
            fam_file_path         => $fam_file_path,
            filehandle            => $filehandle,
            log                   => $log,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Check if "--pedigree" and "--pedigreeValidationType" should be included in analysis
    my %commands = gatk_pedigree_flag(
        {
            fam_file_path => $fam_file_path,
        }
    );

    ## GATK CNNScoreVariants
    say {$filehandle} q{## GATK CNNScoreVariants};

    my $cnn_outfile_path = $outfile_path_prefix . $UNDERSCORE . q{cnn} . $outfile_suffix;
    my $mv_infile_path   = $cnn_outfile_path;
    gatk_cnnscorevariants(
        {
            alignment_infile_paths_ref => \@bam_infiles_paths,
            filehandle                 => $filehandle,
            infile_path                => $infile_path,
            outfile_path               => $cnn_outfile_path,
            referencefile_path         => $referencefile_path,
            temp_directory             => $temp_directory,
            verbosity                  => $active_parameter_href->{gatk_logging_level},
        }
    );
    say {$filehandle} $NEWLINE;

    if ( not $active_parameter_href->{gatk_variantrecalibration_keep_unnormalised} ) {

        ## Bcftools norm, left-align and normalize indels, split multiallelics
        my $norm_outfile_path =
          $outfile_path_prefix . $UNDERSCORE . q{normalized} . $outfile_suffix;
        $mv_infile_path = $norm_outfile_path;
        bcftools_norm(
            {
                filehandle      => $filehandle,
                infile_path     => $cnn_outfile_path,
                multiallelic    => $DASH,
                outfile_path    => $norm_outfile_path,
                output_type     => q{v},
                reference_path  => $referencefile_path,
                stderrfile_path => $outfile_path_prefix
                  . $UNDERSCORE
                  . q{normalized.stderr},
            }
        );
        say {$filehandle} $NEWLINE;
    }

    ## Change name of file to accomodate downstream
    gnu_mv(
        {
            filehandle   => $filehandle,
            infile_path  => $mv_infile_path,
            outfile_path => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    close $filehandle;

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        # Used to find order of samples in qccollect downstream
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => q{pedigree_check},
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command            => $profile_base_command,
                case_id                 => $case_id,
                dependency_method       => q{sample_to_case},
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_chain            => $job_id_chain,
                job_id_href             => $job_id_href,
                job_reservation_name    => $active_parameter_href->{job_reservation_name},
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
