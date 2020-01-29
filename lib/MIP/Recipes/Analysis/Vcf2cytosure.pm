package MIP::Recipes::Analysis::Vcf2cytosure;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
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
use MIP::Constants
  qw{ $AMPERSAND $ASTERISK $DOT $LOG_NAME $NEWLINE $SINGLE_QUOTE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.20;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_vcf2cytosure };

}

sub analysis_vcf2cytosure {

## Function : Convert VCF with structural variations to the “.CGH” format used by the CytoSure Interpret Software
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $bin_size                => Bin size
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
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
    my $bin_size;
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
        bin_size => {
            default     => $arg_href->{active_parameter_href}{tiddit_bin_size},
            strict_type => 1,
            store       => \$bin_size
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

    use MIP::Cluster qw{ get_core_number update_memory_allocation };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Program::Vcf2cytosure qw{ vcf2cytosure_convert };
    use MIP::Processmanagement::Processes qw{ print_wait submit_recipe };
    use MIP::Program::Bcftools qw{ bcftools_view };
    use MIP::Program::Tiddit qw{ tiddit_coverage };
    use MIP::Sample_info
      qw{ get_pedigree_sample_id_attributes set_file_path_to_store set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my $job_id_chain = get_recipe_attributes(
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

    ## Set and get the io files per chain, id and stream
    my %io = parse_io_outfiles(
        {
            chain_id         => $job_id_chain,
            id               => $case_id,
            file_info_href   => $file_info_href,
            outdata_dir      => $active_parameter_href->{outdata_dir},
            file_name_prefix => $case_id,
            iterators_ref    => $active_parameter_href->{sample_ids},
            parameter_href   => $parameter_href,
            recipe_name      => $recipe_name,
        }
    );

    my %outfile_name        = %{ $io{out}{file_name_href} };
    my %outfile_path        = %{ $io{out}{file_path_href} };
    my $outfile_path_prefix = $io{out}{file_path_prefix};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    my $core_number = get_core_number(
        {
            max_cores_per_node   => $active_parameter_href->{max_cores_per_node},
            modifier_core_number => scalar @{ $active_parameter_href->{sample_ids} },
            recipe_core_number   => $recipe_resource{core_number},
        }
    );

    ## Update recipe memory allocation
    my $memory_allocation = update_memory_allocation(
        {
            node_ram_memory           => $active_parameter_href->{node_ram_memory},
            parallel_processes        => $core_number,
            process_memory_allocation => $recipe_resource{memory},
        }
    );

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $case_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            log                             => $log,
            memory_allocation               => $memory_allocation,
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
            temp_directory                  => $temp_directory,
        }
    );

    ## Store file info from ".bam", ".tab" and ".vcf"
    my %vcf2cytosure_file_info;

    ### Get case vcf from sv_anno
    ## Get the io infiles per chain and id
    my %case_io = get_io_files(
        {
            id             => $case_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => q{sv_annotate},
            stream         => q{out},
        }
    );
    my $infile_path_prefix = $case_io{out}{file_path_prefix};
    my $infile_suffix      = $case_io{out}{file_suffix};
    my $infile_path        = $infile_path_prefix . $infile_suffix;
    $vcf2cytosure_file_info{$case_id}{in}{$infile_suffix} =
      $infile_path;

    ## Collect BAM infiles for dependence recipes streams for all sample_ids
    my %recipe_tag_keys = ( gatk_baserecalibration => q{out}, );

    my $process_batches_count = 1;
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
                }
            );
            my $infile_path_prefix_bam = $sample_io{$stream}{file_path_prefix};
            my $infile_suffix_bam      = $sample_io{$stream}{file_suffix};
            my $infile_path_bam        = $infile_path_prefix_bam . $infile_suffix_bam;

            $vcf2cytosure_file_info{$sample_id}{in}{$infile_suffix_bam} =
              $infile_path_bam;
        }
    }

    ### SHELL:

    say {$filehandle} q{## Creating coverage file with tiddit -cov for samples};

  SAMPLE_ID:
    while ( my ( $sample_id_index, $sample_id ) =
        each @{ $active_parameter_href->{sample_ids} } )
    {

        my $tiddit_cov_file_path =
          $outfile_path_prefix . $UNDERSCORE . q{tiddit} . $UNDERSCORE . $sample_id;

        ## Store file for use downstream
        $vcf2cytosure_file_info{$sample_id}{in}{q{.tab}} =
          $tiddit_cov_file_path . q{.tab};

        ## Tiddit coverage
        tiddit_coverage(
            {
                bin_size            => $bin_size,
                filehandle          => $filehandle,
                infile_path         => $vcf2cytosure_file_info{$sample_id}{in}{q{.bam}},
                outfile_path_prefix => $tiddit_cov_file_path,
            }
        );
        say {$filehandle} $AMPERSAND . $SPACE . $NEWLINE;
    }
    say {$filehandle} q{wait}, $NEWLINE;

    # Extract SV from this sample from merged SV VCF file
    say {$filehandle} q{## Using bcftools_view to extract SVs for samples} . $NEWLINE;

  SAMPLE_ID:
    while ( my ( $sample_id_index, $sample_id ) =
        each @{ $active_parameter_href->{sample_ids} } )
    {
        my $bcftools_outfile_path =
            $outfile_path_prefix
          . $UNDERSCORE
          . q{filtered}
          . $UNDERSCORE
          . $sample_id . q{.vcf};

        ## Store file for use downstream
        $vcf2cytosure_file_info{$sample_id}{in}{q{.vcf}} =
          $bcftools_outfile_path;

        # Bcftools view
        bcftools_view(
            {
                exclude      => $active_parameter_href->{vcf2cytosure_exclude_filter},
                filehandle   => $filehandle,
                infile_path  => $vcf2cytosure_file_info{$case_id}{in}{q{.vcf}},
                samples_ref  => [$sample_id],
                outfile_path => $bcftools_outfile_path,
            }
        );
        say {$filehandle} $AMPERSAND . $SPACE . $NEWLINE;
    }
    say {$filehandle} q{wait}, $NEWLINE;

    say {$filehandle}
      q{## Converting sample's SV VCF file into cytosure, using Vcf2cytosure} . $NEWLINE;
  SAMPLE_ID:
    while ( my ( $sample_id_index, $sample_id ) =
        each @{ $active_parameter_href->{sample_ids} } )
    {

        # Get parameter
        my $sample_id_sex = get_pedigree_sample_id_attributes(
            {
                attribute        => q{sex},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        ## Special case for unknown sex
        if ( $sample_id_sex eq q{unknown} ) {

            $sample_id_sex = undef;
        }

        vcf2cytosure_convert(
            {
                blacklist_file_path => $active_parameter_href->{vcf2cytosure_blacklist},
                coverage_file       => $vcf2cytosure_file_info{$sample_id}{in}{q{.tab}},
                filehandle          => $filehandle,
                maxbnd              => $active_parameter_href->{vcf2cytosure_maxbnd},
                outfile_path        => $outfile_path{$sample_id},
                sex                 => $sample_id_sex,
                variant_size        => $active_parameter_href->{vcf2cytosure_var_size},
                vcf_infile_path     => $vcf2cytosure_file_info{$sample_id}{in}{q{.vcf}},
            }
        );
        say {$filehandle} $AMPERSAND . $SPACE . $NEWLINE;

        if ( $recipe_mode == 1 ) {

            set_recipe_outfile_in_sample_info(
                {
                    infile           => $outfile_name{$sample_id},
                    path             => $outfile_path{$sample_id},
                    recipe_name      => q{vcf2cytosure},
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );

            set_file_path_to_store(
                {
                    format           => q{meta},
                    id               => $sample_id,
                    path             => $outfile_path{$sample_id},
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                }
            );
        }
    }
    say {$filehandle} q{wait}, $NEWLINE;

    if ( $recipe_mode == 1 ) {

        submit_recipe(
            {
                base_command            => $profile_base_command,
                case_id                 => $case_id,
                dependency_method       => q{case_to_island},
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
