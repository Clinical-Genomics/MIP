package MIP::Recipes::Analysis::Vcf2cytosure;

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
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_vcf2cytosure };

}

## Constants
Readonly my $AMPERSAND    => q{&};
Readonly my $ASTERISK     => q{*};
Readonly my $DOT          => q{.};
Readonly my $NEWLINE      => qq{\n};
Readonly my $SINGLE_QUOTE => q{'};
Readonly my $SPACE        => q{ };
Readonly my $SV_LENGTH    => 3000;
Readonly my $UNDERSCORE   => q{_};

sub analysis_vcf2cytosure {

## Function : Convert VCF with structural variations to the “.CGH” format used by the CytoSure Interpret Software
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $bin_size                => Bin size
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $sample_info_href;

    ## Default(s)
    my $bin_size;
    my $family_id;
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

    use MIP::Cluster qw{ get_core_number };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter
      qw{ get_module_parameters get_program_attributes get_program_parameters };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Program::Variantcalling::Vcf2cytosure qw{ vcf2cytosure_convert };
    use MIP::Processmanagement::Processes qw{ print_wait };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_family_dead_end };
    use MIP::Program::Variantcalling::Bcftools qw{ bcftools_view };
    use MIP::Program::Variantcalling::Tiddit qw{ tiddit_coverage };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    my $job_id_chain = get_program_attributes(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
            attribute      => q{chain},
        }
    );
    my $program_mode = $active_parameter_href->{$program_name};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Set and get the io files per chain, id and stream
    my %io = parse_io_outfiles(
        {
            chain_id         => $job_id_chain,
            id               => $family_id,
            file_info_href   => $file_info_href,
            outdata_dir      => $active_parameter_href->{outdata_dir},
            file_name_prefix => $family_id,
            iterators_ref    => $active_parameter_href->{sample_ids},
            parameter_href   => $parameter_href,
            program_name     => $program_name,
            temp_directory   => $temp_directory,
        }
    );

    my %outfile_name             = %{ $io{out}{file_name_href} };
    my %outfile_path             = %{ $io{out}{file_path_href} };
    my $temp_outfile_path_prefix = $io{temp}{file_path_prefix};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
            modifier_core_number =>
              scalar @{ $active_parameter_href->{sample_ids} },
            module_core_number =>
              $active_parameter_href->{module_core_number}{$program_name},
        }
    );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $family_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            program_directory               => $program_name,
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    ## Store file info from ".bam", ".tab" and ".vcf"
    my %vcf2cytosure_file_info;

    ### Get family vcf from sv_anno
    ## Get the io infiles per chain and id
    my %family_io = get_io_files(
        {
            id             => $family_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            program_name   => q{sv_annotate},
            stream         => q{out},
            temp_directory => $temp_directory,
        }
    );
    my $infile_path_prefix = $family_io{out}{file_path_prefix};
    my $infile_suffix      = $family_io{out}{file_suffix};
    my $infile_path =
      $infile_path_prefix . substr( $infile_suffix, 0, 2 ) . $ASTERISK;
    my $temp_infile_path_prefix = $family_io{temp}{file_path_prefix};
    my $temp_infile_path        = $temp_infile_path_prefix . $infile_suffix;
    $vcf2cytosure_file_info{$family_id}{in}{$infile_suffix} =
      $temp_infile_path;

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $infile_path,
            outfile_path => $temp_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Collect BAM infiles for dependence programs streams for all sample_ids
    my %program_tag_keys = ( gatk_baserecalibration => q{out}, );

    my $process_batches_count = 1;
    while ( my ( $sample_id_index, $sample_id ) =
        each @{ $active_parameter_href->{sample_ids} } )
    {

      PROGRAM_TAG:
        while ( my ( $program_tag, $stream ) = each %program_tag_keys ) {

            ## Get the io infiles per chain and id
            my %sample_io = get_io_files(
                {
                    id             => $sample_id,
                    file_info_href => $file_info_href,
                    parameter_href => $parameter_href,
                    program_name   => $program_tag,
                    stream         => $stream,
                    temp_directory => $temp_directory,
                }
            );
            my $infile_path_prefix_bam = $sample_io{$stream}{file_path_prefix};
            my $infile_suffix_bam      = $sample_io{$stream}{file_suffix};
            my $infile_path_bam =
                $infile_path_prefix_bam
              . substr( $infile_suffix_bam, 0, 2 )
              . $ASTERISK;
            my $temp_infile_path_prefix_bam =
              $sample_io{temp}{file_path_prefix};
            my $temp_infile_path_bam =
              $temp_infile_path_prefix_bam . $infile_suffix_bam;

            $vcf2cytosure_file_info{$sample_id}{in}{$infile_suffix_bam} =
              $temp_infile_path_bam;

            $process_batches_count = print_wait(
                {
                    FILEHANDLE            => $FILEHANDLE,
                    max_process_number    => $core_number,
                    process_batches_count => $process_batches_count,
                    process_counter       => $sample_id_index,
                }
            );

            ## Copy file(s) to temporary directory
            say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
            migrate_file(
                {
                    FILEHANDLE   => $FILEHANDLE,
                    infile_path  => $infile_path_bam,
                    outfile_path => $temp_directory,
                }
            );
        }
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Excute vcf2cytosure just to get an error message for version
    say {$FILEHANDLE} q{## Log vcf2cytosure version - use dummy parameters}
      . $NEWLINE;
    my $stderrfile_path = $program_info_path . $DOT . q{stderr.txt};
    vcf2cytosure_convert(
        {
            coverage_file   => q{Na},
            FILEHANDLE      => $FILEHANDLE,
            stderrfile_path => $stderrfile_path,
            vcf_infile_path => q{Na},
            version         => 1,
        }
    );

    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Creating coverage file with tiddit -cov for samples};

  SAMPLE_ID:
    while ( my ( $sample_id_index, $sample_id ) =
        each @{ $active_parameter_href->{sample_ids} } )
    {

        my $tiddit_temp_cov_file_path =
            $temp_outfile_path_prefix
          . $UNDERSCORE
          . q{tiddit}
          . $UNDERSCORE
          . $sample_id;

        ## Store file for use downstream
        $vcf2cytosure_file_info{$sample_id}{in}{q{.tab}} =
          $tiddit_temp_cov_file_path . q{.tab};

        ## Tiddit coverage
        tiddit_coverage(
            {
                bin_size    => $bin_size,
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $vcf2cytosure_file_info{$sample_id}{in}{q{.bam}},
                outfile_path_prefix => $tiddit_temp_cov_file_path,
            }
        );
        say {$FILEHANDLE} $AMPERSAND . $SPACE . $NEWLINE;
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    # Extract SV from this sample from merged SV VCF file
    say {$FILEHANDLE} q{## Using bcftools_view to extract SVs for samples}
      . $NEWLINE;

  SAMPLE_ID:
    while ( my ( $sample_id_index, $sample_id ) =
        each @{ $active_parameter_href->{sample_ids} } )
    {
        my $bcftools_temp_outfile_path =
            $temp_outfile_path_prefix
          . $UNDERSCORE
          . q{filtered}
          . $UNDERSCORE
          . $sample_id . q{.vcf};
        ## Store file for use downstream
        $vcf2cytosure_file_info{$sample_id}{in}{q{.vcf}} =
          $bcftools_temp_outfile_path;

        # Bcftools view
        bcftools_view(
            {
                exclude =>
                  $active_parameter_href->{vcf2cytosure_exclude_filter},
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $vcf2cytosure_file_info{$family_id}{in}{q{.vcf}},
                samples_ref => [$sample_id],
                outfile_path => $bcftools_temp_outfile_path,
            }
        );
        say {$FILEHANDLE} $AMPERSAND . $SPACE . $NEWLINE;
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    say {$FILEHANDLE}
      q{## Converting sample's SV VCF file into cytosure, using Vcf2cytosure}
      . $NEWLINE;
  SAMPLE_ID:
    while ( my ( $sample_id_index, $sample_id ) =
        each @{ $active_parameter_href->{sample_ids} } )
    {

        vcf2cytosure_convert(
            {
                coverage_file =>
                  $vcf2cytosure_file_info{$sample_id}{in}{q{.tab}},
                FILEHANDLE   => $FILEHANDLE,
                outfile_path => $outfile_path{$sample_id},
                vcf_infile_path =>
                  $vcf2cytosure_file_info{$sample_id}{in}{q{.vcf}},
            }
        );
        say {$FILEHANDLE} $AMPERSAND . $SPACE . $NEWLINE;

        if ( $program_mode == 1 ) {

            add_program_outfile_to_sample_info(
                {
                    infile           => $outfile_name{$sample_id},
                    sample_id        => $sample_id,
                    path             => $outfile_path{$sample_id},
                    program_name     => q{vcf2cytosure},
                    sample_info_href => $sample_info_href,
                }
            );

            ## For logging version - until present in cgh file
            add_program_outfile_to_sample_info(
                {
                    infile           => $outfile_name{$sample_id},
                    sample_id        => $sample_id,
                    path             => $stderrfile_path,
                    program_name     => q{vcf2cytosure_version},
                    sample_info_href => $sample_info_href,
                }
            );
        }
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    if ( $program_mode == 1 ) {

        slurm_submit_job_sample_id_dependency_family_dead_end(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                sbatch_file_name => $file_path,
            }
        );
    }
    return;
}

1;
