package MIP::Recipes::Analysis::Gatk_asereadcounter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use POSIX qw{ floor };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Check::Cluster qw{ check_max_core_number };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.10;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_gatk_asereadcounter };

}

## Constants
Readonly my $ASTERISK => q{*};
Readonly my $NEWLINE  => qq{\n};

sub analysis_gatk_asereadcounter {

## Function : Gatk asereadcounter analysis for rna recipe
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_id               => Sample id
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
        case_id_ref => {
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

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_parameters get_recipe_attributes };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Alignment::Gatk qw{ gatk_asereadcounter };
    use MIP::Program::Variantcalling::Gatk qw{ gatk_indexfeaturefile };
    use MIP::Program::Variantcalling::Bcftools qw{ bcftools_view };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::QC::Sample_info qw{ set_recipe_outfile_in_sample_info };

    ## Constants
    Readonly my $ALLELES                => 2;
    Readonly my $JAVA_MEMORY_ALLOCATION => 20;

    ### PREPROCESSING

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

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
    my $variant_infile_name_prefix    = $io{in}{file_name_prefix};
    my $variant_infile_path           = ${ $io{in}{file_paths} }[0];
    my $variant_suffix                = $io{in}{file_suffix};
    my $temp_variant_file_path_prefix = $io{temp}{file_path_prefix};
    my $temp_variant_file_path        = $temp_variant_file_path_prefix . $variant_suffix;

    ## Get bam infile from GATK BaseRecalibration
    my %alignment_io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => q{gatk_baserecalibration},
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my $alignment_file_path_prefix      = $alignment_io{out}{file_path_prefix};
    my $alignment_suffix                = $alignment_io{out}{file_suffix};
    my $alignment_file_path             = $alignment_file_path_prefix . $alignment_suffix;
    my $temp_alignment_file_path_prefix = $alignment_io{temp}{file_path_prefix};
    my $temp_alignment_file_path = $temp_alignment_file_path_prefix . $alignment_suffix;

    my $job_id_chain = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $recipe_mode        = $active_parameter_href->{$recipe_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};

    ## Get module parameters
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Outpaths
    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => [$variant_infile_name_prefix],
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
                temp_directory         => $temp_directory,
            }
        )
    );
    my $outfile_name = ${ $io{out}{file_names} }[0];
    my $outfile_path = ${ $io{out}{file_paths} }[0];
    my $outdir_path  = $io{out}{dir_path};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $sample_id,
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

    ### SHELL

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy bam file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $alignment_file_path_prefix
              . substr( $alignment_suffix, 0, 2 )
              . $ASTERISK,
            outfile_path => $temp_directory,
        }
    );
    say {$FILEHANDLE} q{wait} . $NEWLINE;

    ## Restrict analysis to biallelic, heterogenous SNPs
    say {$FILEHANDLE} q{## Bcftools view};
    bcftools_view(
        {
            FILEHANDLE   => $FILEHANDLE,
            genotype     => q{het},
            infile_path  => $variant_infile_path,
            max_alleles  => $ALLELES,
            min_alleles  => $ALLELES,
            outfile_path => $temp_variant_file_path,
            types        => q{snps},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Index VCF
    say {$FILEHANDLE} q{## GATK IndexFeatureFile};
    gatk_indexfeaturefile(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $temp_variant_file_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## GATK ASEReadCounter
    say {$FILEHANDLE} q{## GATK ASEReadCounter};
    gatk_asereadcounter(
        {
            FILEHANDLE           => $FILEHANDLE,
            infile_path          => $temp_alignment_file_path,
            java_use_large_pages => $active_parameter_href->{java_use_large_pages},
            memory_allocation    => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
            outfile_path         => $outfile_path,
            referencefile_path   => $active_parameter_href->{human_genome_reference},
            variant_infile_path  => $temp_variant_file_path,
            verbosity            => $active_parameter_href->{gatk_logging_level},
            temp_directory       => $temp_directory,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    close $FILEHANDLE;

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => $outfile_name,
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command            => $profile_base_command,
                dependency_method       => q{sample_to_sample},
                case_id                 => $case_id,
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
    return 1;
}
1;
