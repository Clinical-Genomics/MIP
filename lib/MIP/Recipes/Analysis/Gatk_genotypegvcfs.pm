package MIP::Recipes::Analysis::Gatk_genotypegvcfs;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
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
use MIP::Constants qw{ $DOT $LOG_NAME $NEWLINE $TAB $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.18;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_gatk_genotypegvcfs };

}

## Constants
Readonly my $INCLUDE_NONVARIANT_SITES_TIME => 50;
Readonly my $JAVA_MEMORY_ALLOCATION        => 8;

sub analysis_gatk_genotypegvcfs {

## Function : GATK GenoTypeGVCFs.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending"
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

    use MIP::File::Format::Pedigree qw{ create_fam_file };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Gnu::Coreutils qw{ gnu_cat gnu_echo gnu_rm };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Variantcalling::Gatk
      qw{ gatk_genomicsdbimport  gatk_genotypegvcfs };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my $job_id_chain            = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $recipe_mode          = $active_parameter_href->{$recipe_name};
    my $recipe_files_tracker = 0;

    ## Gatk genotype is most safely processed in single thread mode, but we need some java heap allocation
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number = $recipe_resource{core_number};
    my $time        = $recipe_resource{time};

    ## If all sites should be included
    if ( $active_parameter_href->{gatk_genotypegvcfs_all_sites} == 1 ) {

        # Including all sites requires longer processing time
        $time = $INCLUDE_NONVARIANT_SITES_TIME;
    }

    ## Set and get the io files per chain, id and stream
    my %io = parse_io_outfiles(
        {
            chain_id         => $job_id_chain,
            id               => $case_id,
            file_info_href   => $file_info_href,
            file_name_prefix => $case_id,
            iterators_ref    => $file_info_href->{contigs_size_ordered},
            outdata_dir      => $active_parameter_href->{outdata_dir},
            parameter_href   => $parameter_href,
            recipe_name      => $recipe_name,
        }
    );
    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my %outfile_path        = %{ $io{out}{file_path_href} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

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

  CONTIG:
    foreach my $contig ( @{ $file_info_href->{contigs} } ) {

        ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
        my ($recipe_file_path) = setup_script(
            {
                active_parameter_href           => $active_parameter_href,
                core_number                     => $core_number,
                directory_id                    => $case_id,
                filehandle                      => $filehandle,
                job_id_href                     => $job_id_href,
                log                             => $log,
                memory_allocation               => $recipe_resource{memory},
                process_time                    => $time,
                recipe_directory                => $recipe_name,
                recipe_name                     => $recipe_name,
                sleep                           => 1,
                source_environment_commands_ref => $recipe_resource{load_env_ref},
                temp_directory                  => $temp_directory,
            }
        );

        ## Collect infiles for all sample_ids (WGS)
        my @genotype_infile_paths;

        ## Collect infiles for all samples_ids (WES)
        my @sample_vcf_path_lines;
        my $sample_name_map_path;

      SAMPLE_ID:
        while ( my ( $sample_id_index, $sample_id ) =
            each @{ $active_parameter_href->{sample_ids} } )
        {

            ## Get the io infiles per chain and id
            my %sample_io = get_io_files(
                {
                    id             => $sample_id,
                    file_info_href => $file_info_href,
                    parameter_href => $parameter_href,
                    recipe_name    => $recipe_name,
                    stream         => q{in},
                }
            );
            if ( $consensus_analysis_type eq q{wes} ) {

                push @sample_vcf_path_lines,
                  $sample_id . $TAB . $sample_io{in}{file_path} . $NEWLINE;
            }
            else {

                push @genotype_infile_paths, $sample_io{in}{file_path};
            }
        }

        my $genomicsdb_file_path =
          $outfile_path_prefix . $UNDERSCORE . q{DB} . $UNDERSCORE . $contig;

        ## Remove potential GenomicsDB from previous analysis as this causes
        ## GenomicsDBImport to crash
        say {$filehandle} q{## Remove potential GenomicsDB from previous analysis};
        gnu_rm(
            {
                filehandle  => $filehandle,
                force       => 1,
                infile_path => $genomicsdb_file_path,
                recursive   => 1,
            }
        );
        say {$filehandle} $NEWLINE;

        ## GATK GenomicsDBImport
        say {$filehandle} q{## GATK GenomicsDBImport};

        ## Files to import into GenomicsDB
        if ( $consensus_analysis_type eq q{wes} ) {

            $sample_name_map_path =
              catfile( $outdir_path_prefix, q{analysis_sample_map.txt} );
            my $echo_outfile_path =
              catfile( $outdir_path_prefix, q{dynamic_sample_map.txt} );
            _merge_sample_name_map_files(
                {
                    echo_outfile_path => $echo_outfile_path,
                    filehandle        => $filehandle,
                    gatk_genotypegvcfs_ref_gvcf =>
                      $active_parameter_href->{gatk_genotypegvcfs_ref_gvcf},
                    outfile_path => $sample_name_map_path,
                    strings_ref  => \@sample_vcf_path_lines,
                }
            );

        }

        gatk_genomicsdbimport(
            {
                filehandle                => $filehandle,
                genomicsdb_workspace_path => $genomicsdb_file_path,
                intervals_ref             => [$contig],
                infile_paths_ref          => \@genotype_infile_paths,
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                memory_allocation    => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
                referencefile_path   => $active_parameter_href->{human_genome_reference},
                sample_name_map_path => $sample_name_map_path,
                temp_directory       => $temp_directory,
                verbosity            => $active_parameter_href->{gatk_logging_level},
            }
        );
        say {$filehandle} $NEWLINE;

        ## GATK GenoTypeGVCFs
        say {$filehandle} q{## GATK GenoTypeGVCFs};

        gatk_genotypegvcfs(
            {
                dbsnp_path =>
                  $active_parameter_href->{gatk_haplotypecaller_snp_known_set},
                filehandle => $filehandle,
                include_nonvariant_sites =>
                  $active_parameter_href->{gatk_genotypegvcfs_all_sites},
                infile_path          => q{gendb://} . $genomicsdb_file_path,
                intervals_ref        => [$contig],
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                memory_allocation    => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
                outfile_path         => $outfile_path{$contig},
                pedigree             => $fam_file_path,
                referencefile_path   => $active_parameter_href->{human_genome_reference},
                temp_directory       => $temp_directory,
                verbosity            => $active_parameter_href->{gatk_logging_level},
                use_new_qual_calculator =>
                  $active_parameter_href->{gatk_use_new_qual_calculator},
            }
        );
        say {$filehandle} $NEWLINE;

        close $filehandle;

        if ( $recipe_mode == 1 ) {

            submit_recipe(
                {
                    base_command            => $profile_base_command,
                    case_id                 => $case_id,
                    dependency_method       => q{sample_to_case_parallel},
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_chain            => $job_id_chain,
                    job_id_href             => $job_id_href,
                    job_reservation_name =>
                      $active_parameter_href->{job_reservation_name},
                    log                  => $log,
                    recipe_file_path     => $recipe_file_path,
                    recipe_files_tracker => $recipe_files_tracker,
                    sample_ids_ref       => \@{ $active_parameter_href->{sample_ids} },
                    submission_profile   => $active_parameter_href->{submission_profile},
                }
            );
        }
        $recipe_files_tracker++;    # Tracks nr of recipe scripts
    }
    return 1;
}

sub _merge_sample_name_map_files {

## Function : Merge sample_name_map files
## Returns  :
## Arguments: $filehandle                  => Filehandle to write to
##          : $echo_outfile_path           => Echo outfile path for dynamic samples
##          : $gatk_genotypegvcfs_ref_gvcf => Merged reference sample name map file path
##          : $outfile_path                => Outfile path
##          : $strings_ref                 => Strings to echo {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $echo_outfile_path;
    my $gatk_genotypegvcfs_ref_gvcf;
    my $outfile_path;
    my $strings_ref;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        echo_outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$echo_outfile_path,
            strict_type => 1,
        },
        gatk_genotypegvcfs_ref_gvcf => {
            defined     => 1,
            required    => 1,
            store       => \$gatk_genotypegvcfs_ref_gvcf,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        strings_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$strings_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    gnu_echo(
        {
            enable_interpretation => 1,
            filehandle            => $filehandle,
            no_trailing_newline   => 1,
            outfile_path          => $echo_outfile_path,
            strings_ref           => $strings_ref,
        }
    );
    say {$filehandle} $NEWLINE;

    gnu_cat(
        {
            filehandle       => $filehandle,
            infile_paths_ref => [ $echo_outfile_path, $gatk_genotypegvcfs_ref_gvcf ],
            stdoutfile_path  => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;
    return;
}

1;
