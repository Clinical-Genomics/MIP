package MIP::Recipes::Analysis::Gatk_variantevalexome;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use autodie qw{ :all };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Spec::Functions qw{ catdir catfile };

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_gatk_variantevalexome };

}

##Constants
Readonly my $ASTERIX    => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $PIPE       => q{|};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_gatk_variantevalexome {

## Function : GATK varianteval for exome variants.
## Returns  :
## Arguments: $parameter_href          => Parameter hash {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $sample_id               => Sample id
##          : $infamily_directory      => In family directory
##          : $outsample_directory     => Out sample directory
##          : $program_name            => Program name {REF}
##          : $family_id               => Family id
##          : $temp_directory          => Temporary directory
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $call_type               => Variant call type

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $sample_id;
    my $infamily_directory;
    my $outsample_directory;
    my $program_name;

    ## Default(s)
    my $family_id;
    my $temp_directory;
    my $outaligner_dir;
    my $call_type;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
        },
        infamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infamily_directory,
        },
        outsample_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outsample_directory,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{get_file_suffix get_merged_infile_prefix };
    use MIP::Gnu::Coreutils qw(gnu_cat gnu_sort);
    use MIP::IO::Files qw(migrate_file);
    use MIP::Language::Java qw{java_core};
    use MIP::Script::Setup_script qw(setup_script);
    use MIP::Set::File qw{set_file_suffix };
    use Program::Variantcalling::Bedtools qw(intersectbed);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_family_dead_end);
    use MIP::Program::Variantcalling::Gatk
      qw(gatk_varianteval gatk_selectvariants);
    use MIP::QC::Record qw(add_program_outfile_to_sample_info);

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $gatk_jar =
      catfile( $active_parameter_href->{gatk_path}, q{GenomeAnalysisTK.jar} );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $sample_id,
            program_name          => $program_name,
            program_directory     => catfile( $outaligner_dir, $program_name ),
            call_type             => $call_type,
            core_number           => $core_number,
            process_time          => $time,
            temp_directory        => $temp_directory
        }
    );

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$family_id}{pgatk_combinevariantcallsets}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{pgatk_combinevariantcallsets}{file_tag};

    ## Files
    my $infile_prefix  = $family_id . $infile_tag . $call_type;
    my $outfile_prefix = $merged_infile_prefix . $outfile_tag;

    ## Paths
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
            jobid_chain =>
              $parameter_href->{pgatk_combinevariantcallsets}{chain},
        }
    );
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{variant_eval_file_suffix},
            job_id_chain   => $job_id_chain,
            file_suffix => $parameter_href->{$mip_program_name}{outfile_suffix},
        }
    );

    #### DEVELOP use BCFTOOLS instead
    my $extract_exonic_regexp =
      q?perl -ne ' if ( ($_=~/exonic/) || ($_=~/splicing/) ) {print $_;}' ?;

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => catfile(
                $infamily_directory, $infile_prefix . $infile_suffix . $ASTERIX
            ),
            outfile_path => $temp_directory
        }
    );

    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Select sample id from family id vcf file

    ## GATK SelectVariants
    say {$FILEHANDLE} q{## GATK SelectVariants};

    gatk_selectvariants(
        {
            memory_allocation => q{Xmx2g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory     => $temp_directory,
            java_jar           => $gatk_jar,
            sample_names_ref   => [$sample_id],
            logging_level      => $active_parameter_href->{gatk_logging_level},
            referencefile_path => $referencefile_path,
            infile_path        => $file_path_prefix . $infile_suffix,
            outfile_path       => $outfile_path_prefix
              . $call_type
              . $UNDERSCORE . q{temp}
              . $infile_suffix,
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Update file_tags
    $infile_tag       = $file_info_href->{$family_id}{prankvariant}{file_tag};
    $infile_prefix    = $family_id . $infile_tag . $call_type;
    $file_path_prefix = catfile( $temp_directory, $infile_prefix );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => catfile(
                $infamily_directory, $infile_prefix . $infile_suffix . $ASTERIX
            ),
            outfile_path => $temp_directory
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Extract exonic variants
    say   {$FILEHANDLE} q{## Extract exonic variants};
    print {$FILEHANDLE} $extract_exonic_regexp;

    ## InFile
    print {$FILEHANDLE} $file_path_prefix . $infile_suffix . $SPACE;
    ## OutFile
    say {$FILEHANDLE} q{>}
      . $SPACE
      . catfile(
        $temp_directory,
        $sample_id
          . $infile_tag
          . $call_type
          . $UNDERSCORE
          . q{exonic_variants}
          . $infile_suffix
      ),
      $NEWLINE;

    ## Include potential select file variants
    if ( $active_parameter_href->{vcfparser_outfile_count} == 2 ) {

        ## SelectFile variants
        my $vcfparser_analysis_type = $DOT . q{selected};

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        migrate_file(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => catfile(
                    $infamily_directory,
                    $infile_prefix
                      . $vcfparser_analysis_type
                      . $infile_suffix
                      . $ASTERIX
                ),
                outfile_path => $temp_directory
            }
        );
        say {$FILEHANDLE} q{wait}, $NEWLINE;

        ## Extract exonic variants
        say   {$FILEHANDLE} q{## Extract exonic variants};
        print {$FILEHANDLE} $extract_exonic_regexp;

        ## InFile
        print {$FILEHANDLE} $file_path_prefix
          . $vcfparser_analysis_type
          . $infile_suffix
          . $SPACE;

        ## OutFile
        say {$FILEHANDLE} q{>}
          . $SPACE
          . catfile(
            $temp_directory,
            $sample_id
              . $infile_tag
              . $call_type
              . $UNDERSCORE
              . q{exonic_variants}
              . $vcfparser_analysis_type
              . $infile_suffix
          ),
          $NEWLINE;

        ## Merge orphans and selectfiles
        say {$FILEHANDLE} q{## Merge orphans and selectfile(s)};
        gnu_cat(
            {
                infile_paths_ref => [
                    catfile(
                        $temp_directory,
                        $sample_id
                          . $infile_tag
                          . $call_type
                          . $UNDERSCORE
                          . q{exonic_variants}
                          . $infile_suffix
                    ),
                    catfile(
                        $temp_directory,
                        $sample_id
                          . $infile_tag
                          . $call_type
                          . $UNDERSCORE
                          . q{exonic_variants}
                          . $vcfparser_analysis_type
                          . $infile_suffix
                    )
                ],
                outfile_path => catfile(
                    $temp_directory,
                    $sample_id
                      . $infile_tag
                      . $call_type
                      . $UNDERSCORE
                      . q{exonic_variants_combined}
                      . $infile_suffix
                ),
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Sort combined file
        say {$FILEHANDLE} q{## Sort combined file};
        gnu_sort(
            {
                keys_ref    => [ q{1,1}, q{2,2n} ],
                infile_path => catfile(
                    $temp_directory,
                    $sample_id
                      . $infile_tag
                      . $call_type
                      . $UNDERSCORE
                      . q{exonic_variants_combined}
                      . $infile_suffix
                ),
                outfile_path => catfile(
                    $temp_directory,
                    $sample_id
                      . $infile_tag
                      . $call_type
                      . $UNDERSCORE
                      . q{exonic_variants}
                      . $infile_suffix
                ),
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    print {$FILEHANDLE} q?perl -ne ' if ($_=~/^#/) {print $_;}' ?;

    ## InFile
    print {$FILEHANDLE} $file_path_prefix . $infile_suffix . $SPACE;

    ## Pipe
    print {$FILEHANDLE} $PIPE . $SPACE;

    gnu_cat(
        {
            infile_paths_ref => [
                q{-},
                catfile(
                    $temp_directory,
                    $sample_id
                      . $infile_tag
                      . $call_type
                      . $UNDERSCORE
                      . q{exonic_variants}
                      . $infile_suffix
                )
            ],
            outfile_path => catfile(
                $temp_directory,
                $sample_id
                  . $infile_tag
                  . $call_type
                  . $UNDERSCORE
                  . q{exonic_variants_head}
                  . $infile_suffix
            ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Intersect exonic variants from created sample_id vcf file (required for gatk_varianteval for exonic variants)
    say {$FILEHANDLE}
      q{## Intersect exonic variants from created sample_id vcf file};
    intersectbed(
        {
            with_header => 1,
            infile_path => $outfile_path_prefix
              . $call_type
              . $UNDERSCORE . q{temp}
              . $infile_suffix,
            intersectfile_path => catfile(
                $temp_directory,
                $sample_id
                  . $infile_tag
                  . $call_type
                  . $UNDERSCORE
                  . q{exonic_variants_head}
                  . $infile_suffix
            ),
            outfile_path => $outfile_path_prefix
              . $call_type
              . $UNDERSCORE
              . q{exome}
              . $infile_suffix,
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## VariantEval
    say {$FILEHANDLE} q{## GATK varianteval};

    gatk_varianteval(
        {
            memory_allocation => q{Xmx2g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory   => $temp_directory,
            java_jar         => $gatk_jar,
            infile_paths_ref => [
                    $outfile_path_prefix
                  . $call_type
                  . $UNDERSCORE
                  . q{exome}
                  . $infile_suffix
            ],
            outfile_path => $outfile_path_prefix
              . $call_type
              . $UNDERSCORE
              . q{exome}
              . $outfile_suffix,
            logging_level      => $active_parameter_href->{gatk_logging_level},
            referencefile_path => $referencefile_path,
            dbsnp_file_path => $active_parameter_href->{gatk_varianteval_dbsnp},
            indel_gold_standard_file_path =>
              $active_parameter_href->{gatk_varianteval_gold},
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            infile_path => $outfile_path_prefix
              . $call_type
              . $UNDERSCORE
              . q{exome}
              . $outfile_suffix,
            outfile_path => $outsample_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $mip_program_mode == 1 ) {

        ## Collect QC metadata info for later use
        my $qc_exome_outfile =
            $outfile_prefix
          . $call_type
          . $UNDERSCORE
          . q{exome}
          . $outfile_suffix;
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                sample_id        => $sample_id,
                program_name     => q{variantevalexome},
                infile           => $merged_infile_prefix,
                outdirectory     => $outsample_directory,
                outfile          => $qc_exome_outfile,
            }
        );

        slurm_submit_job_sample_id_dependency_family_dead_end(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $family_id,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
    return;
}

1;
