package MIP::Recipes::Picardtools_genotypeconcordance;

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
use File::Spec::Functions qw{ catdir catfile devnull };

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_picardtools_genotypeconcordance };

}

##Constants
Readonly my $ASTERIX    => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_picardtools_genotypeconcordance {

## Function : Compare metrics for this analysis run with the NIST reference dataset.
## Returns  :
## Arguments: $parameter_href          => Parameter hash {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $sample_id               => Sample id
##          : $infamily_directory      => In family directory
##          : $outfamily_directory     => Out family directory
##          : $program_name            => Program name
##          : $family_id               => Family id
##          : $temp_directory          => Temporary directory
##          : $reference_dir           => MIP reference directory
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $call_type               => Variant call type

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id;
    my $temp_directory;
    my $reference_dir;
    my $outaligner_dir;
    my $call_type;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $sample_id;
    my $infamily_directory;
    my $outfamily_directory;
    my $program_name;

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
            store       => \$sample_id,
        },
        infamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infamily_directory,
        },
        outfamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfamily_directory,
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
        reference_dir => {
            default     => $arg_href->{active_parameter_href}{reference_dir},
            strict_type => 1,
            store       => \$reference_dir,
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

    use MIP::Gnu::Coreutils qw{ gnu_cat };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_family_dead_end };
    use MIP::Program::Interval::Picardtools qw{ picardtools_intervallisttools };
    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_stats bcftools_rename_vcf_samples };
    use MIP::Program::Variantcalling::Gatk
      qw{ gatk_selectvariants gatk_leftalignandtrimvariants };
    use MIP::Program::Variantcalling::Picardtools
      qw{ picardtools_genotypeconcordance };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $family_id,
            program_name          => $program_name,
            program_directory     => catfile( $outaligner_dir, $program_name ),
            core_number           => $core_number,
            process_time          => $time,
            temp_directory        => $temp_directory,
            set_errexit => 0,    # Special case to allow "vcf.idx" to be created
            error_trap  => 0,    # Special case to allow "vcf.idx" to be created
        }
    );

    ## Used downstream
    $parameter_href->{$mip_program_name}{indirectory} = $outfamily_directory;

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$family_id}{pgatk_combinevariantcallsets}{file_tag};
    my $infile_prefix = $family_id . $infile_tag . $call_type;
    my $file_path_prefix = catfile( $temp_directory, $infile_prefix );
    my $outfile_prefix =
        $family_id
      . $DOT . q{Vs}
      . $DOT
      . $active_parameter_href->{nist_id}
      . q{-NIST_genome};
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix ),
      my $call_file_path    = catfile( $temp_directory, $family_id );
    my $nist_file_path = catfile( $temp_directory, q{NIST} );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => catfile(
                $infamily_directory, $infile_prefix . $DOT . q{vcf} . $ASTERIX
            ),
            outfile_path => $temp_directory
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## Rename vcf samples. The samples array will replace the sample names in the same order as supplied.
    bcftools_rename_vcf_samples(
        {
            sample_ids_ref => [ $active_parameter_href->{nist_id} . q{-NIST} ],
            temp_directory => $temp_directory,
            infile  => $active_parameter_href->{nist_high_confidence_call_set},
            outfile => $nist_file_path . $UNDERSCORE . q{refrm.vcf},
            FILEHANDLE => $FILEHANDLE,
        }
    );

    ## Modify since different ref genomes
    say {$FILEHANDLE} q{## Modify since different ref genomes};

    ## Do not print GL contigs
    print {$FILEHANDLE}
      q?perl -nae 'unless($_=~/##contig=<ID=GL\d+/) {print $_}' ?;

    ## Infile
    print {$FILEHANDLE} $nist_file_path . $UNDERSCORE . q{refrm.vcf} . $SPACE;

    ## Outfile
    print {$FILEHANDLE} q{>}
      . $SPACE
      . $nist_file_path
      . $DOT . q{vcf}
      . $SPACE;
    say {$FILEHANDLE} $NEWLINE;

    ## Bcftools stats
    say {$FILEHANDLE} q{## bcftools stats};
    bcftools_stats(
        {
            infile_path  => $nist_file_path . $DOT . q{vcf},
            outfile_path => $nist_file_path . $DOT . q{vcf.stats},
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Generate ".idx" for downstream Picard by failling this process
    say {$FILEHANDLE}
      q{## Generate '.idx' for downstream Picard by failling this process};

    gatk_selectvariants(
        {
            memory_allocation => q{Xmx2g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $temp_directory,
            java_jar       => catfile(
                $active_parameter_href->{gatk_path},
                q{GenomeAnalysisTK.jar}
            ),
            sample_names_ref   => [ $sample_id . q{XXX} ],
            logging_level      => $active_parameter_href->{gatk_logging_level},
            referencefile_path => $referencefile_path,
            infile_path        => $nist_file_path . $DOT . q{vcf},
            outfile_path       => $nist_file_path . q{XXX.vcf},
            FILEHANDLE         => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Create .interval_list file from NIST union bed
    say {$FILEHANDLE} q{## Prepare .interval_list file from NIST union bed},
      $NEWLINE;

    my $genome_dict_file_path = catfile( $temp_directory,
            $file_info_href->{human_genome_reference_name_prefix}
          . $DOT
          . q{dict} );
    ## Do not print contigs
    print {$FILEHANDLE}
q?perl -nae 'unless($_=~/NC_007605/ || $_=~/hs37d5/ || $_=~/GL\d+/) {print $_}' ?;
    print {$FILEHANDLE} catfile( $reference_dir,
        $file_info_href->{human_genome_reference_name_prefix} . $DOT . q{dict} )
      . $SPACE;
    print {$FILEHANDLE} q{>} . $SPACE . $genome_dict_file_path . $SPACE;
    say   {$FILEHANDLE} $NEWLINE;

    gnu_cat(
        {
            infile_paths_ref => [
                $genome_dict_file_path,
                $active_parameter_href->{nist_high_confidence_call_set_bed}
            ],
            outfile_path => $nist_file_path . $DOT . q{bed.dict_body},
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE}
q{## Remove target annotations, 'track', 'browse' and keep only 5 columns};
    print {$FILEHANDLE}
q?perl  -nae 'if ($_=~/@/) {print $_;} elsif ($_=~/^track/) {} elsif ($_=~/^browser/) {} else {print @F[0], "\t", (@F[1] + 1), "\t", @F[2], "\t", "+", "\t", "-", "\n";}' ?;
    ## Infile
    print {$FILEHANDLE} $nist_file_path . $DOT . q{bed.dict_body} . $SPACE;

    ## Remove unnecessary info and reformat
    print {$FILEHANDLE} q{>}
      . $SPACE
      . $nist_file_path
      . $DOT
      . q{bed.dict_body_col_5.interval_list}
      . $SPACE;
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Create }
      . $active_parameter_href->{nist_high_confidence_call_set_bed}
      . $DOT
      . q{interval_list};

    picardtools_intervallisttools(
        {
            memory_allocation => q{Xmx2g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $active_parameter_href->{temp_directory},
            java_jar       => catfile(
                $active_parameter_href->{picardtools_path},
                q{picard.jar}
            ),
            infile_paths_ref =>
              [ $nist_file_path . $DOT . q{bed.dict_body_col_5.interval_list} ],
            outfile_path       => $nist_file_path . $DOT . q{bed.interval_list},
            referencefile_path => $referencefile_path,
            FILEHANDLE         => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ### MIP data
    ## GATK SelectVariants
    say {$FILEHANDLE} q{## GATK SelectVariants};

    gatk_selectvariants(
        {
            memory_allocation => q{Xmx2g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $temp_directory,
            java_jar       => catfile(
                $active_parameter_href->{gatk_path},
                q{GenomeAnalysisTK.jar}
            ),
            sample_names_ref    => [$sample_id],
            logging_level       => $active_parameter_href->{gatk_logging_level},
            referencefile_path  => $referencefile_path,
            infile_path         => $file_path_prefix . $DOT . q{vcf},
            outfile_path        => $call_file_path . $DOT . q{vcf},
            exclude_nonvariants => 1,
            FILEHANDLE          => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Left align, trim and split allels
    say {$FILEHANDLE} q{## GATK LeftAlignAndTrimVariants};

    gatk_leftalignandtrimvariants(
        {
            memory_allocation => q{Xmx2g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $temp_directory,
            java_jar       => catfile(
                $active_parameter_href->{gatk_path},
                q{GenomeAnalysisTK.jar}
            ),
            infile_path         => $call_file_path . $DOT . q{vcf},
            logging_level       => $active_parameter_href->{gatk_logging_level},
            referencefile_path  => $referencefile_path,
            split_multiallelics => 1,
            outfile_path        => $call_file_path . $UNDERSCORE . q{lts.vcf},
            FILEHANDLE          => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Modify since different ref genomes
    say   {$FILEHANDLE} q{## Modify since different ref genomes};
    print {$FILEHANDLE}
q?perl -nae 'unless($_=~/##contig=<ID=NC_007605,length=171823>/ || $_=~/##contig=<ID=hs37d5,length=35477943>/ || $_=~/##contig=<ID=GL\d+/) {print $_}' ?;

    ## Infile
    print {$FILEHANDLE} $call_file_path . $UNDERSCORE . q{lts.vcf} . $SPACE;

    ## Outfile
    print {$FILEHANDLE} q{>}
      . $SPACE
      . $call_file_path
      . $UNDERSCORE
      . q{lts_refrm.vcf}
      . $SPACE;
    say {$FILEHANDLE} $NEWLINE;

    ## BcfTools Stats
    say {$FILEHANDLE} q{## bcftools stats};
    bcftools_stats(
        {
            infile_path  => $call_file_path . $UNDERSCORE . q{lts_refrm.vcf},
            outfile_path => $call_file_path
              . $UNDERSCORE
              . q{lts_refrm.vcf.stats},
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Generate ".idx" for downstream Picard by failling this process
    say {$FILEHANDLE}
      q{## Generate '.idx' for downstream Picard by failling this process};

    gatk_selectvariants(
        {
            memory_allocation => q{Xmx2g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $temp_directory,
            java_jar       => catfile(
                $active_parameter_href->{gatk_path},
                q{GenomeAnalysisTK.jar}
            ),
            sample_names_ref   => [ $sample_id . q{XXX} ],
            logging_level      => $active_parameter_href->{gatk_logging_level},
            referencefile_path => $referencefile_path,
            infile_path  => $call_file_path . $UNDERSCORE . q{lts_refrm.vcf},
            outfile_path => $call_file_path . q{XXX.vcf},
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE}
q{## Picard GenotypeConcordance - Genome restricted by union - good quality};
    picardtools_genotypeconcordance(
        {
            intervals_ref => [ $nist_file_path . $DOT . q{bed.interval_list} ],
            infile_path   => $call_file_path . $UNDERSCORE . q{lts_refrm.vcf},
            truth_file_path     => $nist_file_path . $DOT . q{vcf},
            outfile_prefix_path => $outfile_path_prefix . $UNDERSCORE . q{bed},
            truth_sample        => $active_parameter_href->{nist_id} . q{-NIST},
            call_sample         => $sample_id,
            min_genotype_quality => 20,
            min_depth            => 10,
            FILEHANDLE           => $FILEHANDLE,
            referencefile_path   => $referencefile_path,
            memory_allocation    => q{Xmx2g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $active_parameter_href->{temp_directory},
            java_jar       => catfile(
                $active_parameter_href->{picardtools_path},
                q{picard.jar}
            ),
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Picard GenotypeConcordance - Genome - good quality};
    picardtools_genotypeconcordance(
        {
            infile_path     => $call_file_path . $UNDERSCORE . q{lts_refrm.vcf},
            truth_file_path => $nist_file_path . $DOT . q{vcf},
            outfile_prefix_path => $outfile_path_prefix,
            truth_sample        => $active_parameter_href->{nist_id} . q{-NIST},
            call_sample         => $sample_id,
            min_genotype_quality => 20,
            min_depth            => 10,
            FILEHANDLE           => $FILEHANDLE,
            referencefile_path   => $referencefile_path,
            memory_allocation    => q{Xmx2g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory => $active_parameter_href->{temp_directory},
            java_jar       => catfile(
                $active_parameter_href->{picardtools_path},
                q{picard.jar}
            ),
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    my @outfiles = ( $outfile_prefix . $ASTERIX, $ASTERIX . q{vcf.stats} );
  OUTFILE:
    foreach my $outfile (@outfiles) {

        migrate_file(
            {
                infile_path  => catfile( $temp_directory, $outfile ),
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {

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
