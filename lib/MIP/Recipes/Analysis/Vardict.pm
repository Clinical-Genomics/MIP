package MIP::Recipes::Analysis::Vardict;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use List::MoreUtils qw { any };
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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_vardict_tn };

}

## Constants
Readonly my $ASTERISK     => q{*};
Readonly my $DOUBLE_QUOTE => q{"};
Readonly my $NEWLINE      => qq{\n};
Readonly my $PIPE         => q{|};
Readonly my $SPACE        => q{ };
Readonly my $STDOUT       => q{>};
Readonly my $UNDERSCORE   => q{_};

sub analysis_vardict_tn {

## Function : Analysis recipe for varDict variant calling
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $insample_directory      => In sample directory
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outsample_directory     => Out sample directory
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $sample_origin           => Info on sample origin (tumor vs normal)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $outfamily_directory;
    my $outsample_directory;
    my $parameter_href;
    my $program_name;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
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
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir,
            strict_type => 1,
        },
        outfamily_directory => {
            defined     => 1,
            required    => 1,
            store       => \$outfamily_directory,
            strict_type => 1,
        },
        outsample_directory => {
            defined     => 1,
            required    => 1,
            store       => \$outsample_directory,
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
            strict_type => 1,
            store       => \$temp_directory,
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Cluster qw{ check_max_core_number };
    use MIP::Cluster qw{ get_core_number };
    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files
      qw{ migrate_file migrate_files xargs_migrate_contig_files };
    use MIP::Program::Variantcalling::Bcftools qw{ bcftools_filter };
    use MIP::Program::Variantcalling::Vardict
      qw{ vardict vardict_var2vcf_single vardict_var2vcf_paired vardict_testsomatic vardict_teststrandbias };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::QC::Record
      qw{ add_program_metafile_to_sample_info add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw { set_file_suffix };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    ## Add get core number not more than max
    $core_number = get_core_number(
        {
            modifier_core_number =>
              scalar( @{ $active_parameter_href->{sample_ids} } ),
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $family_id,
            FILEHANDLE            => $FILEHANDLE,
            job_id_href           => $job_id_href,
            process_time          => $time,
            program_directory     => catfile( $outaligner_dir, q{vardict} ),
            program_name          => $program_name,
            source_environment_commands_ref => [$source_environment_cmd],
            temp_directory                  => $temp_directory,
        }
    );

    ## Get infile_suffix from baserecalibration jobid chain
    my $infile_suffix = get_file_suffix(
        {
            jobid_chain    => $parameter_href->{pgatk_baserecalibration}{chain},
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
        }
    );

    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            file_suffix => $parameter_href->{$mip_program_name}{outfile_suffix},
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{outfile_suffix},
        }
    );

    ## Assign file_tags
    my $outfile_tag =
      $file_info_href->{$family_id}{$mip_program_name}{file_tag};

    my $outfile_prefix = $family_id . $outfile_tag;

    ## Files
    my $outfile_name = $outfile_prefix . $outfile_suffix;

    ## Paths
    my $outfile_path = catfile( $outsample_directory, $outfile_name );

    my %file_path_prefix;

  SAMPLE:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Assign directories
        my $insample_directory = catdir( $active_parameter_href->{outdata_dir},
            $sample_id, $outaligner_dir );

        ## Add merged infile name prefix after merging all BAM files per sample_id
        my $merged_infile_prefix = get_merged_infile_prefix(
            {
                file_info_href => $file_info_href,
                sample_id      => $sample_id,
            }
        );

        ## Assign file_tags
        my $infile_tag =
          $file_info_href->{$sample_id}{pgatk_baserecalibration}{file_tag};

        ## Files
        my $infile_prefix = $merged_infile_prefix . $infile_tag;

        my $infile_name = $infile_prefix . $infile_suffix;

        my $infile_path = catfile( $insample_directory, $infile_name );

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $infile_path,
                outfile_path => $temp_directory,
            }
        );

        ## Paths
        $file_path_prefix{$sample_id} =
          catfile( $temp_directory, $infile_prefix );

    }

    ## Prepare an array of inputs for vardict
    my @infile_path =
      map { $file_path_prefix{$_} . $UNDERSCORE . $infile_suffix }
      @{ $active_parameter_href->{sample_ids} };

    ## Check if the array size is less than two
    ## TODO: sub in vardict module to check for size of array and prepare Tumor vs Normal set.
    ## TODO: sub in vardict module to create pairs of analysis: TvsN and NvsN.
    ## TODO: add merge bam files for more than one normal sample id. This should be added to picardtools_mergesamfiles as a sub.

    ## Vardict
    ## Note: Vardict doesn't runs the analysis based on the input bed file, so if the input bed file points towards the
    ## whole genome, Vardict will run it on whole genome, but mind you, this might take a long time (up to 3 weeks). So
    ## in order to reduce this time two measures need to be taken: 1) use vardict-java to support multiple threads, 2)
    ## run bedtools genomecov and run it on callable regions the same way GATK is doing. After having the callable
    ## regions, by having multiple bed files, it can create a more proper form of parallelization. For now vardict's
    ## recipe is designed to accommodate small bed files with small regions (it might be good to have the regions
    ## overlap with each other).

    ## Create input filename bam file for $infile_paths_ref in vardict module. This has to be an array where first
    ## is tumor and the second one is normal tissue.
    ## Get parameters
    # Unpack
    my $af_threshold       = $active_parameter_href->{vrd_af_threshold};
    my $chrom_start        = $active_parameter_href->{vrd_chrom_start};
    my $input_bed_file     = $active_parameter_href->{vrd_input_bed_file};
    my $max_mm             = $active_parameter_href->{vrd_max_mm};
    my $max_pval           = $active_parameter_href->{vrd_max_pval};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $region_end         = $active_parameter_href->{vrd_region_end};
    my $region_start       = $active_parameter_href->{vrd_region_start};
    my $somatic_only       = $active_parameter_href->{vrd_somatic_only};

    ## Assign family_id to sample_name
    my $sample_name    = $active_parameter_href->{family_id};
    my $segment_annotn = $active_parameter_href->{vrd_segment_annotn};

    say {$FILEHANDLE} q{## varDict variant calling};

    vardict(
        {
            af_threshold           => $af_threshold,
            FILEHANDLE             => $FILEHANDLE,
            infile_paths_ref       => \@infile_path,
            infile_bed_region_info => $input_bed_file,
            max_pval               => $max_pval,
            max_mm                 => $max_mm,
            out_chrom_start        => $chrom_start,
            out_region_start       => $region_start,
            out_region_end         => $region_end,
            out_segment_annotn     => $segment_annotn,
            referencefile_path     => $referencefile_path,
            sample_name            => $sample_name,
            somatic_only           => $somatic_only,
        }
    );

    print {$FILEHANDLE} $PIPE . $SPACE;

    vardict_teststrandbias(
        {
            FILEHANDLE => $FILEHANDLE,
        }
    );

    print {$FILEHANDLE} $PIPE . $SPACE;

    vardict_var2vcf_paired(
        {
            af_threshold    => $af_threshold,
            FILEHANDLE      => $FILEHANDLE,
            max_pval        => $max_pval,
            max_mm          => $max_mm,
            sample_name     => $sample_name,
            somatic_only    => $somatic_only,
            stdoutfile_path => $outfile_path,
        }
    );

    print {$FILEHANDLE} $PIPE . $SPACE;

    ## Filter for low allele frequency with low depth
    bcftools_filter(
        {
            exclude     => _bcftools_low_freq_low_depth(),
            filter_mode => q{+},
            soft_filter => q{LowAlleleDepth},
            FILEHANDLE  => $FILEHANDLE,
        }
    );

    print {$FILEHANDLE} $PIPE . $SPACE;

    ## Filter for low frequency with poor quality
    bcftools_filter(
        {
            exclude     => _bcftools_low_freq_low_qual(),
            filter_mode => q{+},
            soft_filter => q{LowFreqQuality},
            FILEHANDLE  => $FILEHANDLE,
        }
    );

    print {$FILEHANDLE} $PIPE . $SPACE;

    ## Include only those with QUAL > 0
    bcftools_filter(
        {
            include    => $DOUBLE_QUOTE . q{QUAL>=0} . $DOUBLE_QUOTE,
            FILEHANDLE => $FILEHANDLE,
        }
    );

    print {$FILEHANDLE} _tag_somatic() . $PIPE . $SPACE;

    ## Filter for rejected non somatic
    bcftools_filter(
        {
            exclude => q?'STATUS !~ \".*Somatic\"'?,

            #TODO add the following to bcftools
            filter_mode => q{+},
            filter_name => q{REJECT},
            FILEHANDLE  => $FILEHANDLE,
        }
    );

    print {$FILEHANDLE} $PIPE . $SPACE;

    print {$FILEHANDLE} _replace_ambiguous() . $PIPE . $SPACE;

    print {$FILEHANDLE} _filter_duplicate()
      . $SPACE
      . $STDOUT
      . print {$outfile_path}
      . $NEWLINE;

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $outfile_prefix . $outfile_suffix . $ASTERISK,
            outfile_path => $outfamily_directory,
        }
    );
    say {$FILEHANDLE} q{wait} . $NEWLINE;

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});
      or $log->logcroak(q{Could not close XARGSFILEHANDLE});

    if ( $mip_program_mode == 1 ) {

        add_program_outfile_to_sample_info(
            {
                path => catfile(
                    $outfamily_directory, $outfile_prefix . $outfile_suffix
                ),
                program_name     => q{vardict},
                sample_info_href => $sample_info_href,
            }
        );

        slurm_submit_job_sample_id_dependency_add_to_family(
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

sub _bcftools_low_freq_low_qual {

    ## Function: Create filter expression for bcftools.
    ## Return: $bcftools_filter_expr

    ## The filter is: AF < 0.2 && QUAL < 55 && SSF > 0.06
    ## Notes: This filter is aimed at variants with low frequency and low quality, essentially removing false positives (SSF>0.06).

    Readonly my $CONDITON_AND => q{&&};

    my $bcftools_filter_expr;

    my $af_expr = q{AF < 0.2};

    my $pval_expr = q{SSF > 0.06};

    my $qual_expr = q{QUAL < 55};

    $bcftools_filter_expr =
        $DOUBLE_QUOTE
      . $af_expr
      . $CONDITON_AND
      . $qual_expr
      . $CONDITON_AND
      . $pval_expr;

    return $bcftools_filter_expr;
}

sub _bcftools_low_freq_low_depth {

    ## Function: Create filter expression for bcftools. We are using only one aligner in MIP pipeline.
    ## Return: $bcftools_filter_expr

    ## The filter is:
    ## ((AF * DP < 6) &&
    ## ((MQ < 55.0 && NM > 1.0) ||
    ##  (MQ < 60.0 && NM > 2.0) ||
    ##  (DP < 10) ||
    ##  (QUAL < 45)))

    ## Notes:    - AF * DP < 6: According to [ ref 1 ], a read depth of 13x is required to detect heterozygous SNVs 95% of the time.
    ##           - (MQ<55.0 && NM>1.0) || (MQ<60.0 && NM>2.0): It seems bwa mem, mentioned by Brad Chapman in bcbio, and
    ##           also found in [ ref 2] has a maximum mapping quality of 60. This coupled with number of read mismatches (NM) will
    ##           reduce number of false positives
    ##           - DP<10: This is tied to AF*DP above. As it can be seen in [ ref 1 ], at 10 DP, ~90% of SNVs detected correctly
    ##           - QUAL < 45: This is the used threshold in bcbio and ALASCCA pipeline. There is no publication to back this up, although based on
    ##           [ ref 3 ], for rare variants at quality of 60 is the bare minimum. So QUAL < 45 seems to be generalized relax cutoff to filter
    ##           obvious false positives.
    ##           ref 1: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-195
    ##           ref 2: https://github.com/lh3/bwa/blob/b58281621136a0ce2a66837ba509716c727b9387/bwamem.c#L972
    ##           ref 3: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-247

    Readonly my $CONDITON_AND => q{&&};
    Readonly my $CONDITON_OR  => q{||};
    Readonly my $OPEN_PARAN   => q{(};
    Readonly my $CLOSE_PARAN  => q{)};

    my $bcftools_filter_expr;

    my $af_dp_expr = q{(AF * DP < 6)};

    my $mq_nm_expr_1 = q{(MQ < 55.0 && NM > 1.0)};

    my $mq_nm_expr_2 = q{(MQ < 60.0 && NM > 2.0)};

    my $dp_expr = q{(DP < 10)};

    my $qual_expr = q{(QUAL < 45)};

    $bcftools_filter_expr =
        $DOUBLE_QUOTE
      . $OPEN_PARAN
      . $af_dp_expr
      . $CONDITON_AND
      . $OPEN_PARAN
      . $mq_nm_expr_1
      . $CONDITON_OR
      . $mq_nm_expr_2
      . $CONDITON_OR
      . $dp_expr
      . $CONDITON_OR
      . $qual_expr
      . $CLOSE_PARAN
      . $CLOSE_PARAN;

    return $bcftools_filter_expr;
}

sub _filter_duplicate {

    ## Function: Writes a awk expression to open $FILEHANDLE. The awk expression will remove duplicates by comparing ref and alt allele. Adopted from bcbio.
    ## Returns: $awk_filter
    my $awk_filter;

    $awk_filter =
      q?awk -F$'\t' -v OFS='\t' '$1!~/^#/ && $4 == $5 ? . q?{next} {print}'?;

    return $awk_filter;
}

sub _replace_ambiguous {

    ## Function: Writes a awk expression to open $FILEHANDLE. The awk expression will replace ambiguous reference with N. Adopted from bcbio.
    ## Returns: $awk_filter

    ## Notes: - The two awk expressions below will specifically look into column ref and alt alleles, $awk_filter_ref and $awk_filter_alt respectively.
    ##        If the ref/alt alleles have IUPAC nucletotide codes: KMRYSWBVHDX [ https://www.bioinformatics.org/sms/iupac.html ], it will be replaced with N.
    ##        - "if ($0 !~ /^#/)": Lines not starting with "#" essentially the VCF header.
    ##        - "gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", $4)": Replacing IUPAC characters at column 4 with N (any base).
    ##        - "{ print }": Print!

    ## Remove abimguous from ref, column 4 in VCF
    my $awk_filter_ref =
        q?awk -F$'\t' -v OFS='\t' ' ?
      . q?{ if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", $4)} ?
      . q?{ print }'?;

    ## Remove abimguous from alt, column 5 in VCF
    my $awk_filter_alt =
        q?awk -F$'\t' -v OFS='\t' ' ?
      . q?{ if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", $5)} ?
      . q?{ print }'?;

    my $awk_filter =
      $awk_filter_ref . $SPACE . $PIPE . $SPACE . $awk_filter_alt;

    return $awk_filter;
}

sub _tag_somatic {

    ## Function: Writes a sed expression to open $FILEHANDLE. The sed expression will update the somatic. Adopted from bcbio.
    ## Returns: $sed_filter

    my $somatic_replace = q?sed 's/\\\\.*Somatic\\\\/Somatic/'?;

    my $reject_update =
q?sed 's/REJECT,Description=\".*\">/REJECT,Description=\"Not Somatic via VarDict\">/?;

    my $sed_string =
      $somatic_replace . $SPACE . $PIPE . $SPACE . $reject_update;

    return $sed_string;
}

1;
