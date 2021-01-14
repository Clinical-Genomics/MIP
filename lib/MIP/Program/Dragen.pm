package MIP::Program::Dragen;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ dragen_dna_analysis
      dragen_build_hash_table
      dragen_load_hash_table };
}

sub dragen_dna_analysis {

## Function : Perl wrapper for a dragen DNA analysis commands. Dragen version 3.3.5.
## Returns  : @commands
## Arguments: $alignment_output_format       => Alignemnt output format (SAM|BAM|CRAM)
##          : $alt_aware                     => ALT-Aware mapping with a liftover reference
##          : $bam_file_path                 => BAM file path
##          : $cnv_enable_self_normalization => Enable cnv normalization on the fly
##          : $combine_samples_by_name       => Read multiple fastq files by the sample name given in the file name
##          : $dbsnp_file_path               => Dbsnp file path [.vcf|.vcf.gz]
##          : $disable_vcf_compression       => Disable vcf compression (.gz)
##          : $dragen_hash_ref_dir_path      => Dragen reference genome dir path
##          : $enable_bam_indexing           => Enable indexing of BAM file
##          : $enable_cnv                    => Enable cnv hash table creation
##          : $enable_duplicate_marking      => Enable duplication marking in BAM
##          : $enable_combinegvcfs           => Enable generation of a single gVCF file that represents all the input gVCFfiles
##          : $enable_joint_genotyping       => Enable joint genotyping
##          : $enable_map_align              => Enable mapping and alignment
##          : $enable_map_align_output       => Enable map align output
##          : $enable_multi_sample_gvcf      => Enable multi sample gVCF generation
##          : $enable_sampling               => Enable automatic sampling of the insert-length distribution
##          : $enable_sort                   => Enable sorting of alignment file
##          : $enable_variant_caller         => Enable variant caller
##          : $fastq_infile_path             => Infile path (fastq read 1)
##          : $fastq_list_all_samples        => Process all samples together in the same run
##          : $fastq_list_file_path          => Fastq list file path [csv_file]
##          : $fastq_list_sample_id          => Fastq list sample ID
##          : $filehandle                    => Filehandle to write to
##          : $force                         => Force overwrite of existing files
##          : $is_fastq_interleaved          => Smart pairing
##          : $outdirectory_path             => Outdirectory path
##          : $outfile_prefix                => Outfile prefix
##          : $pedigree_file_path            => Pedigree file path
##          : $sample_gvcf_file_paths_ref    => Paths to sample gVCF files
##          : $sample_id                     => Sample ID
##          : $second_fastq_infile_path      => Second infile path (fastq read 2)
##          : $stderrfile_path               => Stderrfile path
##          : $stderrfile_path_append        => Append stderr info to file path
##          : $stdinfile_path                => Stdinfile path
##          : $stdoutfile_path               => Stdoutfile path
##          : $vc_emit_ref_confidence        => Enable gVCF generation
##          : $vc_enable_gatk_acceleration   => Run in GATK mode
##          : $vc_target_bed_file_path       => Restricts processing to regions specified in the BED file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $alignment_output_format;
    my $bam_file_path;
    my $dbsnp_file_path;
    my $dragen_hash_ref_dir_path;
    my $fastq_infile_path;
    my $fastq_list_file_path;
    my $fastq_list_sample_id;
    my $filehandle;
    my $is_fastq_interleaved;
    my $outdirectory_path;
    my $outfile_prefix;
    my $pedigree_file_path;
    my $sample_gvcf_file_paths_ref;
    my $sample_id;
    my $second_fastq_infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;
    my $vc_target_bed_file_path;

    ## Default(s)
    my $alt_aware;
    my $cnv_enable_self_normalization;
    my $combine_samples_by_name;
    my $disable_vcf_compression;
    my $enable_bam_indexing;
    my $enable_cnv;
    my $enable_duplicate_marking;
    my $enable_combinegvcfs;
    my $enable_joint_genotyping;
    my $enable_map_align;
    my $enable_map_align_output;
    my $enable_multi_sample_gvcf;
    my $enable_sampling;
    my $enable_sort;
    my $enable_variant_caller;
    my $fastq_list_all_samples;
    my $force;
    my $vc_emit_ref_confidence;
    my $vc_enable_gatk_acceleration;

    my $tmpl = {
        alignment_output_format => {
            allow       => [ undef, qw{ SAM BAM CRAM } ],
            default     => 0,
            store       => \$alignment_output_format,
            strict_type => 1,
        },
        alt_aware => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$alt_aware,
            strict_type => 1,
        },
        bam_file_path => {
            store       => \$bam_file_path,
            strict_type => 1,
        },
        cnv_enable_self_normalization => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$cnv_enable_self_normalization,
            strict_type => 1,
        },
        combine_samples_by_name => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$combine_samples_by_name,
            strict_type => 1,
        },
        dbsnp_file_path => {
            store       => \$dbsnp_file_path,
            strict_type => 1,
        },
        disable_vcf_compression => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$disable_vcf_compression,
            strict_type => 1,
        },
        dragen_hash_ref_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$dragen_hash_ref_dir_path,
            strict_type => 1,
        },
        enable_bam_indexing => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$enable_bam_indexing,
            strict_type => 1,
        },
        enable_cnv => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$enable_cnv,
            strict_type => 1,
        },
        enable_duplicate_marking => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$enable_duplicate_marking,
            strict_type => 1,
        },
        enable_combinegvcfs => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$enable_combinegvcfs,
            strict_type => 1,
        },
        enable_joint_genotyping => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$enable_joint_genotyping,
            strict_type => 1,
        },
        enable_map_align => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$enable_map_align,
            strict_type => 1,
        },
        enable_map_align_output => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$enable_map_align_output,
            strict_type => 1,
        },
        enable_multi_sample_gvcf => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$enable_multi_sample_gvcf,
            strict_type => 1,
        },
        enable_sampling => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$enable_sampling,
            strict_type => 1,
        },
        enable_sort => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$enable_sort,
            strict_type => 1,
        },
        enable_variant_caller => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$enable_variant_caller,
            strict_type => 1,
        },
        fastq_infile_path => {
            store       => \$fastq_infile_path,
            strict_type => 1,
        },
        fastq_list_all_samples => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$fastq_list_all_samples,
            strict_type => 1,
        },
        fastq_list_file_path => {
            store       => \$fastq_list_file_path,
            strict_type => 1,
        },
        fastq_list_sample_id => {
            store       => \$fastq_list_sample_id,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        force => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$force,
            strict_type => 1,
        },
        is_fastq_interleaved => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$is_fastq_interleaved,
            strict_type => 1,
        },
        outdirectory_path => {
            defined     => 1,
            required    => 1,
            store       => \$outdirectory_path,
            strict_type => 1,
        },
        outfile_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_prefix,
            strict_type => 1,
        },
        pedigree_file_path => {
            store       => \$pedigree_file_path,
            strict_type => 1,
        },
        sample_gvcf_file_paths_ref => {
            default     => [],
            store       => \$sample_gvcf_file_paths_ref,
            strict_type => 1,
        },
        sample_id => {
            store       => \$sample_id,
            strict_type => 1,
        },
        second_fastq_infile_path => {
            store       => \$second_fastq_infile_path,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
        stdinfile_path  => { store => \$stdinfile_path, strict_type => 1, },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        vc_emit_ref_confidence => {
            allow       => [ undef, 0, q{GVCF} ],
            default     => 0,
            store       => \$vc_emit_ref_confidence,
            strict_type => 1,
        },
        vc_enable_gatk_acceleration => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$vc_enable_gatk_acceleration,
            strict_type => 1,
        },
        vc_target_bed_file_path => {
            store       => \$vc_target_bed_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ dragen };

    if ($force) {

        push @commands, q{--force};
    }

    push @commands, q{--ref-dir} . $SPACE . $dragen_hash_ref_dir_path;

    if ($alt_aware) {

        push @commands, q{--alt-aware} . $SPACE . q{true};
    }

    ## FASTQ  options
    if ($fastq_list_file_path) {

        push @commands, q{--fastq-list} . $SPACE . $fastq_list_file_path;
    }

    if ($fastq_list_all_samples) {

        push @commands, q{--fastq-list-all-samples} . $SPACE . q{true};
    }

    if ($fastq_list_sample_id) {

        push @commands, q{--fastq-list-sample-id} . $SPACE . $fastq_list_sample_id;
    }

    if ($fastq_infile_path) {

        push @commands, q{-1} . $SPACE . $fastq_infile_path;
    }
    if ($second_fastq_infile_path) {

        push @commands, q{-2} . $SPACE . $second_fastq_infile_path;
    }

    if ($is_fastq_interleaved) {

        push @commands, q{--interleaved};
    }

    ## BAM options
    if ($bam_file_path) {

        push @commands, q{-b} . $SPACE . $bam_file_path;
    }

    if ($enable_map_align) {

        push @commands, q{--enable-map-align} . $SPACE . q{true};
    }
    if ($enable_map_align_output) {

        push @commands, q{--enable-map-align-output} . $SPACE . q{true};
    }
    if ($enable_bam_indexing) {

        push @commands, q{--enable-bam-indexing} . $SPACE . q{true};
    }

    if ($enable_sampling) {

        push @commands, q{--enable-sampling} . $SPACE . q{true};
    }
    if ($combine_samples_by_name) {

        push @commands, q{--combine-samples-by-name} . $SPACE . q{true};
    }
    if ($enable_sort) {

        push @commands, q{--enable-sort} . $SPACE . q{true};
    }
    if ($enable_duplicate_marking) {

        push @commands, q{--enable-duplicate-marking} . $SPACE . q{true};
    }

    ## Variant calling
    if ($enable_variant_caller) {

        push @commands, q{--enable-variant-caller} . $SPACE . q{true};
    }

    if ($sample_id) {

        push @commands, q{--vc-sample-name} . $SPACE . $sample_id;
    }

    if ($vc_target_bed_file_path) {

        push @commands, q{--vc-target-bed} . $SPACE . $vc_target_bed_file_path;
    }

    if ($vc_enable_gatk_acceleration) {

        push @commands, q{--vc-enable-gatk-acceleration} . $SPACE . q{true};
    }

    if ($dbsnp_file_path) {

        push @commands, q{--dbsnp} . $SPACE . $dbsnp_file_path;
    }
    if ($vc_emit_ref_confidence) {

        push @commands, q{--vc-emit-ref-confidence} . $SPACE . q{GVCF};
    }
    if ($enable_combinegvcfs) {

        push @commands, q{--enable-combinegvcfs} . $SPACE . q{true};
    }
    if ( @{$sample_gvcf_file_paths_ref} ) {

        push @commands, q{--variant} . $SPACE . join $SPACE . q{--variant} . $SPACE,
          @{$sample_gvcf_file_paths_ref};
    }

    if ($enable_joint_genotyping) {

        push @commands, q{--enable-joint-genotyping} . $SPACE . q{true};
    }

    if ($enable_multi_sample_gvcf) {

        push @commands, q{--enable-multi-sample-gvcf} . $SPACE . q{true};
    }

    if ($disable_vcf_compression) {

        push @commands, q{--enable-vcf-compression} . $SPACE . q{false};
    }
    if ($pedigree_file_path) {

        push @commands, q{--vc-pedigree} . $SPACE . $pedigree_file_path;
    }
    ## CNV calling
    if ($enable_cnv) {

        push @commands, q{--enable-cnv} . $SPACE . q{true};
    }
    if ($cnv_enable_self_normalization) {

        push @commands, q{--cnv-enable-self-normalization} . $SPACE . q{true};
    }

    push @commands, q{--output-directory} . $SPACE . $outdirectory_path;

    push @commands, q{--output-file-prefix} . $SPACE . $outfile_prefix;

    if ($alignment_output_format) {

        push @commands, q{--output-format} . $SPACE . $alignment_output_format;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdinfile_path         => $stdinfile_path,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );
    return @commands;
}

sub dragen_build_hash_table {

## Function : Perl wrapper for builing a dragen hash table. Dragen version 3.3.5.
## Returns  : @commands
## Arguments: $build_hash_table           => Build hash table
##          : $enable_cnv                 => Enable cnv hash table creation
##          : $filehandle                 => Filehandle to write to
##          : $ht_alt_liftover_file_path  => Path to lift over file
##          : $ht_decoys_file_path        => Path to decoys file
##          : $outdirectory_path          => Outdirectory path
##          : $reference_genome_file_path => Reference genome file path [.fasta]
##          : $stderrfile_path            => Stderrfile path
##          : $stderrfile_path_append     => Append stderr info to file path
##          : $stdinfile_path             => Stdinfile path
##          : $stdoutfile_path            => Stdoutfile path
##          : $thread_number              => Number of threads

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $ht_alt_liftover_file_path;
    my $ht_decoys_file_path;
    my $outdirectory_path;
    my $reference_genome_file_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;

    ## Default(s)
    my $build_hash_table;
    my $enable_cnv;
    my $thread_number;

    my $tmpl = {
        build_hash_table => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$build_hash_table,
            strict_type => 1,
        },
        enable_cnv => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$enable_cnv,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        ht_alt_liftover_file_path => {
            store       => \$ht_alt_liftover_file_path,
            strict_type => 1,
        },
        ht_decoys_file_path => {
            store       => \$ht_decoys_file_path,
            strict_type => 1,
        },
        outdirectory_path => {
            defined     => 1,
            required    => 1,
            store       => \$outdirectory_path,
            strict_type => 1,
        },
        reference_genome_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$reference_genome_file_path,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
        stdinfile_path  => { store => \$stdinfile_path, strict_type => 1, },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        thread_number => {
            allow       => qr{ \A\d+\z }sxm,
            default     => 32,
            store       => \$thread_number,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ dragen };

    return if ( not $build_hash_table );

    push @commands, q{--build-hash-table} . $SPACE . q{true};

    push @commands, q{--ht-num-threads} . $SPACE . $thread_number;

    ## Should match $thread_number
    push @commands, q{--ht-max-table-chunks} . $SPACE . $thread_number;

    push @commands, q{--ht-reference} . $SPACE . $reference_genome_file_path;

    if ($ht_alt_liftover_file_path) {

        push @commands, q{--ht-alt-liftover} . $SPACE . $ht_alt_liftover_file_path;
    }
    if ($ht_decoys_file_path) {

        push @commands, q{--ht-decoys} . $SPACE . $ht_decoys_file_path;
    }

    if ($enable_cnv) {

        push @commands, q{--enable-cnv} . $SPACE . q{true};
    }
    push @commands, q{--output-directory} . $SPACE . $outdirectory_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdinfile_path         => $stdinfile_path,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );
    return @commands;
}

sub dragen_load_hash_table {

## Function : Perl wrapper for generic commands module.
## Returns  : @commands
## Arguments: $dragen_hash_ref_dir_path => Dragen reference genome dir path
##          : $filehandle               => Filehandle to write to
##          : $force_load_reference     => Force load dragen reference
##          : $stderrfile_path          => Stderrfile path
##          : $stderrfile_path_append   => Append stderr info to file path
##          : $stdinfile_path           => Stdinfile path
##          : $stdoutfile_path          => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $dragen_hash_ref_dir_path;
    my $filehandle;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;

    ## Default(s)
    my $force_load_reference;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        force_load_reference => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$force_load_reference,
            strict_type => 1,
        },
        dragen_hash_ref_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$dragen_hash_ref_dir_path,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
        stdinfile_path  => { store => \$stdinfile_path, strict_type => 1, },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ dragen };

    if ($force_load_reference) {

        push @commands, q{--force-load-reference};
    }

    push @commands, $dragen_hash_ref_dir_path;
    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdinfile_path         => $stdinfile_path,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
