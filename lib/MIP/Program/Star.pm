package MIP::Program::Star;

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
use MIP::Constants qw{ $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.07;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ star_aln star_genome_generate };
}

sub star_aln {

## Function  : Perl wrapper for STAR v2.7.3a.
## Returns   : @commands
## Arguments : $align_intron_max              => Maximum intron size
##           : $align_mates_gap_max           => Maximum gap between two mates
##           : $align_sjdb_overhang_min       => Minimum overhang (i.e. block size) for spliced alignments
##           : $align_sj_stitch_mismatch_nmax => Number of mismatches allowed for each splicing motif (4)
##           : $chim_junction_overhang_min    => Minimum overhang for a chimeric junction
##           : $chim_score_junction_non_gtag  => Penalty for a non-GT/AG chimeric junction
##           : $chim_score_drop_max           => Max drop of total score from the read length
##           : $chim_score_min                => Minimum total score of chimeric alignment
##           : $chim_score_separation         => Minimum difference between chimeric scores
##           : $chim_segment_min              => Minimum length of chimeric segment
##           : $chim_segment_read_gap_max     => Maximum gap in the read sequence between chimeric segments
##           : $filehandle                    => Filehandle to write to
##           : $genome_dir_path               => Directory of the reference genome
##           : $infile_paths_ref              => Fastq file path(s)
##           : $limit_bam_sort_ram            => Memory available for sorting the output bam
##           : $outfile_name_prefix           => Prefix of the output files (remember to end with a ".")
##           : $out_bam_compression           => Compression level of BAM file
##           : $out_filter_mismatch_nmax      => Max number of missmatches allowed in alignment
##           : $out_filter_multimap_nmax      => Max number of loci a read is allowed to map to
##           : $out_sam_attr_rgline           => SAM/BAM read group line
##           : $out_sam_strand_field          => Cufflinks-like strand field flag
##           : $out_sam_type                  => Format of the output aligned reads
##           : $out_sam_unmapped              => Write unmapped reads to main BAM file
##           : $pe_overlap_nbases_min         => Min overlapp to trigger merging and realignment
##           : $quant_mode                    => Types of quantification requested
##           : $read_files_command            => A command which will be applied to the input files
##           : $stderrfile_path               => Stderrfile path
##           : $stderrfile_path_append        => Append stderr info to file path
##           : $stdout_data_type              => Which output will be directed to stdout
##           : $stdoutfile_path               => Stdoutfile path
##           : $thread_number                 => Number of threads
##           : $two_pass_mode                 => Two pass mode setting (None or Basic)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $align_intron_max;
    my $align_mates_gap_max;
    my $align_sjdb_overhang_min;
    my $align_sj_stitch_mismatch_nmax;
    my $chim_junction_overhang_min;
    my $chim_out_type;
    my $chim_score_drop_max;
    my $chim_score_junction_non_gtag;
    my $chim_score_min;
    my $chim_score_separation;
    my $chim_segment_min;
    my $filehandle;
    my $genome_dir_path;
    my $infile_paths_ref;
    my $limit_bam_sort_ram;
    my $outfile_name_prefix;
    my $out_bam_compression;
    my $out_filter_mismatch_nmax;
    my $out_filter_multimap_nmax;
    my $out_sam_attr_rgline;
    my $out_sam_unmapped;
    my $pe_overlap_nbases_min;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdout_data_type;
    my $stdoutfile_path;
    my $thread_number;

    ## Default(s)
    my $chim_segment_read_gap_max;
    my $out_sam_strand_field;
    my $out_sam_type;
    my $quant_mode;
    my $read_files_command;
    my $two_pass_mode;

    my $tmpl = {
        align_intron_max => {
            store       => \$align_intron_max,
            strict_type => 1,
        },
        align_mates_gap_max => {
            store       => \$align_mates_gap_max,
            strict_type => 1,
        },
        align_sjdb_overhang_min => {
            store       => \$align_sjdb_overhang_min,
            strict_type => 1,
        },
        align_sj_stitch_mismatch_nmax => {
            allow => qr{\A (?:  # open non-capture group
                        -?\d+\s # possible negative digit followed by space
                        )       # close grouping
                        {3}     # there is three of them
                        -?\d+   # final digit possibly negative
                        \z }xms,
            store       => \$align_sj_stitch_mismatch_nmax,
            strict_type => 1,
        },
        chim_junction_overhang_min => {
            allow       => [ undef, qr/\A \d+ \z /xms ],
            store       => \$chim_junction_overhang_min,
            strict_type => 1,
        },
        chim_out_type => {
            allow => [
                undef,
                qw{ Junctions SeparateSAMold WithinBAM },
                q{WithinBAM HardClip},
                q{WithinBAM SoftClip}
            ],
            store       => \$chim_out_type,
            strict_type => 1,
        },
        chim_score_drop_max => {
            allow       => qr/ \A \d+ \z /xms,
            store       => \$chim_score_drop_max,
            strict_type => 1,
        },
        chim_score_junction_non_gtag => {
            allow       => qr/ \A -? \d+ \z /xms,
            store       => \$chim_score_junction_non_gtag,
            strict_type => 1,
        },
        chim_score_min => {
            allow       => qr/ \A \d+ \z /xms,
            store       => \$chim_score_min,
            strict_type => 1,
        },
        chim_score_separation => {
            allow       => qr/ \A \d+ \z /xms,
            store       => \$chim_score_separation,
            strict_type => 1,
        },
        chim_segment_min => {
            store       => \$chim_segment_min,
            strict_type => 1,
        },
        chim_segment_read_gap_max => {
            allow       => qr/ \A \d+ \z /xms,
            default     => 3,
            store       => \$chim_segment_read_gap_max,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        genome_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$genome_dir_path,
            strict_type => 1,
        },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        limit_bam_sort_ram => {
            store       => \$limit_bam_sort_ram,
            strict_type => 1,
        },
        outfile_name_prefix => {
            defined     => 1,
            store       => \$outfile_name_prefix,
            strict_type => 1,
        },
        out_bam_compression => {
            allow       => qr/ -1 | [0-10] /xms,
            store       => \$out_bam_compression,
            strict_type => 1,
        },
        out_filter_mismatch_nmax => {
            allow       => qr/ \A \d+ \z /xms,
            store       => \$out_filter_mismatch_nmax,
            strict_type => 1,
        },
        out_filter_multimap_nmax => {
            allow       => qr/ \A \d+ \z /xms,
            store       => \$out_filter_multimap_nmax,
            strict_type => 1,
        },
        out_sam_attr_rgline => {
            store       => \$out_sam_attr_rgline,
            strict_type => 1,
        },
        out_sam_strand_field => {
            allow       => [qw{ None intronMotif }],
            default     => q{intronMotif},
            store       => \$out_sam_strand_field,
            strict_type => 1,
        },
        out_sam_type => {
            default     => q{BAM} . $SPACE . q{SortedByCoordinate},
            store       => \$out_sam_type,
            strict_type => 1,
        },
        out_sam_unmapped => {
            allow       => [ qw{ None Within }, q{Within KeepPairs} ],
            store       => \$out_sam_unmapped,
            strict_type => 1,
        },
        pe_overlap_nbases_min => {
            store       => \$pe_overlap_nbases_min,
            strict_type => 1,
        },
        quant_mode => {
            allow       => [qw{ - GeneCounts }],
            default     => q{GeneCounts},
            store       => \$quant_mode,
            strict_type => 1,
        },
        read_files_command => {
            default     => q{gunzip} . $SPACE . q{-c},
            store       => \$read_files_command,
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
        stdout_data_type => {
            allow       => [qw{ Log SAM BAM_Unsorted BAM_SortedByCoordinate BAM_Quant }],
            store       => \$stdout_data_type,
            strict_type => 1,
        },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        thread_number => {
            allow       => qr/ \A \d+ \z /xms,
            store       => \$thread_number,
            strict_type => 1,
        },
        two_pass_mode => {
            allow       => [qw{ Basic None }],
            default     => q{Basic},
            store       => \$two_pass_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ STAR };

    push @commands, q{--genomeDir} . $SPACE . $genome_dir_path;

    push @commands, q{--readFilesCommand} . $SPACE . $read_files_command;

    push @commands, q{--readFilesIn} . $SPACE . join $SPACE, @{$infile_paths_ref};

    push @commands, q{--outSAMtype} . $SPACE . $out_sam_type;

    if ($align_intron_max) {

        push @commands, q{--alignIntronMax} . $SPACE . $align_intron_max;
    }
    if ($align_mates_gap_max) {

        push @commands, q{--alignMatesGapMax} . $SPACE . $align_mates_gap_max;
    }
    if ($align_sjdb_overhang_min) {

        push @commands, q{--alignSJDBoverhangMin} . $SPACE . $align_sjdb_overhang_min;
    }
    if ($align_sj_stitch_mismatch_nmax) {

        push @commands,
          q{--alignSJstitchMismatchNmax} . $SPACE . $align_sj_stitch_mismatch_nmax;
    }
    if ( defined $chim_junction_overhang_min ) {

        push @commands,
          q{--chimJunctionOverhangMin} . $SPACE . $chim_junction_overhang_min;
    }
    if ($chim_out_type) {

        push @commands, q{--chimOutType} . $SPACE . $chim_out_type;
    }
    if ( defined $chim_score_drop_max ) {

        push @commands, q{--chimScoreDropMax} . $SPACE . $chim_score_drop_max;
    }
    if ( defined $chim_score_junction_non_gtag ) {

        push @commands,
          q{--chimScoreJunctionNonGTAG} . $SPACE . $chim_score_junction_non_gtag;
    }
    if ( defined $chim_score_min ) {

        push @commands, q{--chimScoreMin} . $SPACE . $chim_score_min;
    }
    if ( defined $chim_score_separation ) {

        push @commands, q{--chimScoreSeparation} . $SPACE . $chim_score_separation;
    }
    if ( defined $chim_segment_min ) {

        push @commands, q{--chimSegmentMin} . $SPACE . $chim_segment_min;
    }
    if ( defined $chim_segment_read_gap_max ) {

        push @commands, q{--chimSegmentReadGapMax} . $SPACE . $chim_segment_read_gap_max;
    }
    if ($limit_bam_sort_ram) {

        push @commands, q{--limitBAMsortRAM} . $SPACE . $limit_bam_sort_ram;
    }
    if ( defined $out_bam_compression ) {

        push @commands, q{--outBAMcompression} . $SPACE . $out_bam_compression;
    }
    if ($outfile_name_prefix) {

        push @commands, q{--outFileNamePrefix} . $SPACE . $outfile_name_prefix;
    }
    if ($out_filter_mismatch_nmax) {

        push @commands, q{--outFilterMismatchNmax} . $SPACE . $out_filter_mismatch_nmax;
    }
    if ($out_filter_multimap_nmax) {

        push @commands, q{--outFilterMultimapNmax} . $SPACE . $out_filter_multimap_nmax;
    }
    if ($out_sam_attr_rgline) {

        push @commands, q{--outSAMattrRGline} . $SPACE . $out_sam_attr_rgline;
    }
    if ($out_sam_strand_field) {

        push @commands, q{--outSAMstrandField} . $SPACE . $out_sam_strand_field;
    }
    if ($out_sam_unmapped) {

        push @commands, q{--outSAMunmapped} . $SPACE . $out_sam_unmapped;
    }
    if ($pe_overlap_nbases_min) {

        push @commands, q{--peOverlapNbasesMin} . $SPACE . $pe_overlap_nbases_min;
    }
    if ($quant_mode) {

        push @commands, q{--quantMode} . $SPACE . $quant_mode;
    }
    if ($stdout_data_type) {

        push @commands, q{--outStd} . $SPACE . $stdout_data_type;
    }
    if ($thread_number) {

        push @commands, q{--runThreadN} . $SPACE . $thread_number;
    }
    if ($two_pass_mode) {

        push @commands, q{--twopassMode} . $SPACE . $two_pass_mode;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
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

sub star_genome_generate {

## Function  : Perl wrapper for STAR genomeGenerate, v2.5.3a.
## Returns   : @commands
## Arguments : $fasta_path             => Input reference fasta path
##           : $filehandle             => Filehandle to write to
##           : $genome_dir_path        => Output directory path
##           : $gtf_path               => Input gtf path
##           : $read_length            => Maximum expected readlength
##           : $stderrfile_path        => Stderrfile path
##           : $stderrfile_path_append => Append stderr info to file path
##           : $stdoutfile_path        => Stdoutfile path
##           : $thread_number          => Number of threads

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fasta_path;
    my $filehandle;
    my $genome_dir_path;
    my $gtf_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $read_length;
    my $thread_number;

    my $tmpl = {
        fasta_path => {
            defined     => 1,
            required    => 1,
            store       => \$fasta_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        genome_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$genome_dir_path,
            strict_type => 1,
        },
        gtf_path => {
            defined     => 1,
            required    => 1,
            store       => \$gtf_path,
            strict_type => 1,
        },
        read_length => {
            allow       => qr/ ^\d+$ /xms,
            default     => 150,
            store       => \$read_length,
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
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        thread_number => {
            allow       => qr/ ^\d+$ /xms,
            default     => 36,
            store       => \$thread_number,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ STAR --runMode genomeGenerate };

    # Options
    push @commands, q{--genomeFastaFiles} . $SPACE . $fasta_path;

    push @commands, q{--genomeDir} . $SPACE . $genome_dir_path;

    push @commands, q{--sjdbGTFfile} . $SPACE . $gtf_path;

    if ($read_length) {
        push @commands, q{--sjdbOverhang} . $SPACE . $read_length;

    }
    if ($thread_number) {
        push @commands, q{--runThreadN} . $SPACE . $thread_number;

    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
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
