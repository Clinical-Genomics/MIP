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
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ star_aln star_genome_generate };
}

sub star_aln {

## Function  : Perl wrapper for STAR v2.5.3a.
## Returns   : @commands
## Arguments : $align_intron_max           => Maximum intron size
##           : $align_mates_gap_max        => Maximum gap between two mates
##           : $align_sjdb_overhang_min    => Minimum overhang (i.e. block size) for spliced alignments
##           : $chim_junction_overhang_min => Minimum overhang for a chimeric junction
##           : $chim_segment_min           => Minimum length of chimaeric segment
##           : $chim_segment_read_gap_max  => Maximum gap in the read sequence between chimeric segments
##           : $filehandle                 => Filehandle to write to
##           : $genome_dir_path            => Directory of the reference genome
##           : $infile_paths_ref           => Fastq file path(s)
##           : $limit_bam_sort_ram         => Memory available for sorting the output bam
##           : $outfile_name_prefix        => Prefix of the output files (remember to end with a ".")
##           : $out_sam_strand_field       => Cufflinks-like strand field flag
##           : $out_sam_type               => Format of the output aligned reads
##           : $pe_overlap_nbases_min      => Min overlapp to trigger merging and realignment
##           : $quant_mode                 => Types of quantification requested
##           : $read_files_command         => A command which will be applied to the input files
##           : $stderrfile_path            => Stderrfile path
##           : $stderrfile_path_append     => Append stderr info to file path
##           : $stdoutfile_path            => Stdoutfile path
##           : $thread_number              => Number of threads
##           : $two_pass_mode              => Two pass mode setting (None or Basic)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $chim_out_type;
    my $genome_dir_path;
    my $infile_paths_ref;
    my $outfile_name_prefix;
    my $pe_overlap_nbases_min;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $align_intron_max;
    my $align_mates_gap_max;
    my $align_sjdb_overhang_min;
    my $chim_junction_overhang_min;
    my $chim_segment_min;
    my $chim_segment_read_gap_max;
    my $limit_bam_sort_ram;
    my $out_sam_strand_field;
    my $out_sam_type;
    my $quant_mode;
    my $read_files_command;
    my $thread_number;
    my $two_pass_mode;

    my $tmpl = {
        align_intron_max => {
            default     => 100_000,
            store       => \$align_intron_max,
            strict_type => 1,
        },
        align_mates_gap_max => {
            default     => 100_000,
            store       => \$align_mates_gap_max,
            strict_type => 1,
        },
        align_sjdb_overhang_min => {
            default     => 10,
            store       => \$align_sjdb_overhang_min,
            strict_type => 1,
        },
        chim_junction_overhang_min => {
            default     => 12,
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
        chim_segment_min => {
            default     => 12,
            store       => \$chim_segment_min,
            strict_type => 1,
        },
        chim_segment_read_gap_max => {
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
            default     => 315_321_372_30,
            store       => \$limit_bam_sort_ram,
            strict_type => 1,
        },
        outfile_name_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_name_prefix,
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
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        thread_number => {
            allow       => qr/ ^\d+$ /xms,
            default     => 16,
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

    ## Stores commands depending on input parameters
    my @commands = qw{ STAR };

    push @commands, q{--genomeDir} . $SPACE . $genome_dir_path;

    push @commands, q{--readFilesCommand} . $SPACE . $read_files_command;

    push @commands, q{--readFilesIn} . $SPACE . join $SPACE, @{$infile_paths_ref};

    push @commands, q{--outFileNamePrefix} . $SPACE . $outfile_name_prefix;

    push @commands, q{--outSAMtype} . $SPACE . $out_sam_type;

    ## Options
    if ($align_intron_max) {
        push @commands, q{--alignIntronMax} . $SPACE . $align_intron_max;

    }
    if ($align_mates_gap_max) {
        push @commands, q{--alignMatesGapMax} . $SPACE . $align_mates_gap_max;

    }
    if ($align_sjdb_overhang_min) {
        push @commands, q{--alignSJDBoverhangMin} . $SPACE . $align_sjdb_overhang_min;

    }
    if ($chim_segment_min) {
        push @commands, q{--chimSegmentMin} . $SPACE . $chim_segment_min;

    }
    if ($chim_junction_overhang_min) {
        push @commands,
          q{--chimJunctionOverhangMin} . $SPACE . $chim_junction_overhang_min;

    }
    if ($chim_out_type) {
        push @commands, q{--chimOutType} . $SPACE . $chim_out_type;

    }
    if ($chim_segment_read_gap_max) {
        push @commands, q{--chimSegmentReadGapMax} . $SPACE . $chim_segment_read_gap_max;

    }
    if ($limit_bam_sort_ram) {
        push @commands, q{--limitBAMsortRAM} . $SPACE . $limit_bam_sort_ram;

    }
    if ($pe_overlap_nbases_min) {
        push @commands, q{--peOverlapNbasesMin} . $SPACE . $pe_overlap_nbases_min;

    }
    if ($quant_mode) {
        push @commands, q{--quantMode} . $SPACE . $quant_mode;

    }
    if ($out_sam_strand_field) {
        push @commands, q{--outSAMstrandField} . $SPACE . $out_sam_strand_field;
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
            default     => 16,
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
