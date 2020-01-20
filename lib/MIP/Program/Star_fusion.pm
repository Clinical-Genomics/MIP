package MIP::Program::Star_fusion;

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
    our $VERSION = 1.08;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ star_fusion star_fusion_gtf_file_to_feature_seqs star_fusion_prep_genome_lib };
}

Readonly my $READLENGTH => 150;

sub star_fusion {

## Function : Splice site analysis using STAR-fusion. Based on STAR-Fusion v1.8.0.
## Returns  :
## Arguments: $cpu                    => Number of threads for running STAR
##          : $examine_coding_effect  => Append coding effect to fusion
##          : $fastq_r1_path          => The path of the R1 fastq
##          : $fastq_r2_path          => The path of the R2 fastq
##          : $filehandle             => Filehandle to write to
##          : $fusion_inspector       => Run FusionInspector in either inspect or validate mode
##          : $genome_lib_dir_path    => Path to the directory containing the genome library
##          : $min_junction_reads     => Minimum number of reads spanning the junction {REF}
##          : $output_directory_path  => output directory path
##          : $samples_file_path      => Sample file path
##          : $sjdb_path              => Splice junction database file path (the junctions tab file)
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cpu;
    my $fastq_r1_path;
    my $fastq_r2_path;
    my $filehandle;
    my $fusion_inspector;
    my $genome_lib_dir_path;
    my $output_directory_path;
    my $min_junction_reads;
    my $samples_file_path;
    my $sjdb_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $examine_coding_effect;

    my $tmpl = {
        cpu => {
            store       => \$cpu,
            strict_type => 1,
        },
        examine_coding_effect => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$examine_coding_effect,
            strict_type => 1,
        },
        fastq_r1_path => {
            store       => \$fastq_r1_path,
            strict_type => 1,
        },
        fastq_r2_path => {
            store       => \$fastq_r2_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        fusion_inspector => {
            allow       => [ undef, qw{ inspect validate } ],
            store       => \$fusion_inspector,
            strict_type => 1,
        },
        genome_lib_dir_path => {
            required    => 1,
            store       => \$genome_lib_dir_path,
            strict_type => 1,
        },
        min_junction_reads => {
            allow       => [ undef, qr/ \A \d+ \z /xms ],
            store       => \$min_junction_reads,
            strict_type => 1,
        },
        output_directory_path => {
            required    => 1,
            store       => \$output_directory_path,
            strict_type => 1,
        },
        samples_file_path => {
            store       => \$samples_file_path,
            strict_type => 1,
        },
        sjdb_path => {
            store       => \$sjdb_path,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{STAR-Fusion};

    push @commands, q{--genome_lib_dir} . $SPACE . $genome_lib_dir_path;

    if ($samples_file_path) {

        push @commands, q{--samples_file} . $SPACE . $samples_file_path;
    }
    elsif ( $sjdb_path and not $fastq_r1_path ) {

        push @commands, q{-J} . $SPACE . $sjdb_path;
    }
    elsif ( $fastq_r1_path and $fastq_r2_path ) {

        push @commands, q{--right_fq} . $SPACE . $fastq_r1_path;
        push @commands, q{--left_fq} . $SPACE . $fastq_r2_path;
    }
    else {
        croak(
q{Error: You must either specify the fastq file paths or a splice junction database (SJDB) file.}
        );
    }

    if ($cpu) {
        push @commands, q{--CPU} . $SPACE . $cpu;
    }

    if ($examine_coding_effect) {
        push @commands, q{--examine_coding_effect};
    }

    if ($fusion_inspector) {
        push @commands, q{--FusionInspector} . $SPACE . $fusion_inspector,;
    }

    if ( defined $min_junction_reads ) {
        push @commands, q{--min_junction_reads} . $SPACE . $min_junction_reads;
    }

    push @commands, q{--output_dir} . $SPACE . $output_directory_path;

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

sub star_fusion_gtf_file_to_feature_seqs {

## Function : Perl wrapper for writing gtf_file_to_feature_seqs.pl (part of CTAT Genome Lib which is bundled with STAR-Fusion) command to $filehandle or return commands array. Based on STAR-Fusion v1.8.0
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $gtf_path               => Input gtf path
##          : $referencefile_path     => Reference sequence file
##          : $seq_type               => Sequence type
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $gtf_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $seq_type;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        gtf_path => {
            defined     => 1,
            required    => 1,
            store       => \$gtf_path,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        seq_type => {
            allow       => [qw{ cDNA CDS prot }],
            default     => q{cDNA},
            store       => \$seq_type,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ gtf_file_to_feature_seqs.pl };

    # Transcripts file
    push @commands, q{--gtf_file} . $SPACE . $gtf_path;

    # Reference sequence file
    push @commands, q{--genome_fa} . $SPACE . $referencefile_path;

    # Sequence type
    push @commands, q{--seqType} . $SPACE . $seq_type;

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

sub star_fusion_prep_genome_lib {

## Function : Perl wrapper for writing prep_genome_lib.pl (part of CTAT Genome Lib which is bundled with STAR-Fusion) command to $filehandle or return commands array. Based on STAR-Fusion v1.8.0
## Returns  : @commands
## Arguments: $dfam_db_path           => DNA transposable element database
##          : $filehandle             => Filehandle to write to
##          : $fusion_annot_lib_path  => Fusion annotation library
##          : $gtf_path               => Input gtf path
##          : $human_gencode_filter   => Customized prep operations for human/gencode genome and annotation data
##          : $output_dir_path        => Output directory path
##          : $pfam_db_path           => Pfam database
##          : $read_length            => Read length
##          : $referencefile_path     => Reference sequence file
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $tempdir_path           => Temporary directory
##          : $thread_number          => Number of threads (CPUs)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $fusion_annot_lib_path;
    my $gtf_path;
    my $output_dir_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $tempdir_path;

    ## Default(s)
    my $dfam_db_path;
    my $human_gencode_filter;
    my $pfam_db_path;
    my $read_length;
    my $thread_number;

    my $tmpl = {
        dfam_db_path => {
            default     => q{human},
            store       => \$dfam_db_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        fusion_annot_lib_path => {
            store       => \$fusion_annot_lib_path,
            strict_type => 1,
        },
        gtf_path => {
            defined     => 1,
            required    => 1,
            store       => \$gtf_path,
            strict_type => 1,
        },
        human_gencode_filter => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$human_gencode_filter,
            strict_type => 1,
        },
        output_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$output_dir_path,
            strict_type => 1,
        },
        pfam_db_path => {
            default     => $arg_href->{pfam_db_path} ||= q{current},
            store       => \$pfam_db_path,
            strict_type => 1,
        },
        read_length => {
            default     => $READLENGTH,
            store       => \$read_length,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
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
        tempdir_path => {
            store       => \$tempdir_path,
            strict_type => 1,
        },
        thread_number => {
            allow       => [ undef, qr/ \A \d+ \z /sxm ],
            default     => 16,
            store       => \$thread_number,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ prep_genome_lib.pl };

    push @commands, q{--dfam_db} . $SPACE . $dfam_db_path;

    if ($fusion_annot_lib_path) {

        push @commands, q{--fusion_annot_lib} . $SPACE . $fusion_annot_lib_path;
    }

    push @commands, q{--genome_fa} . $SPACE . $referencefile_path;

    push @commands, q{--gtf} . $SPACE . $gtf_path;

    if ($human_gencode_filter) {

        push @commands, q{--human_gencode_filter};
    }

    push @commands, q{--output_dir} . $SPACE . $output_dir_path;

    if ($pfam_db_path) {

        push @commands, q{--pfam_db} . $SPACE . $pfam_db_path;
    }

    if ($read_length) {

        push @commands, q{--max_readlength} . $SPACE . $read_length;
    }

    if ($tempdir_path) {

        push @commands, q{--outTmpDir} . $SPACE . $tempdir_path;
    }

    if ($thread_number) {

        push @commands, q{--CPU} . $SPACE . $thread_number;
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
