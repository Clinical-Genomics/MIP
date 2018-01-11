package MIP::Program::Alignment::Star;

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
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ space separated subroutines };
}

## Constants
Readonly my $SPACE => q{ };

sub star_aln {

## Function : Perl wrapper for STAR.
#J35P312 was here!
## Returns  : @commands
## Arguments: $FILEHANDLE                 => Filehandle to write to
##          : $stderrfile_path            => Stderrfile path
##          : $stderrfile_path_append     => Append stderr info to file path
##          : $stdoutfile_path            => Stdoutfile path
##          : $thread_number              => Number of threads
##          : $first_infile_path          => first infile path (read 1)
##          : $second_infile_path         => Second infile path (read 2)
##          : $genomeDir                  => Directory of the reference genome
##          : $twopassMode                => two pass mode setting (None or Basic)
##          : $chimSegmentMin             => minimum length of chimaeric segment
##          : $chimJunctionOverhangMin    => minimum overhang for a chimeric junction
##          : $chimSegmentReadGapMax      => maximum gap in the read sequence between chimeric segments
##          : $alignSJDBoverhangMin       => minimum overhang (i.e. block size) for spliced alignments
##          : $alignMatesGapMax           => maximum gap between two mates
##          : $alignIntronMax             => maximum intron size
##          : $limitBAMsortRAM            => memory available for sorting the output bam
##          : $outFileNamePrefix          => prefix of the output files (remember to end with a ".")
##          : $quantMode                  => types of quantification requested
##          : $outSAMstrandField          => Cufflinks-like strand field flag

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## STAR arguments
    my $second_infile_path;
    my $first_infile_path;
    my $thread_number;
    my $genomeDir;
    my $twopassMode;
    my $chimSegmentMin;
    my $chimJunctionOverhangMin;
    my $chimSegmentReadGapMax; 
    my $alignSJDBoverhangMin;
    my $alignMatesGapMax;
    my $alignIntronMax;
    my $limitBAMsortRAM;
    my $outFileNamePrefix;
    my $quantMode;
    my $outSAMstrandField;

    ## Default(s)
    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
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
        first_infile_path => {
            required     => 1,
            defined      => 1,
            strict_type  => 1,
            store        => \$first_infile_path
        },
        second_infile_path => {
            strict_type  => 1,
            store        => \$second_infile_path
        },
        thread_number => {
            allow       => qr/ ^\d+$ /xms,
            strict_type => 1,
            default     => 16,
            store       => \$thread_number
        },
        genomeDir => {
            required     => 1,
            defined      => 1,
            strict_type  => 1,
            store        => \$genomeDir
        },
        twopassMode => {
            allow       => [qw{ Basic None }],
            default      => q{Basic},
            strict_type  => 1,
            store        => \$twopassMode
        },
        chimSegmentMin => {
            default      => 12,
            strict_type  => 1,
            store        => \$chimSegmentMin
        },
        chimJunctionOverhangMin => {
            default      => 12,
            strict_type  => 1,
            store        => \$chimJunctionOverhangMin
        },
        chimSegmentReadGapMax => {
            default      => 3,
            strict_type  => 1,
            store        => \$chimSegmentReadGapMax
        },
        alignSJDBoverhangMin => {
            default      => 10,
            strict_type  => 1,
            store        => \$alignSJDBoverhangMin
        },
        alignMatesGapMax => {
            default      => 100000,
            strict_type  => 1,
            store        => \$alignMatesGapMax
        },
        alignIntronMax   => {
            default      => 100000,
            strict_type  => 1,
            store        => \$alignIntronMax
        },
        limitBAMsortRAM   => {
            default      => 31532137230,
            strict_type  => 1,
            store        => \$limitBAMsortRAM
        },
        outFileNamePrefix   => {
            required     => 1,
            defined      => 1,
            strict_type  => 1,
            store        => \$outFileNamePrefix
        },
        quantMode => {
            allow       => [qw{ - GeneCounts }],
            default      => q{GeneCounts},
            strict_type  => 1,
            store        => \$quantMode
        },
        outSAMstrandField => {
            allow       => [qw{ None intronMotif }],
            default      => q{intronMotif},
            strict_type  => 1,
            store        => \$outSAMstrandField
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{STAR};

    #required input
    if ($genomeDir) {
        push @commands, q{--genomeDir} . $SPACE . $genomeDir;

    }if ($first_infile_path){
        if ($second_infile_path){
            push @commands, q{--readFilesIn} . $SPACE . $first_infile_path . $SPACE . $second_infile_path;
        }else{
            push @commands, q{--readFilesIn} . $SPACE . $first_infile_path;    
        }
    }

    #Assume gzipped fastqs
    push @commands, q{--readFilesCommand gzip -c};
    #We only want sorted bam output
    push @commands, q{--outSAMtype BAM SortedByCoordinate};


    #Options
    if ($thread_number) {
        push @commands, q{--runThreadN} . $SPACE . $thread_number;

    }if ($twopassMode){
        push @commands, q{--twopassMode} . $SPACE . $twopassMode;

    }if ($chimSegmentMin){
        push @commands, q{--chimSegmentMin} . $SPACE . $chimSegmentMin;

    }if ($chimJunctionOverhangMin){
        push @commands, q{--chimJunctionOverhangMin} . $SPACE . $chimJunctionOverhangMin;

    }if ($chimSegmentReadGapMax){
        push @commands, q{--chimSegmentReadGapMax} . $SPACE . $chimSegmentReadGapMax;

    }if ($alignSJDBoverhangMin){
        push @commands, q{--alignSJDBoverhangMin} . $SPACE . $alignSJDBoverhangMin;

    }if ($alignMatesGapMax){
        push @commands, q{--alignMatesGapMax} . $SPACE . $alignMatesGapMax;

    }if ($alignIntronMax){
        push @commands, q{--alignIntronMax} . $SPACE . $alignIntronMax;

    }if ($limitBAMsortRAM){
        push @commands, q{--limitBAMsortRAM} . $SPACE . $limitBAMsortRAM;

    }if ($outFileNamePrefix){
        push @commands, q{--outFileNamePrefix} . $SPACE . $outFileNamePrefix;

    }if ($quantMode){
        push @commands, q{--quantMode} . $SPACE . $quantMode;

    }if ($outSAMstrandField){
        push @commands, q{--outSAMstrandField} . $SPACE . $outSAMstrandField;
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
            FILEHANDLE   => $FILEHANDLE,
            commands_ref => \@commands,
            separator    => $SPACE,

        }
    );
    return @commands;
}

sub star_index{

## Function : Perl wrapper for STAR genomeGenerate.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $thread_number          => Number of threads
##          : $genomeDir              => output directory
##          : $fasta                  => input reference fasta
##          : $gtf                    => input gtf
##          : $readlen                => maximum expected readlength

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $thread_number;
    my $genomeDir;
    my $fasta;
    my $gtf;
    my $readlen;

    ## Default(s)

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
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
            strict_type => 1,
            default     => 16,
            store       => \$thread_number
        },
        genomeDir => {
            required     => 1,
            defined      => 1,
            strict_type  => 1,
            store        => \$genomeDir
        },
        fasta => {
            required     => 1,
            defined      => 1,
            strict_type  => 1,
            store        => \$fasta
        },
        gtf => {
            required     => 1,
            defined      => 1,
            strict_type  => 1,
            store        => \$gtf
        },
        readlen => {
            allow       => qr/ ^\d+$ /xms,
            strict_type => 1,
            default     => 150,
            store       => \$readlen
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{STAR --runMode genomeGenerate};

    #options
    if ($thread_number) {
        push @commands, q{--runThreadN} . $SPACE . $thread_number;

    }if ($genomeDir) {
        push @commands, q{--genomeDir} . $SPACE . $genomeDir;

    }if ($fasta) {
        push @commands, q{--genomeFastaFiles} . $SPACE . $fasta;

    }if ($gtf) {
        push @commands, q{--sjdbGTFfile} . $SPACE . $gtf;

    }if ($readlen) {
        push @commands, q{--sjdbOverhang} . $SPACE . $readlen;

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
            FILEHANDLE   => $FILEHANDLE,
            commands_ref => \@commands,
            separator    => $SPACE,

        }
    );
    return @commands;

}

1;
