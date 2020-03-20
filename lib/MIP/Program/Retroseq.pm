package MIP::Program::Retroseq;

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
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ retroseq_call retroseq_discover };
}

sub retroseq_call {

## Function : Perl wrapper for Retroseq version 1.5 mobile element detector
## Returns  : @commands
## Arguments: $filehandle              => Filehandle to write to
##          : $infile_path             => Path to input bam file
##          : $outputfile_path         => Path to the output file
##          : $reference_fasta_path    => Reference genome path
##          : $retroseq_bed_path       => Path to the retroseq discover bedfile
##          : $stderrfile_path         => Stderrfile path
##          : $stderrfile_path_append  => Append stderr info to file path
##          : $stdoutfile_path         => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outputfile_path;
    my $reference_fasta_path;
    my $retroseq_bed_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outputfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outputfile_path,
            strict_type => 1,
        },
        reference_fasta_path => {
            defined     => 1,
            required    => 1,
            store       => \$reference_fasta_path,
            strict_type => 1,
        },
        retroseq_bed_path => {
            defined     => 1,
            required    => 1,
            store       => \$retroseq_bed_path,
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
    my @commands = qw{ retroseq.pl };
    push @commands, q{-call -soft};
    push @commands, q{-bam} . $SPACE . $infile_path;

    push @commands, q{-input} . $SPACE . $retroseq_bed_path;

    push @commands, q{-ref} . $SPACE . $reference_fasta_path;

    push @commands, q{-output} . $SPACE . $outputfile_path;

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

sub retroseq_discover {

## Function : Perl wrapper for Retroseq version 1.5 mobile element detector
## Returns  : @commands
## Arguments: $filehandle              => Filehandle to write to
##          : $infile_path             => Path to input bam file
##          : $mobile_element_tsv_path => Tab separated file containing the name and path of mobile elements
##          : $outputfile_path         => path to the output file
##          : $stderrfile_path         => Stderrfile path
##          : $stderrfile_path_append  => Append stderr info to file path
##          : $stdoutfile_path         => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $mobile_element_tsv_path;
    my $outputfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        mobile_element_tsv_path => {
            defined     => 1,
            required    => 1,
            store       => \$mobile_element_tsv_path,
            strict_type => 1,
        },
        outputfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outputfile_path,
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
    my @commands = qw{ retroseq.pl };
    push @commands, q{-discover};

    push @commands, q{-bam} . $SPACE . $infile_path;

    push @commands, q{-refTEs} . $SPACE . $mobile_element_tsv_path;

    push @commands, q{-output} . $SPACE . $outputfile_path;

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
