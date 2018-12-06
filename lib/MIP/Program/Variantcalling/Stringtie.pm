package MIP::Program::Variantcalling::Stringtie;

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
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ stringtie };
}

## Constants
Readonly my $SPACE => q{ };

sub stringtie {

## Function : Perl wrapper for StringTie. Based on version 1.3.4.
## Returns  : @commands
## Arguments: $cov_ref_transcripts_outfile_path => Fully covered reference transcripts
##          : $FILEHANDLE                       => Filehandle to write to
##          : $gene_abundance_outfile_path      => Gene aboundances
##          : $gtf_reference_path               => Input GTF refrence file
##          : $infile_path                      => Input bam file path
##          : $library_type                     => Orientation and strandedness of the library
##          : $outfile_path                     => Path to output GTF
##          : $stderrfile_path                  => Stderrfile path
##          : $stderrfile_path_append           => Append stderr info to file path
##          : $stdoutfile_path                  => Stdoutfile path
##          : $threads                          => Number of threads

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cov_ref_transcripts_outfile_path;
    my $FILEHANDLE;
    my $gene_abundance_outfile_path;
    my $gtf_reference_path;
    my $infile_path;
    my $library_type;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $threads;

    my $tmpl = {
        cov_ref_transcripts_outfile_path => {
            store       => \$cov_ref_transcripts_outfile_path,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        gene_abundance_outfile_path => {
            store       => \$gene_abundance_outfile_path,
            strict_type => 1,
        },
        gtf_reference_path => {
            defined     => 1,
            required    => 1,
            store       => \$gtf_reference_path,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        library_type => {
            defined     => 1,
            store       => \$library_type,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
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
        threads => {
            store       => \$threads,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{stringtie};

    push @commands, $infile_path;

    if ($cov_ref_transcripts_outfile_path) {
        push @commands, q{-C} . $SPACE . $cov_ref_transcripts_outfile_path;
    }

    if ($gene_abundance_outfile_path) {
        push @commands, q{-A} . $SPACE . $gene_abundance_outfile_path;
    }

    push @commands, q{-G} . $SPACE . $gtf_reference_path;

    if ( $library_type and $library_type eq q{forward_stranded} ) {
        push @commands, q{--fr};
    }
    elsif ( $library_type and $library_type eq q{reverse_stranded} ) {
        push @commands, q{--rf};
    }

    if ($threads) {
        push @commands, q{-p} . $SPACE . $threads;
    }

    push @commands, q{-o} . $SPACE . $outfile_path;

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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,

        }
    );
    return @commands;

}

1;
