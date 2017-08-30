package MIP::Program::Alignment::Bedtools;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use Carp;
use utf8;    #Allow unicode characters in this script
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };

use FindBin qw{$Bin};    # Find directory of script
use File::Basename qw{dirname};
use File::Spec::Functions qw{catdir};

BEGIN {
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Inherit from Exporter to export functions and variables
    use base qw {Exporter};

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{bedtools_genomecov};

}

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Unix::Standard_streams qw{unix_standard_streams};
use MIP::Unix::Write_to_file qw{unix_write_to_file};

use Params::Check qw{check allow last_error};
use Readonly;

## Constants
Readonly my $SPACE => q{ };

sub bedtools_genomecov {

## bedtools_genomecov

## Function : Perl wrapper for writing bedtools genomecov recipe to $FILEHANDLE. Based on bedtools 2.26.0.
## Returns  : "@commands"
## Arguments: $infile_path, $outfile_path, $referencefile_path, $stderrfile_path, $FILEHANDLE, $max_coverage
##          : $infile_path        => Infile paths
##          : $outfile_path       => Outfile path
##          : $referencefile_path => Genome reference file
##          : $stderrfile_path    => Stderrfile path
##          : $FILEHANDLE         => Sbatch filehandle to write to
##          : $max_coverage       => Combine all positions with a depth >= max into a single bin in the histogram

    my ($arg_href) = @_;

    ## Default(s)
    my $max_coverage;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $stderrfile_path_append;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path       => { strict_type => 1, store => \$outfile_path },
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE   => { store => \$FILEHANDLE },
        max_coverage => {
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$max_coverage
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ## Array @commands stores commands depending on input parameters
    my @commands = qw{bedtools genomecov};

    ## Options
    if ( defined $max_coverage ) {

        push @commands, q{-max} . $SPACE . $max_coverage;
    }

    ## Infile
    push @commands, q{-ibam} . $SPACE . $infile_path;

    if ($referencefile_path) {

        push @commands, q{-g} . $SPACE . $referencefile_path;
    }

    ## Output
    if ($outfile_path) {

        # Specify output filename
        push @commands, q{>} . $SPACE . $outfile_path;
    }

    # Redirect stderr output to program specific stderr file
    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

1;
