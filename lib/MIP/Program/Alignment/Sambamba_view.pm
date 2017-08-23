package MIP::Program::Alignment::Sambamba_view;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use Carp;
use utf8;    #Allow unicode characters in this script
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };

use FindBin qw{$Bin};    #Find directory of script
use File::Basename qw{dirname};
use File::Spec::Functions qw{catdir};

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Unix::Standard_streams qw{unix_standard_streams};
use MIP::Unix::Write_to_file qw{unix_write_to_file};

BEGIN {
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Inherit from Exporter to export functions and variables
    use base qw {Exporter};

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{view};

}

use Params::Check qw{check allow last_error};
use Readonly;
Readonly my $SPACE => q{ };

sub view {

##view

##Function : Perl wrapper for writing sambamba view recipe to $FILEHANDLE. Based on sambamba 0.6.5
##Returns  : "@commands"
##Arguments: $regions_ref, $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $with_header, $show_progress, $output_format, $referencefile_path
##         : $regions_ref        => The regions to process {REF}
##         : $FILEHANDLE         => Sbatch filehandle to write to
##         : $infile_path        => Infile path
##         : $outfile_path       => Outfile path
##         : $stderrfile_path    => Stderrfile path
##         : $referencefile_path => Reference for writing CRAM
##         : $with_header        => Include header
##         : $show_progress      => Show progress
##         : $output_format      => Output format

    my ($arg_href) = @_;

    ## Default(s)
    my $with_header;
    my $show_progress;
    my $output_format;

    ## Flatten argument(s)
    my $regions_ref;
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $referencefile_path;
    my $stderrfile_path_append;

    my $tmpl = {
        regions_ref =>
          { default => [], strict_type => 1, store => \$regions_ref },
        FILEHANDLE  => { required => 1,  store => \$FILEHANDLE },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        referencefile_path =>
          { strict_type => 1, store => \$referencefile_path },
        with_header => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$with_header
        },
        show_progress => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$show_progress
        },
        output_format => {
            default     => qw{bam},
            allow       => [qw{sam bam cram json}],
            strict_type => 1,
            store       => \$output_format
        },
    };

    check( $tmpl, $arg_href, 1 )
      or croak qw{Could not parse arguments!};

    ## Sambamba
    my @commands =
      qw{sambamba view};    #Stores commands depending on input parameters

    if ($with_header) {     #Include header

        push @commands, q{--with-header};
    }

    if ($output_format) {

        push @commands, q{--format} . $SPACE . $output_format;    #Output format
    }

    if ($referencefile_path) {

        push @commands,
          q{--ref-filename=} . $referencefile_path;  #Reference for writing CRAM
    }

    if ($show_progress) {

        push @commands, q{--show-progress}
          ; #Show progressbar in STDERR (works only for BAM files with no regions specified)
    }

    if ($outfile_path) {

        push @commands,
          q{--output-filename=} . $outfile_path;    #Specify output filename
    }

    ## Infile
    push @commands, $infile_path;

    if (@$regions_ref) {                            #Limit output to regions

        push @commands, join( $SPACE, @{$regions_ref} );
    }

    if ($stderrfile_path) {

        push @commands,
          unix_standard_streams(
            {
                stderrfile_path        => $stderrfile_path,
                stderrfile_path_append => $stderrfile_path_append,
            }
          );
    }
    if ($FILEHANDLE) {

        unix_write_to_file(
            {
                commands_ref => \@commands,
                separator    => $SPACE,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
    }
    return @commands;
}

1;
