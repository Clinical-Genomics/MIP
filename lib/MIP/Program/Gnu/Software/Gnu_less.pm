package MIP::Program::Gnu::Software::Gnu_less;

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;    #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use autodie;
use Params::Check qw[check allow last_error];

use FindBin qw($Bin);                 #Find directory of script
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir);

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Unix::Standard_streams qw(unix_standard_streams);
use MIP::Unix::Write_to_file qw(unix_write_to_file);

BEGIN {
    use base qw (Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(gnu_less);

}

sub gnu_less {

##gnu_less

##Function : Perl wrapper for writing less recipe to already open $filehandle or return commands array. Based on less 436
##Returns  : "@commands"
##Arguments: $filehandle, $infile_path, $outfile_path, $stderrfile_path, $stderrfile_path_append
##         : $filehandle                => Filehandle to write to
##         : $infile_path               => Infile path
##         : $outfile_path              => Outfile path
##         : $stderrfile_path           => Stderrfile path
##         : $stderrfile_path_append    => Append stderrinfo to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    my $tmpl = {
        filehandle => {
            store => \$filehandle
        },
        infile_path => {
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path => {
            strict_type => 1,
            store       => \$outfile_path
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $SPACE = q{ };

    ## less
    my @commands = qw(less);    #Stores commands depending on input parameters

    ## Options

    ## Infile
    if ($infile_path) {
        push @commands, $infile_path;
    }

    ## Outfile
    if ($outfile_path) {
        push @commands, '> ' . $outfile_path;
    }

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
            filehandle   => $filehandle,
        }
    );
    return @commands;
}

1;
