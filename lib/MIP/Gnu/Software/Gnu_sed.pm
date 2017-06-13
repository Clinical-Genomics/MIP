package MIP::Gnu::Software::Gnu_sed;

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;    #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use autodie;
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

use FindBin qw($Bin);                 #Find directory of script
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir);

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Unix::Standard_streams qw(unix_standard_streams);
use MIP::Unix::Write_to_file qw(unix_write_to_file);

BEGIN {

    use base qw(Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(gnu_sed);

}

sub gnu_sed {

##gnu_sed

##Function : Perl wrapper for writing sed recipe to already open $FILEHANDLE or return commands array. Based on sed 4.2.1.
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $stderrfile_path_append, $script
##         : $FILEHANDLE       => Filehandle to write to
##         : $infile_path      => Infile path
##         : $outfile_path     => Outfile path
##         : $stderrfile_path_append => Append stderr info to file
##         : $stderrfile_path  => Stderrfile path
##         : $script           => Script to edit infile stream

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $script;
    my $stderrfile_path_append;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE
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
        script => {
            strict_type => 1,
            store       => \$script
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $SPACE = q{ };

    ### sed
    ##Stores commands depending on input parameters
    my @commands = qw(sed);

    ## Options
    if ($script) {

        push @commands, $script;
    }

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
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

1;
