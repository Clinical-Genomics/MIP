package MIP::Unix::Standard_streams;

use strict;
use warnings;
use warnings qw(FATAL utf8);
use utf8;    #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use autodie;

BEGIN {

    use base qw(Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(unix_standard_streams);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

sub unix_standard_streams {

##unix_standard_streams

##Function : Perl wrapper for writing unix standard_streams recipe to already open $FILEHANDLE or return commands array.
##Returns  : "@commands"
##Arguments: $stdoutfile_path, $stderrfile_path, stderrfile_path_append, $FILEHANDLE
##         : $stdoutfile_path        => Directory path
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append stderr info to file
##         : $FILEHANDLE             => Filehandle to write to

    my ($arg_href) = @_;

    ## Default(s)
    my $stderrfile_path_append;

    ## Flatten argument(s)
    my $stdoutfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = {
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE => { store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    my $SPACE = q{ };

    ### Standard streams helper module
    ##Stores commands depending on input parameters
    my @commands;

    ## Options
    if ($stdoutfile_path) {

        # Redirect stdout to program specific stdout file
        push @commands, '1> ' . $stdoutfile_path;
    }
    if ($stderrfile_path) {

        # Redirect stderr output to program specific stderr file
        push @commands, '2> ' . $stderrfile_path;
    }
    if ($stderrfile_path_append) {

        # Redirect and append stderr output to program specific stderr file
        push @commands, '2>> ' . $stderrfile_path_append;
    }
    if ($FILEHANDLE) {

        print {$FILEHANDLE} join( $SPACE, @commands ) . $SPACE;
    }
    return @commands;
}

1;
