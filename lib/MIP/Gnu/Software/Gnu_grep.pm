package MIP::Gnu::Software::Gnu_grep;

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;    #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use autodie;

BEGIN {

    use base qw(Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = '1.00';

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(gnu_grep);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case
sub gnu_grep {

##gnu_grep

##Function : Perl wrapper for writing grep recipe to already open $FILEHANDLE or return commands array. Based on grep 2.6.3
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $filter_file_path, $invert_match
##         : $FILEHANDLE       => Filehandle to write to
##         : $infile_path      => Infile path
##         : $outfile_path     => Outfile path
##         : $stderrfile_path  => Stderrfile path
##         : $filter_file_path => Obtain patterns from file, one per line
##         : $invert_match     => Invert the sense of matching, to select non-matching lines

    my ($arg_href) = @_;

    ## Default(s)
    my $invert_match;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $filter_file_path;

    my $tmpl = {
        infile_path      => { strict_type => 1, store => \$infile_path },
        FILEHANDLE       => { store       => \$FILEHANDLE },
        outfile_path     => { strict_type => 1, store => \$outfile_path },
        stderrfile_path  => { strict_type => 1, store => \$stderrfile_path },
        filter_file_path => { strict_type => 1, store => \$filter_file_path },
        invert_match => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$invert_match
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    ## grep
    my @commands = qw(grep);    #Stores commands depending on input parameters

    ## Options
    if ($invert_match) {

        push @commands, '--invert-match';
    }
    if ($filter_file_path) {

        push @commands, '--file=' . $filter_file_path;
    }

    ## Infile

    push @commands, $infile_path;

    ## Outfile
    if ($outfile_path) {

        push @commands, '> ' . $outfile_path;
    }
    if ($stderrfile_path) {

        # Redirect stderr output to program specific stderr file
        push @commands, '2> ' . $stderrfile_path;
    }
    if ($FILEHANDLE) {

        my $SPACE = q{ };
        print {$FILEHANDLE} join( $SPACE, @commands ) . $SPACE;
    }
    return @commands;
}


1;
