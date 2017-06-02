package MIP::Gnu::Software::Gnu_grep;

use strict;
use warnings;
use warnings qw(FATAL utf8);
use utf8;    #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use autodie;
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

use FindBin qw($Bin);  #Find directory of script
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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(gnu_grep);
}

sub gnu_grep {

##gnu_grep

##Function : Perl wrapper for writing grep recipe to already open $FILEHANDLE or return commands array. Based on grep 2.6.3
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $filter_file_path, $stderrfile_path_append, $invert_match
##         : $FILEHANDLE             => Filehandle to write to
##         : $infile_path            => Infile path
##         : $outfile_path           => Outfile path
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append stderr info to file
##         : $filter_file_path       => Obtain patterns from file, one per line
##         : $invert_match           => Invert the sense of matching, to select non-matching lines

    my ($arg_href) = @_;

    ## Default(s)
    my $stderrfile_path_append;
    my $invert_match;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $filter_file_path;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        FILEHANDLE       => { store       => \$FILEHANDLE },
        outfile_path     => { strict_type => 1, store => \$outfile_path },
        stderrfile_path  => { strict_type => 1, store => \$stderrfile_path },
        filter_file_path => { strict_type => 1, store => \$filter_file_path },
		stderrfile_path_append =>
		{ strict_type => 1, store => \$stderrfile_path_append },
        invert_match => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$invert_match
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $SPACE = q{ };

    ### grep
    ## Stores commands depending on input parameters
    my @commands = qw(grep);

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

    push @commands, unix_standard_streams({stderrfile_path => $stderrfile_path,
					   stderrfile_path_append => $stderrfile_path_append,
					  });

    unix_write_to_file({commands_ref => \@commands,
			separator => $SPACE,
			FILEHANDLE => $FILEHANDLE,
		       });
    return @commands;
}

1;
