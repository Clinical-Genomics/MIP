package MIP::Gnu::Bash;

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
    our @EXPORT_OK = qw(gnu_cd);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

sub gnu_cd {

##gnu_cd

##Function : Perl wrapper for writing cd recipe to already open $FILEHANDLE or return commands array. Based on cd 4.0
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $directory_path, $stderrfile_path, append_stderr_info
##         : $FILEHANDLE      => Filehandle to write to
##         : $directory_path  => Directory path
##         : $stderrfile_path => Stderrfile path
##         : $append_stderr_info => Append stderr info to file

    my ($arg_href) = @_;

    ## Default(s)
    my $append_stderr_info;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $directory_path;
    my $stderrfile_path;

    my $tmpl = {
        FILEHANDLE      => { store       => \$FILEHANDLE },
        directory_path  => { strict_type => 1, store => \$directory_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	append_stderr_info => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$append_stderr_info
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    my $SPACE = q{ };

    ### cd
    ##Stores commands depending on input parameters
    my @commands = qw(cd);

    ## Options
    if ($directory_path) {

        push @commands, $directory_path;
    }
    if ($stderrfile_path) {

        if ($append_stderr_info) {

            # Redirect and append stderr output to program specific stderr file
            push @commands, '2>> ' . $stderrfile_path;
        }
        else {

            # Redirect stderr output to program specific stderr file
            push @commands, '2> ' . $stderrfile_path;
        }
    }
    if ($FILEHANDLE) {

        print {$FILEHANDLE} join( $SPACE, @commands ) . $SPACE;
    }
    return @commands;
}

1;
