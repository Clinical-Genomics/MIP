package MIP::Unix::Standard_streams;

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
use MIP::Unix::Write_to_file qw(unix_write_to_file);

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ unix_standard_streams };
}

## Constants
Readonly my $SPACE => q{ };

$Params::Check::PRESERVE_CASE =
  1;    #Do not convert to lower case - required to pass $FILEHANDLE

sub unix_standard_streams {

## Function : Perl wrapper for writing unix standard_streams recipe to already open $FILEHANDLE or return commands array.
## Returns  : @commands
## Arguments: $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        FILEHANDLE      => { store       => \$FILEHANDLE, },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path, },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append, },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ### Standard streams helper module
    ##Stores commands depending on input parameters
    my @commands;

    ## Options
    if ($stdoutfile_path) {

        # Redirect stdout to program specific stdout file
        push @commands, q{1>} . $SPACE . $stdoutfile_path;
    }
    if ($stderrfile_path) {

        # Redirect stderr output to program specific stderr file
        push @commands, q{2>} . $SPACE . $stderrfile_path;
    }
    if ($stderrfile_path_append) {

        # Redirect and append stderr output to program specific stderr file
        push @commands, q{2>>} . $SPACE . $stderrfile_path_append;
    }
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
