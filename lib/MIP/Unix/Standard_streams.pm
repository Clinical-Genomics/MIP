package MIP::Unix::Standard_streams;

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
use MIP::Constants qw{ $SPACE };
use MIP::Unix::Write_to_file qw(unix_write_to_file);

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ unix_standard_streams };
}

sub unix_standard_streams {

## Function : Perl wrapper for writing unix standard_streams recipe to already open $filehandle or return commands array.
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path
##          : $stdoutfile_path_append => Append stdout info to file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;
    my $stdoutfile_path_append;

    my $tmpl = {
        filehandle      => { store => \$filehandle, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdinfile_path  => { store => \$stdinfile_path,  strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
        stdoutfile_path_append =>
          { store => \$stdoutfile_path_append, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ### Standard streams helper module
    ## Stores commands depending on input parameters
    my @commands;

    ## Options
    if ($stdinfile_path) {

        # Set stdin to program
        push @commands, q{<} . $SPACE . $stdinfile_path;
    }
    if ($stdoutfile_path) {

        # Redirect stdout to program specific stdout file
        push @commands, q{1>} . $SPACE . $stdoutfile_path;
    }
    if ($stdoutfile_path_append) {

        # Redirect and append stdout to program specific stdout file
        push @commands, q{1>>} . $SPACE . $stdoutfile_path_append;
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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

1;
