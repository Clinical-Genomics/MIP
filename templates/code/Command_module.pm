package MIP::PATH::TO::MODULE;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;    #Allow unicode characters in this script
use open qw{  :encoding(UTF-8) :std};
use charnames qw{ :full :short };
use Carp;
use English qw{-no_match_vars};
use Params::Check qw{check allow last_error};

use Readonly;

use FindBin qw{$Bin};    #Find directory of script
use File::Basename qw{dirname};
use File::Spec::Functions qw{catdir};

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Unix::Standard_streams qw{unix_standard_streams};
use MIP::Unix::Write_to_file qw{unix_write_to_file};

BEGIN {
    require Exporter;
    use base qw{Exporter};

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ space separated subroutines };
}

## Constants
Readonly my $SPACE => q{ };

sub name_of_subroutine {

## name_of_subroutine

## Function : Perl wrapper for generic commands module.
## Returns  : "@commands"
## Arguments: $stdoutfile_path, $stderrfile_path, stderrfile_path_append, $FILEHANDLE
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Filehandle to write to

    my ($arg_href) = @_;

    ## Default(s)

    ## Flatten argument(s)
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;

    my $tmpl = {
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE => { store => \$FILEHANDLE },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ BASE_COMMAND };

    push @commands,
      unix_standard_streams(
        {
            stdoutfile_path        => $stdoutfile_path,
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
