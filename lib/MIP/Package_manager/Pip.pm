package MIP::Package_manager::Pip;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use Readonly;

## MIPs lib/
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ pip_install };
}

## Constants
Readonly my $SPACE => q{ };

sub pip_install {

## pip_install

## Function : Perl wrapper for writing pip install command. Based on pip v9.0.1.
## Returns  : "@commands"

## Arguments: $packages_ref, $quiet, $requirement, $editable, $stdoutfile_path, $stderrfile_path, $stderrfile_path_append, $FILEHANDLE
##          : $packages_ref           => Array of packages to be installed {REF}
##          : $quiet                  => Quiet output
##          : $requirement            => Install from a requirement file
##          : $editable               => Install in editable mode from a local project path or a VCS url.
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $packages_ref;
    my $quiet;
    my $requirement;
    my $editable;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;

    ## Default(s)

    my $tmpl = {
        packages_ref => {
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$packages_ref
        },
        quiet => {
            store       => \$quiet
        },
        requirement => {
            strict_type => 1,
            store       => \$requirement
        },
        editable => {
            strict_type => 1,
            store       => \$editable
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{pip install};

    if ($quiet) {
        push @commands, q{--quiet};
    }

    if ($requirement) {
        push @commands, q{--requirement} . $SPACE . $requirement;
    }
    elsif ($editable) {
        push @commands, q{--editable} . $SPACE . $editable;
    }
    else {
        push @commands, join $SPACE, @{$packages_ref};
    }

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
