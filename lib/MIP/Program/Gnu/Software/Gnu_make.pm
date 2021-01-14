package MIP::Program::Gnu::Software::Gnu_make;

use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

## CPANM
use Readonly;

## MIPs lib/
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ gnu_make };
}

## Constants
Readonly my $SPACE => q{ };

sub gnu_make {

## Function : Perl wrapper for writing make commands. Based on gnu make version 3.81
## Returns  : @commands
## Arguments: makefile_dir            => Path to dir with make file
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $filehandle             => Filehandle to write to
##          : $test                   => Run test

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $makefile_dir;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $filehandle;
    my $test;

    ## Default(s)

    my $tmpl = {
        makefile_dir => {
            strict_type => 1,
            store       => \$makefile_dir,
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path,
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        filehandle => {
            store => \$filehandle,
        },
        test => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$test,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{make};

    if ($makefile_dir) {
        push @commands, q{--directory=} . $makefile_dir;
    }

    if ($test) {
        push @commands, q{test};
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
            filehandle   => $filehandle,
        }
    );

    return @commands;
}

1;
