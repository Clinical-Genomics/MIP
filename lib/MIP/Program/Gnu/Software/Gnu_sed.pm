package MIP::Program::Gnu::Software::Gnu_sed;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

## MIPs lib/
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ gnu_sed };

}

## Constants
Readonly my $SPACE => q{ };

sub gnu_sed {

##Function : Perl wrapper for writing sed recipe to already open $filehandle or return commands array. Based on sed 4.2.1.
##Returns  : "@commands"
##Arguments: $filehandle             => Filehandle to write to
##         : $infile_path            => Infile path
##         : $inplace_edit           => Edit file in place
##         : $script                 => Script to edit infile stream
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append stderr info to file
##         : $stdoutfile_path        => Stdoutfile path
##         : $stdoutfile_path_append => Append stdout info to file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $inplace_edit;
    my $script;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $stdoutfile_path_append;

    my $tmpl = {
        filehandle => {
            store => \$filehandle
        },
        infile_path => {
            strict_type => 1,
            store       => \$infile_path
        },
        inplace_edit => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$inplace_edit,
            strict_type => 1,
        },
        script => {
            strict_type => 1,
            store       => \$script
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path
        },
        stdoutfile_path_append => {
            strict_type => 1,
            store       => \$stdoutfile_path_append
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ##Stores commands depending on input parameters
    my @commands = q{sed};

    if ($inplace_edit) {

        push @commands, q{-i};
    }
    if ($script) {

        push @commands, $script;
    }
    if ($infile_path) {

        push @commands, $infile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
            stdoutfile_path_append => $stdoutfile_path_append,
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
