package MIP::Program::Gnu::Software::Gnu_sed;

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
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ gnu_sed };

}

## Constants
Readonly my $SPACE => q{ };

sub gnu_sed {

##Function : Perl wrapper for writing sed recipe to already open $filehandle or return commands array. Based on sed 4.2.1.
##Returns  : "@commands"
##Arguments: $filehandle       => Filehandle to write to
##         : $infile_path      => Infile path
##         : $inplace_edit     => Edit file in place
##         : $outfile_path     => Outfile path
##         : $stderrfile_path_append => Append stderr info to file
##         : $stderrfile_path  => Stderrfile path
##         : $script           => Script to edit infile stream

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $inplace_edit;
    my $outfile_path;
    my $stderrfile_path;
    my $script;
    my $stderrfile_path_append;

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
        outfile_path => {
            strict_type => 1,
            store       => \$outfile_path
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        script => {
            strict_type => 1,
            store       => \$script
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ##Stores commands depending on input parameters
    my @commands = q{sed};

    ## Options
    if ($inplace_edit) {

        push @commands, q{-i};
    }

    if ($script) {

        push @commands, $script;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    ## Outfile
    if ($outfile_path) {

        push @commands, q{>} . $SPACE . $outfile_path;
    }

    push @commands,
      unix_standard_streams(
        {
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
