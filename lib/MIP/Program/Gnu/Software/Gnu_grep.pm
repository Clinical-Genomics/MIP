package MIP::Program::Gnu::Software::Gnu_grep;

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
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ gnu_grep };
}

sub gnu_grep {

##Function : Perl wrapper for writing grep recipe to already open $filehandle or return commands array. Based on grep 2.6.3
##Returns  : @commands
##Arguments: $count                  => Count
##         : $filehandle             => Filehandle to write to
##         : $filter_file_path       => Obtain patterns from file, one per line
##         : $infile_path            => Infile path
##         : $invert_match           => Invert the sense of matching, to select non-matching lines
##         : $pattern                => Pattern to match
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append stderr info to file
##         : $stdoutfile_path        => Stdoutfile path
##         : $word_regexp            => Select only those lines containing matches that form whole words

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $count;
    my $filehandle;
    my $filter_file_path;
    my $infile_path;
    my $pattern;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $invert_match;
    my $word_regexp;

    my $tmpl = {
        count => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$count,
            strict_type => 1,
        },
        filehandle       => { store => \$filehandle, },
        filter_file_path => { store => \$filter_file_path, strict_type => 1, },
        infile_path => {
            store       => \$infile_path,
            strict_type => 1,
        },
        invert_match => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$invert_match,
            strict_type => 1,
        },
        pattern         => { store => \$pattern,         strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        word_regexp => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$word_regexp,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ### grep
    ## Stores commands depending on input parameters
    my @commands = qw{ grep };

    ## Options
    if ($count) {

        push @commands, q{--count};
    }
    if ($invert_match) {

        push @commands, q{--invert-match};
    }
    if ($word_regexp) {

        push @commands, q{--word-regexp};
    }
    if ($filter_file_path) {

        push @commands, q{--file=} . $filter_file_path;
    }

    if ($pattern) {

        push @commands, $pattern;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

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
