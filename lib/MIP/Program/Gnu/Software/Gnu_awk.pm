package MIP::Program::Gnu::Software::Gnu_awk;

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
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COLON $COMMA $EQUALS $SEMICOLON $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ gnu_awk };
}

sub gnu_awk {

##Function : Perl wrapper for writing awk recipe to already open $filehandle or return commands array. Based on awk 4.0.2
##Returns  : @commands
##Arguments: $field_separator        => Field separator
##         : $filehandle             => Filehandle to write to
##         : $infile_path            => Infile path
##         : $outfile_path           => Outfile path
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append stderr info to file
##         : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $field_separator;
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        field_separator => {
            allow => [ $COLON, $COMMA, $SPACE, $SEMICOLON ],
            store => \$field_separator,
        },
        filehandle  => { store => \$filehandle, },
        infile_path => {
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path => {
            store       => \$outfile_path,
            strict_type => 1,
        },
        stderrfile_path        => { store => \$stderrfile_path,        strict_type => 1, },
        stderrfile_path_append => { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path        => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ### awk
    ## Stores commands depending on input parameters
    my @commands = qw{ awk };

    ## Options
    if ($field_separator) {

        push @commands, q{'-F[} . $field_separator . q{]'};
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    ## Outfile
    if ($outfile_path) {

        push @commands, q{>} . $outfile_path;
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
