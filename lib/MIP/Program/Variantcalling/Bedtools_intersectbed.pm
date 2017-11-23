package MIP::Program::Variantcalling::Bedtools_intersectbed;

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
    our @EXPORT_OK = qw{ bedtools_intersectbed };
}

## Constants
Readonly my $SPACE => q{ };

sub bedtools_intersectbed {

## Function : Perl wrapper for writing Bedtools intersectbed recipe to already open $FILEHANDLE or return commands array. Based on Bedtools 2.26.0.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $intersectfile_path     => Intersect file (-b)
##          : $outfile_path           => Outfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $with_header            => Include header from infile in output

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $intersectfile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $with_header;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        infile_path        => { strict_type => 1, store => \$infile_path },
        intersectfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$intersectfile_path
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path,
        },
        with_header => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$with_header
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{intersectBed};

    ## Options
    # Include header
    if ($with_header) {

        push @commands, q{-header};
    }

    ## Infile
    if ($infile_path) {

        push @commands, q{-a} . $SPACE . $infile_path;
    }

    if ($intersectfile_path) {

        push @commands, q{-b} . $SPACE . $intersectfile_path;
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
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            commands_ref => \@commands,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
