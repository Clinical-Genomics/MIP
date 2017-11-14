package MIP::Program::Variantcalling::Peddy;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
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
    our @EXPORT_OK = qw{ peddy };
}

## Constants
Readonly my $SPACE => q{ };

sub peddy {

## Function : Perl wrapper for writing peddy recipe to already open $FILEHANDLE or return commands array. Based on peddy 0.2.9.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $infile_path            => Infile path
##          : $outfile_prefix_path    => Outfile path
##          : $family_file_path       => Family file path
##          : $processor_number       => Number of processors to use
##          : $plot                   => Generate plots

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $infile_path;
    my $outfile_prefix_path;
    my $family_file_path;

    ## Default(s)
    my $processor_number;
    my $plot;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
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
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_prefix_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_prefix_path
        },
        family_file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_file_path
        },
        processor_number => {
            default     => 4,
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$processor_number
        },
        plot => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$plot
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{python -m peddy};

    ## Options
    if ($processor_number) {

        push @commands, q{--procs} . $SPACE . $processor_number;
    }

    if ($plot) {

        push @commands, q{--plot};
    }

    ## Outfile
    if ($outfile_prefix_path) {

        push @commands, q{--prefix} . $SPACE . $outfile_prefix_path;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    if ($family_file_path) {

        push @commands, $family_file_path;
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
