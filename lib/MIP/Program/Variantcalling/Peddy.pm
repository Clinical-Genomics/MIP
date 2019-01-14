package MIP::Program::Variantcalling::Peddy;

use 5.026;
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
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ peddy };
}

## Constants
Readonly my $SPACE => q{ };

sub peddy {

## Function : Perl wrapper for writing peddy recipe to already open $FILEHANDLE or return commands array. Based on peddy 0.2.9.
## Returns  : @commands
## Arguments: $case_file_path         => Family file path
##          : $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $outfile_prefix_path    => Outfile path
##          : $plot                   => Generate plots
##          : $processor_number       => Number of processors to use
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_file_path;
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_prefix_path;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $plot;
    my $processor_number;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_prefix_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_prefix_path,
            strict_type => 1,
        },
        case_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$case_file_path,
            strict_type => 1,
        },
        processor_number => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 4,
            store       => \$processor_number,
            strict_type => 1,
        },
        plot => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$plot,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ python -m peddy };

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

    if ($case_file_path) {

        push @commands, $case_file_path;
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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

1;
