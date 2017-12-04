package MIP::Program::Trimming::Cutadapt;

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
    our @EXPORT_OK = qw{ cutadapt };
}

## Constants
Readonly my $SPACE => q{ };

sub cutadapt {

## Function : Perl wrapper for cutadapt.
## Returns  : @commands
## Arguments: $adapter_3_prime        => Sequence of an adapter ligated to the 3' end (paired data: of the first read).
##          : $adapter_5_prime        => Sequence of an adapter ligated to the 5' end (paired data: of the first read).
##          : $FILEHANDLE             => Filehandle to write to
##          : $outputfile_path        => Write trimmed reads to FILE
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $adapter_3_prime;
    my $adapter_5_prime;
    my $FILEHANDLE;
    my $outputfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)

    my $tmpl = {
        adapter_3_prime => {
            default     => undef,
            allow       => [ undef, qr/ ^\w+$ /xsm ],
            strict_type => 1,
            store       => \$adapter_3_prime
        },
        adapter_5_prime => {
            default     => undef,
            allow       => [ undef, qr/ ^\w+$ /xsm ],
            strict_type => 1,
            store       => \$adapter_5_prime
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        outputfile_path => {
            strict_type => 1,
            store       => \$outputfile_path,
        },
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{cutadapt};

    if ($adapter_3_prime) {

        push @commands, q{--adapter=} . $adapter_3_prime;
    }
    if ($adapter_5_prime) {

        push @commands, q{--front=} . $adapter_5_prime;
    }

    if ($outputfile_path) {

        push @commands, q{--output=} . $outputfile_path;
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
