package MIP::Program::Variantcalling::Perl;

use strict;
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

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ replace_iupac };
}

## Constants
Readonly my $SPACE => q{ };
Readonly my $PIPE  => q{|};

sub replace_iupac {

## Function : Replace the IUPAC code in alternative allels with N for input stream and writes to stream.
## Returns  :
## Arguments: $filehandle      => Sbatch filehandle to write to
##          : $stderrfile_path => Stderr path to errors write to
##          : $stdoutfile_path => Stdoutfile path
##          : $xargs           => Write on xargs format

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $stderrfile_path;
    my $stdoutfile_path;

    ## Default(s)
    my $xargs;

    my $tmpl = {
        filehandle => {
            defined  => 1,
            required => 1,
            store    => \$filehandle,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        xargs => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$xargs,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Compose $regexp
    # Execute perl
    my $regexp = q{perl -nae} . $SPACE;

    ## Substitute IUPAC code with N to not break vcf specifications (grch38)
    if ($xargs) {

  # Print comment lines as they are but add escape char at the beginning of the expression
        $regexp .= q{\'if($_=~/^#/) {print $_;}} . $SPACE;

        # Escape chars are needed in front of separators
        $regexp .=
          q?else { @F[4] =~ s/W|K|Y|R|S|M/N/g; print join(\"\\\t\", @F), \"\\\n\"; }\'?
          . $SPACE;
    }
    else {

        # Print comment lines as they are
        $regexp .= q{'if($_=~/^#/) {print $_;}} . $SPACE;

        # Escape chars are NOT needed in front of separators
        $regexp .=
          q?else { @F[4] =~ s/W|K|Y|R|S|M/N/g; print join("\t", @F), "\n"; }'? . $SPACE;
    }

    print {$filehandle} $regexp;

    unix_standard_streams(
        {
            filehandle      => $filehandle,
            stderrfile_path => $stderrfile_path,
            stdoutfile_path => $stdoutfile_path,

        }
    );
    return;
}

1;
