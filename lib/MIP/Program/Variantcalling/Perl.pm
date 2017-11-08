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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ replace_iupac };
}

## Constants
Readonly my $SPACE => q{ };
Readonly my $PIPE  => q{|};

sub replace_iupac {

## replace_iupac

## Function : Replace the IUPAC code in alternative allels with N for input stream and writes to stream.
## Returns  :
## Arguments: $FILEHANDLE      => Sbatch filehandle to write to
##          : $stderrfile_path => Stderr path to errors write to
##          : $xargs           => Write on xargs format

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $stderrfile_path;
    my $FILEHANDLE;

    ## Default(s)
    my $xargs;

    my $tmpl = {
        FILEHANDLE => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$FILEHANDLE
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        xargs           => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$xargs
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Pipe
    print {$FILEHANDLE} $PIPE . $SPACE;

    ## Execute perl
    print {$FILEHANDLE} q{perl -nae} . $SPACE;

    #Add escape char
    if ($xargs) {

        #substitute IUPAC code with N to not break vcf specifications (GRCh38)
        print {$FILEHANDLE}
q?\'if($_=~/^#/) {print $_;} else { @F[4] =~ s/W|K|Y|R|S|M/N/g; print join(\"\\\t\", @F), \"\\\n\"; }\' ?;
    }
    else {

        #substitute IUPAC code with N to not break vcf specifications (GRCh38)
        print {$FILEHANDLE}
q?'if($_=~/^#/) {print $_;} else { @F[4] =~ s/W|K|Y|R|S|M/N/g; print join("\t", @F), "\n"; }' ?;
    }

    #Redirect output to program specific stderr file
    if ($stderrfile_path) {

        print {$FILEHANDLE} q{2>>} . $SPACE . $stderrfile_path . $SPACE;
    }

    return;
}

1;
