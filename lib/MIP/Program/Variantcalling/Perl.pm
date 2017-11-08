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

    my @commands;

    ## Compose $regexp
    # Don't apply to comments
    my $regexp = q{if($_=~/^#/) {print $_;}} . $SPACE;

    # Substitute IUPAC code with N to not break vcf specifications (GRCh38)
    $regexp .= q?else { @F[4] =~ s/W|K|Y|R|S|M/N/g;? . $SPACE;

    if ($xargs) {

        # Escape chars are needed in front of separators
        $regexp = q?\'? . $regexp . q?print join(\"\\\t\", @F), \"\\\n\"; }\' ?;
    }
    else {

        $regexp = q?'? . $regexp . q?print join("\t", @F), "\n"; }' ?;
    }

    # Execute perl over expression
    $regexp = q{perl -nae} . $SPACE . $regexp;

    # Add Pipe in front of regexp
    $regexp = $PIPE . $SPACE . $regexp;

    push @commands, $regexp;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return;
}

1;
