package MIP::File::Format::Pli;

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
use MIP::Constants qw{ $COLON $NEWLINE $SPACE $TAB };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ load_pli_file };
}

sub load_pli_file {

## Function : Load plI file values
## Returns  :
## Arguments: $infile_path    => Infile path
##          : $log            => Log object
##          : $pli_score_href => Pli scores hash

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $log;
    my $pli_score_href;

    my $tmpl = {
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        pli_score_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pli_score_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $filehandle = IO::Handle->new();

    open $filehandle, q{<}, $infile_path
      or $log->logdie( q{Cannot open } . $infile_path . $COLON . $OS_ERROR, $NEWLINE );

  LINE:
    while (<$filehandle>) {

        chomp;

        ## Unpack line
        my $line = $_;

        ## Get hgnc symbol and pli score
        my ( $hgnc_symbol, $pli_score ) = split $TAB, $line;

        ## Skip header
        next LINE if ( $pli_score eq q{pLI} );

        ## Set rounded pli score to hash
        $pli_score_href->{$hgnc_symbol} = sprintf q{%.2f}, $pli_score;
    }
    close $filehandle;
    return 1;
}

1;
