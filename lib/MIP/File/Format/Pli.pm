package MIP::File::Format::Pli;

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

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $TAB };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ load_pli_file };
}

sub load_pli_file {

## Function : Load plI file values
## Returns  :
## Arguments: $infile_path    => Infile path
##          : $pli_score_href => Pli scores hash

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $pli_score_href;

    my $tmpl = {
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
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

    use MIP::Io::Read qw{ read_from_file };

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    $log->info(qq{Loading pli value file: $infile_path});

    my @lines = read_from_file(
        {
            chomp  => 1,
            format => q{line_by_line},
            path   => $infile_path,
        }
    );

  LINE:
    foreach my $line (@lines) {

        my ( $hgnc_symbol, $pli_score ) = split $TAB, $line;

        ## Skip header
        next LINE if ( $pli_score eq q{pLI} );

        ## Set rounded pli score to hash
        $pli_score_href->{$hgnc_symbol} = sprintf q{%.2f}, $pli_score;
    }

    $log->info(q{Loading pli value file: Done});

    return 1;
}

1;
