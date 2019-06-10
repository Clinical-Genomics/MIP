package MIP::File::Format::Dragen;

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
use MIP::Constants qw{ $COMMA $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ create_dragen_fastq_list_sample_id };
}

sub create_dragen_fastq_list_sample_id {

## Function : Create dragen fastq list sample id file for supplying input fastq files to dragen
## Returns  :
## Arguments: $FILEHANDLE            => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;

    ## Default(s)

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @fastq_list_headers = join $COMMA, qw{ RGID RGSM RGLB Lane Read1File Read2File };

    return;
}

1;
