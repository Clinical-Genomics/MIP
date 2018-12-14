package MIP::Check::File;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use List::Util qw { any sum };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_interleaved };
}

## Constants
Readonly my $SPACE => q{ };
Readonly my $THREE => q{3};

sub check_interleaved {

## Function : Detect if fastq file is interleaved
## Returns  : "1(=interleaved)"
## Arguments: $file_path         => File to parse
##          : $log               => Log object
##          : $read_file_command => Command used to read file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $log;
    my $read_file_command;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        read_file_command => {
            defined     => 1,
            required    => 1,
            store       => \$read_file_command,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Casava qw{ casava_header_regexp };

    my %casava_header_regexp = casava_header_regexp();

    ## Select relevant regexps from hash
    my @regexps = (
        $casava_header_regexp{q{1.4_interleaved}},
        $casava_header_regexp{q{1.8_interleaved}},
    );

    ## Store return from regexp
    my $fastq_read_direction;

  REGEXP:
    foreach my $regexp (@regexps) {

        my $fastq_info_headers_cmd = qq{$read_file_command $file_path | $regexp;};

        $fastq_read_direction = `$fastq_info_headers_cmd`;
        last REGEXP if ($fastq_read_direction);
    }

    if ( not $fastq_read_direction ) {

        $log->fatal( q{Malformed fastq file: } . $file_path );
        $log->fatal(q{Could not find a read direction });
        exit 1;
    }

    my @fastq_read_directions = split //sxm, $fastq_read_direction;

    if ( any { /[^123]/sxm } @fastq_read_directions ) {

        $log->fatal(q{Malformed fastq file!});
        $log->fatal(
            q{Read direction is: } . join q{ and },
            @fastq_read_directions
              . q{, allowed entries are '1', '2', '3'. Please check fastq file}
              . $file_path
        );
        exit 1;
    }
    if ( sum(@fastq_read_directions) == $THREE ) {

        $log->info( q{Found interleaved fastq file: } . $file_path );
        return 1;
    }
    return;
}

1;
