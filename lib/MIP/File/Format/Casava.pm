package MIP::File::Format::Casava;

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

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ casava_header_regexp };
}

## Constants
Readonly my $SPACE => q{ };

sub casava_header_regexp {

## Function : Define casava fastq file header formats matching regexp
## Returns  : %casava_header_regexp
## Arguments:

    my ($arg_href) = @_;

    my %casava_header_regexp = (
        q{1.4} =>
q&perl -nae 'chomp; my ($instrument_id, $run_number, $flowcell, $lane, $tile, $x_pos, $y_pos, $direction) = /^(@[^:]*):(\d+):(\w+):(\d+):(\d+):(\d+):(\d+)[\/](\d+)/; if($instrument_id) { print join " ", ($instrument_id, $run_number, $flowcell, $lane, $tile, $x_pos, $y_pos, $direction);} if($.=1) {last;}'&,
        q{1.4_header_features} =>
          [qw{ instrument_id run_number flowcell lane tile x_pos y_pos direction }],
        q{1.8} =>
q&perl -nae 'chomp; my ($instrument_id, $run_number, $flowcell, $lane, $tile, $x_pos, $y_pos, $direction, $filtered, $control_bit, $index,) = /^(@[^:]*):(\d+):(\w+):(\d+):(\d+):(\d+):(\d+)\s(\d+):(\w+):(\d+):(\w+)/; if($instrument_id) { print join " ", ($instrument_id, $run_number, $flowcell, $lane, $tile, $x_pos, $y_pos, $direction, $filtered, $control_bit, $index); } if($.=1) {last;}'&,
        q{1.8_header_features} => [
            qw{ instrument_id run_number flowcell lane tile x_pos y_pos direction filtered control_bit index }
        ],
        q{1.4_interleaved} =>
q&perl -nae 'chomp; if($.==1 or $.==5) { my ($instrument_id, $run_number, $flowcell, $lane, $tile, $x_pos, $y_pos, $direction) = /^(@[^:]*):(\d+):(\w+):(\d+):(\d+):(\d+):(\d+)[\/](\d+)/; if($instrument_id) { print $direction;} } elsif ($.==6) {last;}'&,
        q{1.8_interleaved} =>
q&perl -nae 'chomp; if($.==1 or $.==5) { my ($instrument_id, $run_number, $flowcell, $lane, $tile, $x_pos, $y_pos, $direction, $filtered, $control_bit, $index,) = /^(@[^:]*):(\d+):(\w+):(\d+):(\d+):(\d+):(\d+)\s(\d+):(\w+):(\d+):(\w+)/; if($instrument_id) { print $direction;} } elsif ($.==6) {last;}'&,
    );

    return %casava_header_regexp;
}

1;
