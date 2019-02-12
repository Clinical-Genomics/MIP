package MIP::Constants;

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

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ $AMPERSAND %ANALYSIS $ASTERISK $COLON $COMMA $DASH $DOT $EMPTY_STR $MIP_VERSION $NEWLINE $SPACE $SEMICOLON $UNDERSCORE };
}

## Constants
## Set MIP version
## Constants
Readonly our $MIP_VERSION => q{v7.0.1};

## Symbols
Readonly our $AMPERSAND  => q{&};
Readonly our $ASTERISK   => q{*};
Readonly our $COMMA      => q{,};
Readonly my $COLON       => q{:};
Readonly our $DASH       => q{-};
Readonly our $DOT        => q{.};
Readonly our $EMPTY_STR  => q{};
Readonly our $NEWLINE    => qq{\n};
Readonly our $SPACE      => q{ };
Readonly our $SEMICOLON  => q{;};
Readonly our $UNDERSCORE => q{_};

## Analysis
Readonly our %ANALYSIS => (
    ANNOTATION_DISTANCE    => 5000,
    ANNOTATION_DISTANCE_MT => 0,
);

1;
