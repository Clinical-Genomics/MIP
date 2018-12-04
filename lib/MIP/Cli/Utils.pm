package MIP::Cli::Utils;

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
    our @EXPORT_OK = qw{ run };
}

## Constants
Readonly my $SPACE => q{ };

sub run {

## Function : Provide broadcast that there are more subcommands to go. MooseX::App required sub.
## Returns  :

    my ($arg_href) = @_;

    say {*STDERR} q{Please choose an subcommand to start the analysis};
    return;
}

1;
