package MIP::Cli::Utils;

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

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ run };
}

sub run {

## Function : Provide broadcast that there are more subcommands to go. MooseX::App required sub.
## Returns  :

    my ($arg_href) = @_;

    say {*STDERR} q{Please choose an subcommand to start the analysis};
    return;
}

1;
