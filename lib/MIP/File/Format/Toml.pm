package MIP::File::Format::Toml;

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
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ load_toml };
}

sub load_toml {

## Function : Loads a TOML file into an arbitrary hash and returns it.
## Returns  : %toml
## Arguments: $path => Toml file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $path;

    my $tmpl = {
        path => {
            defined     => 1,
            required    => 1,
            store       => \$path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use TOML::Parser;

    # Create new parser
    my $parser_toml = TOML::Parser->new;

    # Load file info into hash
    my %toml = %{ $parser_toml->parse_file($path) };

    return %toml;
}

1;
