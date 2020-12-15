package MIP::Toml;

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
use MIP::Constants qw{ $COLON $DOUBLE_QUOTE $NEWLINE $SINGLE_QUOTE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ load_toml write_toml };
}

sub load_toml {

## Function : Loads a TOML file into an arbitrary hash and returns reference to hash
## Returns  : \%toml
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

    return \%toml;
}

sub write_toml {

## Function : Writes data hash to file in toml format
## Returns  :
## Arguments: $data_href => Data hash to dump {REF}
##          : $path      => Yaml file path to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $data_href;
    my $path;

    my $tmpl = {
        data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$data_href,
            strict_type => 1,
        },
        path => {
            defined     => 1,
            required    => 1,
            store       => \$path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use TOML::Tiny qw{ to_toml };

    my $toml_str = to_toml($data_href);

    open my $TOML, q{>}, $path
      or croak q{Cannot open}
      . $SPACE
      . $SINGLE_QUOTE
      . $DOUBLE_QUOTE
      . $path
      . $DOUBLE_QUOTE
      . $SINGLE_QUOTE
      . $COLON
      . $OS_ERROR
      . $NEWLINE;

    say {$TOML} $toml_str;
    close $TOML;
    return;
}

1;
