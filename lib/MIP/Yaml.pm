package MIP::Yaml;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie;
use Readonly;
use YAML qw{ Dump LoadFile };

## MIPs lib/
use MIP::Constants qw{ $COLON $DOUBLE_QUOTE $NEWLINE $SINGLE_QUOTE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ load_yaml write_yaml };
}

sub load_yaml {

## Function : Loads a YAML file into an arbitrary hash and returns it.
## Returns  : %yaml
## Arguments: $path => Yaml file path to load

    my ($arg_href) = @_;

    ##Flatten argument(s)
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

    my %yaml;

    open my $YAML, q{<}, $path
      or croak q{cannot open} . $SPACE . $path . $COLON . $OS_ERROR,
      $NEWLINE;

    ## Load hashreference as hash
    %yaml = %{ LoadFile($path) };

    close $YAML;

    return %yaml;
}

sub write_yaml {

## Function : Writes data hash to file
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

    # Localize the current value
    local $YAML::QuoteNumericStrings = 1;

    open my $YAML, q{>}, $path
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

    say {$YAML} Dump($data_href);
    close $YAML;
    return;
}

1;
