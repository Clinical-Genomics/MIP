package MIP::File::Format::Yaml;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

## CPANM
use Readonly;
use YAML;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ load_yaml write_yaml };
}

## Constants
Readonly my $DOUBLE_QUOTE => q{"};
Readonly my $SINGLE_QUOTE => q{'};
Readonly my $SPACE        => q{ };
Readonly my $NEWLINE      => qq{\n};

# Localize the current value
local $YAML::QuoteNumericStrings = $YAML::QuoteNumericStrings;

# Force numeric values to strings in YAML representation
$YAML::QuoteNumericStrings = 1;

sub load_yaml {

## load_yaml

## Function : Loads a YAML file into an arbitrary hash and returns it.
## Returns  : %yaml
## Arguments: $yaml_file
##          : $yaml_file => The yaml file to load

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $yaml_file;

    my $tmpl = {
        yaml_file => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$yaml_file
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %yaml;

    open my $YAML, q{<}, $yaml_file
      or croak q{cannot open} . $SPACE . $yaml_file . q{:} . $OS_ERROR,
      $NEWLINE;

    ## Load hashreference as hash
    %yaml = %{ YAML::LoadFile($yaml_file) };

    close $YAML;

    return %yaml;
}

sub write_yaml {

## write_yaml

## Function : Writes a YAML hash to file
## Returns  : ""
## Arguments: $yaml_href, $yaml_file_path
##          : $yaml_href      => The hash to dump {REF}
##          : $yaml_file_path => The yaml file to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $yaml_href;
    my $yaml_file_path;

    my $tmpl = {
        yaml_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$yaml_href
        },
        yaml_file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$yaml_file_path
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    open my $YAML, q{>}, $yaml_file_path
      or croak q{Cannot open}
      . $SPACE
      . $SINGLE_QUOTE
      . $DOUBLE_QUOTE
      . $yaml_file_path
      . $DOUBLE_QUOTE
      . $SINGLE_QUOTE . q{:}
      . $OS_ERROR
      . $NEWLINE;

    say {$YAML} Dump($yaml_href);
    close $YAML;
    return;
}

1;
