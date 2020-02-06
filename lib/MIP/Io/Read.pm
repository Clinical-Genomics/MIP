package MIP::Io::Read;

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

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ read_from_file };
}

sub read_from_file {

## Function : Read from file and return perl data structure
## Returns  : %hash
## Arguments: $format => File format
##          : $path   => File path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $format;
    my $path;

    my $tmpl = {
        path => {
            defined     => 1,
            required    => 1,
            store       => \$path,
            strict_type => 1,
        },
        format => {
            allow       => [qw{ toml yaml }],
            default     => 1,
            store       => \$format,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Toml qw{ load_toml };
    use MIP::File::Format::Yaml qw{ load_yaml };

    my %file_api = (
        yaml => {
            method   => \&load_yaml,
            arg_href => { path => $path, },
        },
        toml => {
            method   => \&load_toml,
            arg_href => { path => $path },
        },
    );

    my %hash = $file_api{$format}{method}->( { %{ $file_api{$format}{arg_href} } } );

    return %hash;
}

1;
