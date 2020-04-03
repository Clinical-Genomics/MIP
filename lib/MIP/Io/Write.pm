package MIP::Io::Write;

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
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ write_to_file };
}

sub write_to_file {

## Function : Write to file from perl data structure
## Returns  :
## Arguments: $data_href => Data to write to file
##          : $format    => File format
##          : $path      => File path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $data_href;
    my $format;
    my $path;

    my $tmpl = {
        data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$data_href,
            strict_type => 1,
        },
        format => {
            allow       => [qw{ toml yaml }],
            default     => 1,
            store       => \$format,
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

    use MIP::Yaml qw{ write_yaml };
    use MIP::Toml qw{ write_toml };

    my %file_api = (
        toml => {
            arg_href => {
                path      => $path,
                data_href => $data_href,
            },
            method => \&write_toml,

        },
        yaml => {
            arg_href => {
                path      => $path,
                data_href => $data_href,
            },
            method => \&write_yaml,
        },
    );

    $file_api{$format}{method}->( { %{ $file_api{$format}{arg_href} } } );

    return;
}

1;
