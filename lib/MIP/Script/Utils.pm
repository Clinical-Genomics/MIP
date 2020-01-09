package MIP::Script::Utils;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ basename };
use File::Spec::Functions qw{ catdir };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## CPANM
use Readonly;
use autodie;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.09;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ help print_parameter_defaults };
}

sub help {

## Function : Print help text and exit with supplied exit code
## Returns  :
## Arguments: $exit_code => Exit code
##          : $usage     => Help text

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $usage;

    ## Default(s)
    my $exit_code;

    my $tmpl = {
        exit_code => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$exit_code,
            strict_type => 1,
        },
        usage => {
            defined     => 1,
            required    => 1,
            store       => \$usage,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    say {*STDOUT} $usage;
    exit $exit_code;
}

sub print_parameter_defaults {

## Function : Print all parameters and their default values
## Returns  :
## Arguments: $colored                 => Colorize output
##          : $index                   => Display array indices
##          : $output                  => Output stream
##          : $parameter_href          => Holds all parameters {REF}
##          : $print_parameter_default => Print parameter default
##          : $scalar_quotes           => Quote symbols to enclose scalar values

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;

    ## Default
    my $colored;
    my $index;
    my $output;
    my $print_parameter_default;
    my $scalar_quotes;

    my $tmpl = {
        colored => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$colored,
            strict_type => 1,
        },
        index => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$index,
            strict_type => 1,
        },
        output => {
            default     => q{stderr},
            store       => \$output,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        print_parameter_default => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$print_parameter_default,
            strict_type => 1,
        },
        scalar_quotes => {
            allow       => [ undef, q{"} ],
            default     => q{"},
            store       => \$scalar_quotes,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Data::Printer;

    if ($print_parameter_default) {

        ## Print parameter hash an exit
        say {*STDERR} q{Default values from config file:};
        p(
            %{$parameter_href},
            colored       => $colored,
            index         => $index,
            scalar_quotes => $scalar_quotes,
            output        => $output,
        );

        exit 0;
    }
    return;
}

1;
