package MIP::Cli::Mip;

use Carp;
use open qw{ :encoding(UTF-8) :std };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use MooseX::App qw{ BashCompletion Color Typo };

## MIPs lib/
use MIP::Cli::Mip;

our $VERSION = 0.01;

## Enable strict mode
app_strict 1;

## Allows one to specify multiple values with one key
app_permute 1;

## Define, check and get Cli supplied parameters
_build_usage();

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

    option(
        q{config_file} => (
            cmd_aliases   => [qw{ config }],
            documentation => q{YAML config file for analysis parameters},
            is            => q{rw},
            isa           => q{Str},
            required      => 1,
        )
    );

   option(
        q{log_file} => (
            cmd_aliases   => [qw{ log }],
            documentation => q{Log file},
            is            => q{rw},
            isa           => q{Str},
        )
    );

    option(
        q{version} => (
            cmd_aliases   => q{v},
            documentation => q{Show version},
            is            => q{rw},
            isa           => q{Bool},
        )
    );
    return;
}

1;
