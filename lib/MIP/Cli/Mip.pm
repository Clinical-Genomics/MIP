package MIP::Cli::Mip;

use 5.026;
use Carp;
use open qw{ :encoding(UTF-8) :std };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use MooseX::App qw{ BashCompletion Color MutexGroup Typo Version };
use MooseX::Types::Moose qw{ Str Int HashRef Num Bool ArrayRef };
use Moose::Util::TypeConstraints;

## MIPs lib/
use MIP::Cli::Mip;
use MIP::Constants qw{ $MIP_VERSION $MOOSEX_APP_SCEEN_WIDTH};

our $VERSION = $MIP_VERSION;

## Set screen width for usage message
$MooseX::App::Utils::SCREEN_WIDTH = $MOOSEX_APP_SCEEN_WIDTH;

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
        q{conda_path} => (
            documentation => q{Conda path},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{container_manager} => (
            cmd_tags      => [q{Default: singularity}],
            documentation => q{Container manager},
            is            => q{rw},
            isa           => enum( [qw{ docker singularity }] ),
        ),
    );

    option(
        q{log_file} => (
            cmd_aliases   => [qw{ log }],
            documentation => q{Log file},
            is            => q{rw},
            isa           => Str,
        )
    );

    ## Special case:Enable/activate MIP. Cannot be changed from cmd or config
    has(
        q{mip} => (
            default => 1,
            is      => q{rw},
            isa     => q{Int},
        )
    );

    option(
        q{verbose} => (
            cmd_aliases   => [qw{ vb }],
            cmd_flag      => q{verbose},
            documentation => q{Turn on chatty output},
            is            => q{rw},
            isa           => Bool,
            required      => 0,
        ),
    );

    return;
}

1;
