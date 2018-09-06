package MIP::Cli::Mip::Download;

use 5.026;
use Carp;
use open qw{ :encoding(UTF-8) :std };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use MooseX::App::Command;
use Moose::Util::TypeConstraints;
use MooseX::Types::Moose qw{ Str Int HashRef Bool };

our $VERSION = 0.0.1;

extends(qw{ MIP::Cli::Mip });

command_short_description(q{MIP download command});

command_long_description(
    q{Entry point for generating MIP download references scripts});

command_usage(q{download <pipeline>});

## Define, check and get Cli supplied parameters
_build_usage();

sub run {
    my ($arg_href) = @_;

    say {*STDERR}
      q{Please choose a pipeline to start generation of the download script};
    return;
}

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

    option(
        q{print_parameter_default} => (
            cmd_aliases   => [qw{ ppd }],
            cmd_flag      => q{print_parameter_default},
            documentation => q{print the default parameters},
            is            => q{rw},
            isa           => Bool,
        ),
    );

    option(
        q{reference_dir} => (
            cmd_aliases   => [qw{ rd }],
            cmd_flag      => q{reference_dir},
            cmd_tags      => [q{Default: ""}],
            documentation => q{Download references to this dir},
            is            => q{rw},
            isa           => Str,
        ),
    );

    return;
}

1;
