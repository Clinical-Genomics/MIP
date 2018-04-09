package MIP::Cli::Mip::Install;

use 5.018;
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

our $VERSION = 0.0.5;

extends(qw{ MIP::Cli::Mip });

command_short_description(q{MIP install command});

command_long_description(
    q{Entry point for generating MIP installation scripts});

command_usage(q{install <pipeline>});

## Define, check and get Cli supplied parameters
_build_usage();

sub run {
    my ($arg_href) = @_;

    say {*STDERR}
q{Please choose a pipeline to start generation of the installation script};
    return;
}

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

    option(
        q{conda_dir_path} => (
            cmd_aliases   => [qw{ cdp }],
            cmd_flag      => q{conda_dir_path},
            documentation => q{Path to conda_directory},
            is            => q{rw},
            isa           => Str,
            required      => 0,
        ),
    );

    option(
        q{conda_update} => (
            cmd_aliases   => [qw{ cdu }],
            cmd_flag      => q{conda_update},
            documentation => q{Update conda},
            is            => q{rw},
            isa           => Bool,
            required      => 0,
        ),
    );

    option(
        q{disable_env_check} => (
            cmd_aliases   => [qw{ dec }],
            cmd_flag      => q{disable_env_check},
            documentation => q{Disable source environment check},
            is            => q{rw},
            isa           => Bool,
            required      => 0,
        ),
    );

    option(
        q{noupdate} => (
            cmd_aliases   => [qw{ nup }],
            cmd_flag      => q{noupdate},
            documentation => q{Do not update existing shell programs},
            is            => q{rw},
            isa           => Bool,
            required      => 0,
        ),
    );

    option(
        q{prefer_shell} => (
            cmd_aliases => [qw{ psh }],
            cmd_flag    => q{prefer_shell},
            documentation =>
q{Shell will be used for overlapping shell and biconda installations},
            is       => q{rw},
            isa      => Bool,
            required => 0,
        ),
    );

    option(
        q{print_parameter_default} => (
            cmd_aliases   => [qw{ ppd }],
            cmd_flag      => q{print_parameter_default},
            documentation => q{print the default parameters},
            is            => q{rw},
            isa           => Bool,
            required      => 0,
        ),
    );

    option(
        q{quiet} => (
            cmd_aliases   => [qw{ q }],
            cmd_flag      => q{quiet},
            documentation => q{Limit output from programs},
            is            => q{rw},
            isa           => Bool,
            required      => 0,
        ),
    );

    return;
}

1;
