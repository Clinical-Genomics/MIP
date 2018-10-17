package MIP::Cli::Mip::Install;

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
use MooseX::Types::Moose qw{ ArrayRef Bool HashRef Int Str };

our $VERSION = 1.01;

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
        q{conda_no_update_dep} => (
            cmd_aliases   => [qw{ cnud }],
            cmd_flag      => q{conda_no_update_dep},
            documentation => q{Do not update dependencies},
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
        q{email} => (
            cmd_aliases   => [qw{ em }],
            documentation => q{E-mail},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{email_types} => (
            cmd_aliases   => [qw{ emt }],
            cmd_tags      => [q{Default: FAIL}],
            documentation => q{E-mail type},
            is            => q{rw},
            isa           => ArrayRef [ enum( [qw{ FAIL BEGIN END }] ), ],
        )
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
        q{project_id} => (
            cmd_aliases   => [qw{ pro }],
            documentation => q{Project id},
            is            => q{rw},
            isa           => Str,
        )
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

    option(
        q{sbatch_mode} => (
            documentation => q{Write install script for sbatch submisson},
            is            => q{rw},
            isa           => Bool,
            required      => 0,

        ),
    );

    option(
        q{sbatch_process_time} => (
            cmd_aliases   => [qw{ spt }],
            documentation => q{Time limit for sbatch job},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{slurm_quality_of_service} => (
            cmd_aliases   => [qw{ qos }],
            documentation => q{SLURM quality of service},
            is            => q{rw},
            isa           => enum( [qw{ low normal high }] ),
        )
    );

    return;
}

1;
