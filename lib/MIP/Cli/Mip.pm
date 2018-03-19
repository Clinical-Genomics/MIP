package MIP::Cli::Mip;

use Carp;
use open qw{ :encoding(UTF-8) :std };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use MooseX::App qw{ BashCompletion Color MutexGroup Typo };
use MooseX::Types::Moose qw{ Str Int HashRef Num Bool ArrayRef };
use Moose::Util::TypeConstraints;

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
        q{bash_set_errexit} => (
            cmd_aliases   => [qw{ bse }],
            documentation => q{Set errexit in bash scripts},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{bash_set_nounset} => (
            cmd_aliases   => [qw{ bsu }],
            documentation => q{Set nounset in bash scripts},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{bash_set_pipefail} => (
            cmd_aliases   => [qw{ bsp }],
            documentation => q{Set pipefail in bash scripts},
            is            => q{rw},
            isa           => Bool,
        )
    );
    option(
        q{config_file} => (
            cmd_aliases => [qw{ config c }],
            documentation =>
              q{File with configuration parameters in YAML format},
            is       => q{rw},
            isa      => Str,
            required => 1,
        )
    );

    option(
        q{log_file} => (
            cmd_aliases   => [qw{ log }],
            documentation => q{Log file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{temp_directory} => (
            cmd_aliases   => [qw{ tmd }],
            cmd_tags      => [q{Default: "/scratch/$SLURM_JOB_ID"}],
            documentation => q{Set the temporary directory for all programs},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{version} => (
            cmd_aliases   => q{v},
            documentation => q{Show version},
            is            => q{rw},
            isa           => Bool,
        )
    );
    return;
}

1;
