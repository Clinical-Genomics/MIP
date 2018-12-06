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

our $VERSION = 0.03;

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
        q{conda_path} => (
            cmd_aliases   => [qw{ conp }],
            documentation => q{Conda path},
            is            => q{rw},
            isa           => Str,
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

    ## Special case:Enable/activate MIP. Cannot be changed from cmd or config
    has(
        q{mip} => (
            default => 1,
            is      => q{rw},
            isa     => q{Int},
        )
    );

    option(
        q{temp_directory} => (
            cmd_aliases   => [qw{ tmd }],
            cmd_tags      => [q{Default: "/scratch/$SLURM_JOB_ID"}],
            documentation => q{Set the temporary directory for all recipes},
            is            => q{rw},
            isa           => Str,
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
