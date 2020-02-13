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
        q{job_reservation_name} => (
            cmd_aliases   => [qw{ job_res_name }],
            documentation => q{Allocate node resources from named reservation},
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

    option(
        q{max_cores_per_node} => (
            cmd_aliases   => [qw{ mcpn }],
            cmd_tags      => [q{Default: 16}],
            documentation => q{Maximum number of processor cores per node},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{core_ram_memory} => (
            cmd_aliases   => [qw{ crm }],
            cmd_tags      => [q{Default: 5}],
            documentation => q{RAM memory size of the core(s) in GigaBytes},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{node_ram_memory} => (
            cmd_aliases   => [qw{ nrm }],
            cmd_tags      => [q{Default: 128}],
            documentation => q{RAM memory size of the node(s) in GigaBytes},
            is            => q{rw},
            isa           => Int,
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
        q{slurm_quality_of_service} => (
            cmd_aliases   => [qw{ qos }],
            documentation => q{SLURM quality of service},
            is            => q{rw},
            isa           => enum( [qw{ low normal high }] ),
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
