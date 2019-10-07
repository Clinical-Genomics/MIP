package MIP::Cli::Mip::Install;

use 5.026;
use Carp;
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Cli::Utils qw{ run };

our $VERSION = 1.11;

extends(qw{ MIP::Cli::Mip });

command_short_description(q{MIP install command});

command_long_description(q{Entry point for generating MIP installation scripts});

command_usage(q{install <pipeline>});

## Define, check and get Cli supplied parameters
_build_usage();

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

    option(
        q{add_environment_date} => (
            cmd_aliases   => [qw{ aed }],
            cmd_flag      => q{add_environment_date},
            documentation => q{Add creation date to environment},
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
        q{container_dir_path} => (
            cmd_aliases   => [qw{ cntdp }],
            cmd_tags      => [q{Default: /<MIPs_conda_env>/share/containers}],
            documentation => q{Pull containers to this directory },
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{core_number} => (
            cmd_tags      => [q{Default: 1}],
            documentation => q{Number of tasks in sbatch allocation},
            is            => q{rw},
            isa           => Int,
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
        q{environment_prefix} => (
            cmd_aliases => [qw{ ep }],
            documentation =>
              q{Prepend this to environment names. Separated by underscore},
            is  => q{rw},
            isa => Str,
        )
    );

    option(
        q{environment_suffix} => (
            cmd_aliases   => [qw{ es }],
            documentation => q{Append this to environment names. Separated by underscore},
            is            => q{rw},
            isa           => Str,
        )
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
        q{program_test_file} => (
            cmd_aliases   => [qw{ ptf }],
            documentation => q{File with test commands in YAML format},
            is            => q{rw},
            isa           => Str,
        )
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
        q{reference_dir} => (
            cmd_aliases   => [qw{ rd }],
            cmd_tags      => [q{Default: ""}],
            documentation => q{Reference directory},
            is            => q{rw},
            isa           => Str,
        )
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
        q{update_config} => (
            cmd_aliases   => [qw{ uc }],
            cmd_flag      => q{update_config},
            documentation => q{Path to existing config},
            is            => q{rw},
            isa           => Str,
            required      => 0,
        ),
    );

    option(
        q{vep_assemblies} => (
            cmd_aliases   => [qw{ vea }],
            cmd_tags      => [q{Default: GRCh37, GRCh38}],
            cmd_flag      => q{vep_assemblies},
            documentation => q{VEP assemblies to download},
            is            => q{rw},
            isa           => ArrayRef,
            required      => 0,
        ),
    );

    option(
        q{vep_auto_flag} => (
            cmd_aliases   => [qw{ veaf }],
            cmd_flag      => q{vep_auto_flag},
            cmd_tags      => [q{Default: cf}],
            documentation => q{VEP's --AUTO flag},
            is            => q{rw},
            isa           => Str,
            required      => 0,
        ),
    );

    option(
        q{vep_cache_dir} => (
            cmd_aliases => [qw{ vecd }],
            cmd_flag    => q{vep_cache_dir},
            cmd_tags =>
              [q{Default: <reference_dir>/ensembl-tools-release-<version>/cache}],
            documentation => q{VEP's cache directory},
            is            => q{rw},
            isa           => Str,
            required      => 0,
        ),
    );

    option(
        q{vep_plugins} => (
            cmd_aliases   => [qw{ vepl }],
            cmd_flag      => q{vep_plugins},
            documentation => q{VEP plugins to install},
            is            => q{rw},
            isa           => ArrayRef,
            required      => 0,
        ),
    );

    option(
        q{vep_species} => (
            cmd_aliases   => [qw{ ves }],
            cmd_tags      => [q{Default: homo_sapiens_merged}],
            cmd_flag      => q{vep_species},
            documentation => q{VEP species},
            is            => q{rw},
            isa           => ArrayRef,
            required      => 0,
        ),
    );

    option(
        q{write_config} => (
            cmd_aliases   => [qw{ wc }],
            cmd_flag      => q{write_config},
            documentation => q{Generate config from template},
            is            => q{rw},
            isa           => Bool,
            required      => 0,
        ),
    );

    return;
}

1;
