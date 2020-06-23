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
use MooseX::Types::Moose qw{ ArrayRef Bool HashRef Int Str };

## MIPs lib/
## MooseX::App required sub. Called internally by MooseX::App
use MIP::Cli::Utils qw{ run };

our $VERSION = 1.05;

extends(qw{ MIP::Cli::Mip });

command_short_description(q{MIP download command});

command_long_description(q{Entry point for generating MIP download references scripts});

command_usage(q{download <pipeline>});

## Define, check and get Cli supplied parameters
_build_usage();

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

    option(
        q{config_file} => (
            cmd_aliases   => [qw{ config c }],
            documentation => q{File with configuration parameters in YAML format},
            is            => q{rw},
            isa           => Str,
            required      => 1,
        )
    );

    option(
        q{dry_run_all} => (
            cmd_aliases => [qw{ dra }],
            documentation =>
              q{Sets all recipes to dry run mode i.e. no sbatch submission},
            is  => q{rw},
            isa => Bool,
        )
    );

    option(
        q{print_parameter_default} => (
            cmd_flag      => q{print_parameter_default},
            documentation => q{print the default parameters},
            is            => q{rw},
            isa           => Bool,
        ),
    );

    option(
        q{project_id} => (
            documentation => q{Project id},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{reference} => (
            cmd_aliases   => [qw{ ref }],
            cmd_flag      => q{reference},
            documentation => q{References to download},
            is            => q{rw},
            isa           => HashRef,
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
        q{reference_genome_versions} => (
            cmd_aliases   => [qw{ rgv }],
            cmd_flag      => q{reference_genome_versions},
            cmd_tags      => [q{Default: grch37, grch38}],
            documentation => q{Reference genomes to download},
            is            => q{rw},
            isa           => ArrayRef [ enum( [qw{ grch37 grch38 }] ), ],
        ),
    );

    option(
        q{sacct_format_fields} => (
            cmd_tags => [
q{Default: jobid, jobname%50, account, partition, alloccpus, TotalCPU, elapsed, start, end, state, exitcode}
            ],
            documentation => q{Format and fields of sacct output},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{submission_profile} => (
            documentation => q{Submission profile},
            is            => q{rw},
            isa           => enum( [qw{ slurm }] ),
        )
    );

    option(
        q{temp_directory} => (
            cmd_aliases   => [qw{ temp }],
            cmd_tags      => [q{Default: "/Current_dir/mip_download/$SLURM_JOB_ID"}],
            documentation => q{Set the temporary directory for all recipes},
            is            => q{rw},
            isa           => Str,
        )
    );

    return;
}

1;
