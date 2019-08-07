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
use MIP::Cli::Utils qw{ run }
  ;    # MooseX::App required sub. Called internally by MooseX::App

our $VERSION = 1.01;

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
            cmd_aliases   => [qw{ ppd }],
            cmd_flag      => q{print_parameter_default},
            documentation => q{print the default parameters},
            is            => q{rw},
            isa           => Bool,
        ),
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
        q{reference_genome_versions} => (
            cmd_aliases   => [qw{ rg }],
            cmd_flag      => q{reference_genome_versions},
            cmd_tags      => [q{Default: grch37, grch38}],
            documentation => q{Reference genomes to download},
            is            => q{rw},
            isa           => ArrayRef [ enum( [qw{ grch37 grch38 }] ), ],
        ),
    );

    option(
        q{sacct_format_fields} => (
            cmd_aliases => [qw{ sacfrf }],
            cmd_tags    => [
q{Default: jobid, jobname%50, account, partition, alloccpus, TotalCPU, elapsed, start, end, state, exitcode}
            ],
            documentation => q{Format and fields of sacct output},
            is            => q{rw},
            isa           => ArrayRef [Str],
        )
    );

    option(
        q{submission_profile} => (
            cmd_aliases   => [qw{ sbp }],
            documentation => q{Submission profile},
            is            => q{rw},
            isa           => enum( [qw{ slurm }] ),
        )
    );

    return;
}

1;
