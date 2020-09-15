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
use MIP::Definition qw{ get_parameter_from_definition_files };
use MIP::Main::Download qw{ mip_download };

our $VERSION = 1.05;

extends(qw{ MIP::Cli::Mip });

command_short_description(q{Generate bash or sbatch for download of references});
command_long_description(
q{Generates download script(s), which is used for downloading reference(s) for the Mutation Identification Pipeline (MIP).}
);
command_usage(q{mip <download> [options]});

## Define, check and get Cli supplied parameters
_build_usage();

sub run {

    my ($arg_href) = @_;

    ## Remove Moose::App extra variable
    delete $arg_href->{extra_argv};

    ## Input from Cli
    my %active_parameter = %{$arg_href};

    ## %parameter holds all defined parameters for MIP download rd_dna
    ## CLI commands inheritance level
    my %parameter = get_parameter_from_definition_files( { level => q{download}, } );

    ## Start generating the installation script
    mip_download(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
        }
    );
    return;
}

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

    has(
        q{download_pipeline} => (
            default => 1,
            is      => q{ro},
            isa     => Bool,
        )
    );

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
        q{environment_name} => (
            cmd_aliases   => [qw{ envn }],
            cmd_flag      => q{environment_name},
            cmd_tags      => [q{Default: mip_rd_dna}],
            documentation => q{Set environment name},
            is            => q{rw},
            isa           => Str,
        ),
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
