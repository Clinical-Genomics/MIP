package MIP::Cli::Mip::Analyse;

use 5.026;
use Carp;
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use MooseX::App::Command;
use MooseX::Types::Moose qw{ Str Int HashRef Num Bool ArrayRef };
use Moose::Util::TypeConstraints;

## MIPs lib/
use MIP::Cli::Utils qw{ run }
  ;    # MooseX::App required sub. Called internally by MooseX::App

# Set the version for version checking
our $VERSION = 1.09;

extends(qw{ MIP::Cli::Mip });

command_short_description(q{MIP analyse command});

command_long_description(q{Entry point for performing MIP analysis});

command_usage(q{analyse <pipeline>});

## Define, check and get Cli supplied parameters
_build_usage();

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

    my ($arg_href) = @_;

    option(
        q{analysis_constant_path} => (
            cmd_aliases   => [qw{ acp }],
            documentation => q{Set the analysis constant path},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{analysisrunstatus} => (
            cmd_aliases => [qw{ ars }],
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
q{Check analysis output and sets the analysis run status flag to finished in sample_info_file},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{analysis_type} => (
            cmd_aliases   => [qw{ at }],
            cmd_tags      => [q{sample_id=analysis_type}],
            documentation => q{Type of analysis},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{cluster_constant_path} => (
            cmd_aliases   => [qw{ ccp }],
            documentation => q{Set the cluster constant path},
            is            => q{rw},
            isa           => Str,
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
        q{config_file_analysis} => (
            cmd_aliases   => [qw{ cfa }],
            cmd_tags      => [q{YAML}],
            documentation => q{Write YAML configuration file the analysis parameters},
            is            => q{rw},
            isa           => Str,
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
        q{exclude_contigs} => (
            cmd_aliases   => [qw{ exc }],
            documentation => q{Exclude contigs from analysis},
            is            => q{rw},
            isa           => ArrayRef,
        )
    );

    parameter(
        q{case_id} => (
            documentation => q{Group id of samples to be compared},
            is            => q{rw},
            isa           => Str,
            required      => 1,
        )
    );

    option(
        q{gatk_logging_level} => (
            cmd_aliases   => [qw{ gll }],
            cmd_tags      => [q{Default: INFO}],
            documentation => q{Set the GATK log level},
            is            => q{rw},
            isa           => enum( [qw{ DEBUG INFO ERROR FATAL }] ),
        )
    );

    option(
        q{gatk_path} => (
            cmd_aliases   => [qw{ gtp }],
            documentation => q{Path to GATK},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{java_use_large_pages} => (
            cmd_aliases   => [qw{ jul }],
            documentation => q{Use large page memory},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{outdata_dir} => (
            cmd_aliases   => [qw{ odd }],
            documentation => q{Data output directory},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{outscript_dir} => (
            cmd_aliases   => [qw{ osd }],
            documentation => q{Script files (.sh) output directory},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{pedigree_file} => (
            cmd_aliases   => [qw{ ped }],
            cmd_tags      => [q{YAML}],
            documentation => q{Meta data on samples},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{platform} => (
            cmd_aliases   => [qw{ pla }],
            cmd_tags      => [q{Default: ILLUMINA}],
            documentation => q{Platform/technology used to produce the reads},
            is            => q{rw},
            isa           => enum( [qw{ ILLUMINA }] ),
        )
    );

    option(
        q{print_recipe} => (
            cmd_aliases   => [qw{ pr }],
            documentation => q{Print all recipes that are supported},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{print_recipe_mode} => (
            cmd_aliases   => [qw{ prm }],
            documentation => q{Print all recipes that are supported in [mode]},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{sacct} => (
            cmd_aliases => [qw{ sac }],
            cmd_tags    => [q{Analysis recipe switch}],
            documentation =>
              q{Generating sbatch script for SLURM info on each submitted job},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
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
        q{sample_ids} => (
            cmd_aliases   => [qw{ spi }],
            documentation => q{Sample ids},
            is            => q{rw},
            isa           => q{ArrayRef[Str]},
        )
    );

    option(
        q{sample_info_file} => (
            cmd_aliases   => [qw{ sif }],
            cmd_tags      => [q{YAML}],
            documentation => q{File for sample info used in the analysis},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{singularity_container} => (
            cmd_aliases   => [qw{ sic }],
            cmd_tags      => [q{executable=path/to/container}],
            documentation => q{Specify path to singularity container for executable},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{start_with_recipe} => (
            cmd_aliases   => [qw{ swr }],
            documentation => q{Start analysis with recipe},
            is            => q{rw},
            isa           => Str,
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

    option(
        q{supported_capture_kit} => (
            cmd_aliases   => [qw{ sck }],
            cmd_tags      => [q{acronym=file.bed}],
            documentation => q{Set the capture kit acronym shortcut in pedigree file},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{with_singularity} => (
            cmd_aliases => [qw{ wsi }],
            documentation =>
              q{Run programs inside a singularity container where available},
            is  => q{rw},
            isa => Bool,
        )
    );

    return;
}

1;
