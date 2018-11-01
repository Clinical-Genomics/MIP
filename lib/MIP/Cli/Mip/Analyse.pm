package MIP::Cli::Mip::Analyse;

use 5.026;
use Carp;
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

our $VERSION = 1.04;

extends(qw{ MIP::Cli::Mip });

command_short_description(q{MIP analyse command});

command_long_description(q{Entry point for performing MIP analysis});

command_usage(q{analyse <pipeline>});

## Define, check and get Cli supplied parameters
_build_usage();

sub run {
    my ($arg_href) = @_;

    say {*STDERR} q{Please choose an subcommand to start the analysis};
    return;
}

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
        q{node_ram_memory} => (
            cmd_aliases   => [qw{ nrm }],
            cmd_tags      => [q{Default: 128}],
            documentation => q{RAM memory size of the node(s) in GigaBytes},
            is            => q{rw},
            isa           => Int,
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
        q{module_core_number} => (
            cmd_aliases   => [qw{ mcn }],
            cmd_tags      => [q{recipe_name=X(cores)}],
            documentation => q{Set the number of cores for each module},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{module_time} => (
            cmd_aliases   => [qw{ mot }],
            cmd_tags      => [q{recipe_name=time(hours)}],
            documentation => q{Set the time allocation for each module},
            is            => q{rw},
            isa           => HashRef,
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
        q{print_recipes} => (
            cmd_aliases   => [qw{ pp }],
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
        q{project_id} => (
            cmd_aliases   => [qw{ pro }],
            documentation => q{Project id},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{reference_dir} => (
            cmd_aliases   => [qw{ rd }],
            documentation => q{Reference(s) directory},
            is            => q{rw},
            isa           => Str,
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
        q{slurm_quality_of_service} => (
            cmd_aliases   => [qw{ qos }],
            cmd_flag      => q{slurm_quly_sri},
            documentation => q{SLURM quality of service},
            is            => q{rw},
            isa           => enum( [qw{ low normal high }] ),
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

    return;
}

1;
