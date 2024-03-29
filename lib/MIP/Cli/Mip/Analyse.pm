package MIP::Cli::Mip::Analyse;

use 5.026;
use Carp;
use File::Basename qw{ dirname };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use MooseX::App::Command;
use MooseX::Types::Moose qw{ Str Int HashRef Num Bool ArrayRef };
use Moose::Util::TypeConstraints;

## MIPs lib/
# MooseX::App required sub. Called internally by MooseX::App
use MIP::Cli::Utils qw{ run };

# Set the version for version checking
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
            documentation => q{Set the analysis constant path},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{analysisrunstatus} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation =>
q{Check analysis output and sets the analysis run status flag to finished in sample_info_file},
            is  => q{rw},
            isa => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{analysis_type} => (
            cmd_tags      => [q{sample_id=analysis_type}],
            documentation => q{Type of analysis},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{cluster_constant_path} => (
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
            cmd_tags      => [q{YAML}],
            documentation => q{Write YAML configuration file the analysis parameters},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{core_ram_memory} => (
            cmd_tags      => [q{Default: 5}],
            documentation => q{RAM memory size of the core(s) in GigaBytes},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{dry_run_all} => (
            cmd_aliases   => [qw{ dra }],
            documentation => q{Sets all recipes to dry run mode i.e. no sbatch submission},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{email} => (
            documentation => q{E-mail},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{email_types} => (
            cmd_tags      => [q{Default: FAIL}],
            documentation => q{E-mail type},
            is            => q{rw},
            isa           => ArrayRef [ enum( [qw{ FAIL BEGIN END }] ), ],
        )
    );

    option(
        q{exclude_contigs} => (
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
            cmd_tags      => [q{Default: INFO}],
            documentation => q{Set the GATK log level},
            is            => q{rw},
            isa           => enum( [qw{ DEBUG INFO ERROR FATAL }] ),
        )
    );

    option(
        q{container_config_file} => (
            documentation => q{File with install configuration parameters in YAML format},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{java_use_large_pages} => (
            documentation => q{Use large page memory},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{job_reservation_name} => (
            documentation => q{Allocate node resources from named reservation},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{max_cores_per_node} => (
            cmd_tags      => [q{Default: 16}],
            documentation => q{Maximum number of processor cores per node},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{node_ram_memory} => (
            cmd_tags      => [q{Default: 128}],
            documentation => q{RAM memory size of the node(s) in GigaBytes},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{outdata_dir} => (
            documentation => q{Data output directory},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{outscript_dir} => (
            documentation => q{Script files (.sh) output directory},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{pedigree_file} => (
            cmd_aliases   => [qw{ pedigree }],
            cmd_tags      => [q{YAML}],
            documentation => q{Meta data on samples},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{platform} => (
            cmd_tags      => [q{Default: ILLUMINA}],
            documentation => q{Platform/technology used to produce the reads},
            is            => q{rw},
            isa           => enum( [qw{ ILLUMINA }] ),
        )
    );

    option(
        q{print_recipe} => (
            documentation => q{Print all recipes that are supported},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{print_recipe_mode} => (
            documentation => q{Print all recipes that are supported in [mode]},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{project_id} => (
            documentation => q{Project id},
            is            => q{rw},
            isa           => Str,
        )
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
        q{reference_info_file} => (
            cmd_tags      => [q{YAML}],
            documentation => q{File for reference info used in the analysis},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{sacct} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Generating sbatch script for SLURM info on each submitted job},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
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
        q{sample_ids} => (
            documentation => q{Sample ids},
            is            => q{rw},
            isa           => q{ArrayRef[Str]},
        )
    );

    option(
        q{sample_info_file} => (
            cmd_tags      => [q{YAML}],
            documentation => q{File for sample info used in the analysis},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{set_recipe_core_number} => (
            cmd_tags      => [q{recipe_name=X(cores)}],
            documentation => q{Set the number of cores for specific recipe(s)},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{set_recipe_memory} => (
            cmd_aliases   => [qw{ srm }],
            cmd_tags      => [q{recipe_name=X(G)}],
            documentation => q{Set the memory for specific recipe(s)},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{set_recipe_time} => (
            cmd_tags      => [q{recipe_name=time(hours)}],
            documentation => q{Set the time allocation for specific recipe(s)},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{slurm_quality_of_service} => (
            cmd_aliases   => [qw{ qos }],
            documentation => q{SLURM quality of service},
            is            => q{rw},
            isa           => enum( [qw{ low normal high express }] ),
        )
    );

    option(
        q{start_after_recipe} => (
            cmd_aliases   => [qw{ sar }],
            documentation => q{Start analysis after recipe},
            is            => q{rw},
            isa           => Str,
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
            documentation => q{Submission profile},
            is            => q{rw},
            isa           => enum( [qw{ slurm }] ),
        )
    );

    option(
        q{supported_capture_kit} => (
            cmd_tags      => [q{acronym=file.bed}],
            documentation => q{Set the capture kit acronym shortcut in pedigree file},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{temp_directory} => (
            cmd_aliases   => [qw{ temp }],
            cmd_tags      => [q{Default: "$outdata_dir/$SLURM_JOB_ID"}],
            documentation => q{Set the temporary directory for all recipes},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{version_collect} => (
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Collects executable versions across the analysis},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    return;
}

1;
