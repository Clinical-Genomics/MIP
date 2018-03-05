package MIP::Cli::Mip::Analyse;

use 5.018;
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

our $VERSION = 0.01;

extends(qw{ MIP::Cli::Mip });

command_short_description(q{MIP analyse command});

command_long_description(q{Entry point for performing MIP analysis});

command_usage(q{analyse <pipeline>});

## Define, check and get Cli supplied parameters
_build_usage();

sub run {
    my ($arg_href) = @_;

    # do something
    say STDERR q{Please choose an subcommand to start the analysis};
    use Data::Dumper;
    say STDERR $arg_href->{pbwa_mem};
    foreach my $sample ( @{ $arg_href->{sample_ids} } ) {
        say STDERR $sample;
    }
    print Dumper($arg_href);
    return;
}

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

    option(
        q{analysis_constant_path} => (
            cmd_aliases   => [qw{ acp }],
            documentation => q{Set the analysis constant path},
            is            => q{rw},
            isa           => q{Str},
        )
    );

    option(
        q{analysis_type} => (
            cmd_aliases   => [qw{ at }],
            cmd_tags      => [q{sample_id=analysis_type}],
            documentation => q{Type of analysis},
            is            => q{rw},
            isa           => q{HashRef},
        )
    );

    option(
        q{cluster_constant_path} => (
            cmd_aliases   => [qw{ ccp }],
            documentation => q{Set the cluster constant path},
            is            => q{rw},
            isa           => q{Str},
        )
    );

    option(
        q{config_file_analysis} => (
            cmd_aliases => [qw{ cfa }],
            cmd_tags    => [q{YAML}],
            documentation =>
              q{Write YAML configuration file the analysis parameters},
            is  => q{rw},
            isa => q{Str},
        )
    );

    option(
        q{dry_run_all} => (
            cmd_aliases => [qw{ dra }],
            documentation =>
              q{Sets all programs to dry run mode i.e. no sbatch submission},
            is  => q{rw},
            isa => q{Bool},
        )
    );

    option(
        q{email} => (
            cmd_aliases   => [qw{ em }],
            documentation => q{E-mail},
            is            => q{rw},
            isa           => q{Str},
        )
    );

    option(
        q{email_types} => (
            cmd_aliases   => [qw{ emt }],
            cmd_tags      => [q{SLURM}],
            documentation => q{E-mail type},
            is            => q{rw},
            isa           => q{ArrayRef},
        )
    );

    option(
        q{expected_coverage} => (
            cmd_aliases   => [qw{ ec }],
            cmd_tags      => [q{sample_id=expected_coverage}],
            documentation => q{Expected mean target coverage for analysis},
            is            => q{rw},
            isa           => q{HashRef},
        )
    );

    option(
        q{family_id} => (
            cmd_aliases   => [qw{ fam }],
            documentation => q{Group id of samples to be compared},
            is            => q{rw},
            isa           => q{Str},
            required      => 1,
        )
    );

    option(
        q{human_genome_reference} => (
            cmd_aliases   => [qw{ hgr }],
            cmd_tags      => [q{Fasta}],
            documentation => q{Human genome reference},
            is            => q{rw},
            isa           => q{Str},
        )
    );

    option(
        q{infile_dirs} => (
            cmd_aliases   => [qw{ ifd }],
            cmd_tags      => [q{infile_dirs=sample_id}],
            documentation => q{Infile directory(s)},
            is            => q{rw},
            isa           => q{HashRef},
        )
    );

    option(
        q{java_use_large_pages} => (
            cmd_aliases   => [qw{ jul }],
            documentation => q{Use large page memory},
            is            => q{rw},
            isa           => q{Bool},
        )
    );

    option(
        q{module_core_number} => (
            cmd_aliases   => [qw{ mcn }],
            cmd_tags      => [q{program_name=X(cores)}],
            documentation => q{Set the number of cores for each module},
            is            => q{rw},
            isa           => q{HashRef},
        )
    );

    option(
        q{module_time} => (
            cmd_aliases   => [qw{ mot }],
            cmd_tags      => [q{program_name=time(hours)}],
            documentation => q{Set the time allocation for each module},
            is            => q{rw},
            isa           => q{HashRef},
        )
    );

    option(
        q{outaligner_dir} => (
            cmd_aliases => [qw{ ald }],
            documentation =>
q{Sets which aligner out directory was used for alignment in previous analysis},
            is  => q{rw},
            isa => enum( [qw{ bwa star }] ),
        )
    );

    option(
        q{outdata_dir} => (
            cmd_aliases   => [qw{ odd }],
            documentation => q{Data output directory},
            is            => q{rw},
            isa           => q{Str},
        )
    );

    option(
        q{outscript_dir} => (
            cmd_aliases   => [qw{ osd }],
            documentation => q{Script files (.sh) output directory},
            is            => q{rw},
            isa           => q{Str},
        )
    );

    option(
        q{pedigree_file} => (
            cmd_aliases   => [qw{ ped }],
            cmd_tags      => [q{YAML}],
            documentation => q{Meta data on samples},
            is            => q{rw},
            isa           => q{Str},
        )
    );

    option(
        q{platform} => (
            cmd_aliases   => [qw{ pla }],
            documentation => q{Platform/technology used to produce the reads},
            is            => q{rw},
            isa           => enum( [qw{ ILLUMINA }] ),
        )
    );

    option(
        q{print_program_mode} => (
            cmd_aliases   => [qw{ ppm }],
            documentation => q{Print all programs that are supported},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{project_id} => (
            cmd_aliases   => [qw{ pro }],
            documentation => q{Project id},
            is            => q{rw},
            isa           => q{Str},
        )
    );

    option(
        q{print_programs} => (
            cmd_aliases   => [qw{ pp }],
            documentation => q{Print all programs that are supported},
            is            => q{rw},
            isa           => q{Bool},
        )
    );

    option(
        q{reference_dir} => (
            cmd_aliases   => [qw{ rd }],
            documentation => q{Reference(s) directory},
            is            => q{rw},
            isa           => q{Str},
        )
    );

    option(
        q{reduce_io} => (
            cmd_aliases   => [qw{ rio }],
            documentation => q{Run consecutive models at node},
            is            => q{rw},
            isa           => q{Bool},
        )
    );

    option(
        q{replace_iupac} => (
            cmd_aliases => [qw{ riu }],
            documentation =>
              q{Replace IUPAC code in alternative alleles with N},
            is  => q{rw},
            isa => q{Bool},
        )
    );

    option(
        q{sample_info_file} => (
            cmd_aliases   => [qw{ sif }],
            cmd_tags      => [q{YAML}],
            documentation => q{File for sample info used in the analysis},
            is            => q{rw},
            isa           => q{Str},
        )
    );

    option(
        q{start_with_program} => (
            cmd_aliases   => [qw{ swp }],
            documentation => q{Start analysis with program},
            is            => q{rw},
            isa           => q{Str},
        )
    );

    option(
        q{supported_capture_kit} => (
            cmd_aliases => [qw{ sck }],
            cmd_tags    => [q{acronym=file.bed}],
            documentation =>
              q{Set the capture kit acronym shortcut in pedigree file},
            is  => q{rw},
            isa => q{HashRef},
        )
    );

    option(
        q{sample_ids} => (
            cmd_aliases => [qw{ spi }],
            documentation =>
              q{Sets all programs to dry run mode i.e. no sbatch submission},
            is  => q{rw},
            isa => q{ArrayRef},
        )
    );

    option(
        q{pbwa_mem} => (
            cmd_aliases   => [qw{ pmem }],
            cmd_flag      => q{pbwa_mem},
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Align reads using Bwa Mem},
            is            => q{rw},
            isa           => enum( [ 1, 2 ] ),
        )
    );

    return;
}

1;
