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
        q{analysis_type} => (
            cmd_aliases => [qw{ at }],
            cmd_flag    => q{analysis_type},
            documentation =>
q{Type of analysis (defaults to 'wgs'; Valid entries: 'wgs', 'wes', 'wts', 'cancer'; sample_id=analysis_type)},
            is       => q{rw},
            isa      => q{HashRef},
        )
    );

    option(
        q{dry_run_all} => (
            cmd_aliases => [qw{ dra }],
            cmd_flag    => q{dry_run_all},
            documentation =>
              q{Sets all programs to dry run mode i.e. no sbatch submission},
            is       => q{rw},
            isa      => q{Bool},
        )
    );

    option(
        q{config_file_analysis} => (
            cmd_aliases => [qw{ cfa }],
            cmd_flag    => q{config_file_analysis},
            documentation =>
              q{Write YAML configuration file the analysis parameters (defaults to "")},
            is       => q{rw},
            isa      => q{Str},
        )
    );

    option(
        q{family_id} => (
            cmd_aliases => [qw{ fam }],
            cmd_flag    => q{family_id},
            documentation =>
              q{Group id of samples to be compared (defaults to "")},
            is       => q{rw},
            isa      => q{Str},
            required => 1,
        )
    );

    option(
        q{infile_dirs} => (
            cmd_aliases   => [qw{ ifd }],
            cmd_flag      => q{infile_dirs},
            documentation => q{Infile directory(s); infile_dirs=sample_id},
            is            => q{rw},
            isa           => q{HashRef},
        )
    );

    option(
        q{outdata_dir} => (
            cmd_aliases => [qw{ odd }],
            cmd_flag    => q{outdata_dir},
            documentation =>
              q{Data output directory},
            is       => q{rw},
            isa      => q{Str},
        )
    );

    option(
        q{outscript_dir} => (
            cmd_aliases => [qw{ osd }],
            cmd_flag    => q{outscript_dir},
            documentation =>
              q{Script files (.sh) output directory},
            is       => q{rw},
            isa      => q{Str},
        )
    );

    option(
        q{project_id} => (
            cmd_aliases => [qw{ pro }],
            cmd_flag    => q{project_id},
            documentation =>
              q{Project id},
            is       => q{rw},
            isa      => q{Str},
        )
    );

    option(
        q{reference_dir} => (
            cmd_aliases => [qw{ rd }],
            cmd_flag    => q{reference_dir},
            documentation =>
              q{Reference(s) directory},
            is       => q{rw},
            isa      => q{Str},
        )
    );

    option(
        q{supported_capture_kit} => (
            cmd_aliases   => [qw{ sck }],
            cmd_flag      => q{supported_capture_kit},
            documentation => q{Set the capture kit acronym shortcut in pedigree file},
            is            => q{rw},
            isa           => q{HashRef},
        )
    );

    option(
        q{sample_ids} => (
            cmd_aliases => [qw{ spi }],
            cmd_flag    => q{sample_ids},
            documentation =>
              q{Sets all programs to dry run mode i.e. no sbatch submission},
            is       => q{rw},
            isa      => q{ArrayRef},
        )
    );
    return;
}

1;
