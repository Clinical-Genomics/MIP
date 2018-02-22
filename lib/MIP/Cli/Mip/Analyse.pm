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

option(
    q{analysis_type} => (
        cmd_aliases => [qw{ at }],
        cmd_flag    => q{analysis_type},
        documentation =>
q{Type of analysis (defaults to 'wgs'; Valid entries: 'wgs', 'wes', 'wts', 'cancer'; sample_id=analysis_type)},
        is       => q{rw},
        isa      => q{HashRef},
        required => 0,
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
        required => 0,
    )
);

option(
    q{family_id} => (
        cmd_aliases   => [qw{ fam }],
        cmd_flag      => q{family_id},
        documentation => q{Group id of samples to be compared (defaults to "")},
        is            => q{rw},
        isa           => q{Str},
        required      => 1,
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
        required => 0,
    )
);

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

1;
