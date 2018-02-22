package MIP::Cli::Mip::Analyse::Rare_disease;

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

extends(qw{ MIP::Cli::Mip::Analyse });

command_short_description(q{Rare disease analysis});

command_long_description(
    q{Rare disease analysis on wes, wgs or mixed sequence data});

command_usage(q{mip <analyse> <rare_disease> -pbwa_mem INT});

option(
    q{pbwa_mem} => (
        cmd_aliases   => [qw{ mem }],
        cmd_flag      => q{pbwa_mem},
        cmd_tags      => [q{Analysis recipe switch}],
        documentation => q{Align reads using Bwa Mem (defaults to "0" (=no))},
        is            => q{rw},
        isa           => q{Int},
        required      => 0,
    )
);

sub run {
    my ($arg_href) = @_;

    # do something
    use Data::Dumper;
    say STDERR $arg_href->{pbwa_mem};
    foreach my $sample ( @{ $arg_href->{sample_ids} } ) {
        say STDERR $sample;
    }
    print Dumper($arg_href);
    return;
}

1;
