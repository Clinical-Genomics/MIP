package MIP::Cli::Mip::Analyse::Cancer;

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

command_short_description(q{Cancer analysis});
command_long_description(q{Cancer analysis on panel, wes or wgs sequence data});

command_usage(q{mip <analyse> <cancer> -pbwa_mem INT});

option(
    q{pbwa_mem} => (
        cmd_aliases   => [qw{ mem }],
        cmd_flag      => q{pbwa_mem},
        cmd_tags      => [qw{ Analysis recipe switch! }],
        documentation => q{Align reads using Bwa Mem (defaults to "0" (=no))},
        is            => q{rw},
        isa           => q{Int},
        required      => 0,
    )
);

sub run {
    my ($self) = @_;

    # do something
    say STDERR 'HELLO WORLD';
    return;
}

1;
