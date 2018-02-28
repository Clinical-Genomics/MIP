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

command_usage(q{mip <analyse> <cancer> --config <config_file> --fam <family_id>});

## Define, check and get Cli supplied parameters
_build_usage();

sub run {
    my ($self) = @_;

    # do something
    say STDERR 'HELLO WORLD';
    return;
}

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

option(
        q{sample_origin} => (
            cmd_aliases   => [qw{ sao }],
            cmd_tags      => [q{sample_id=sample_origin}],
            documentation => q{Sample origin for analysis},
            is            => q{rw},
            isa           => q{HashRef},
        )
    );

return;
}

1;
