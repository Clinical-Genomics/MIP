package MIP::Cli::Mip::Install;

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

command_short_description(q{MIP install command});

command_long_description(
    q{Entry point for generating MIP installation scripts});

command_usage(q{install <pipeline>});

## Define, check and get Cli supplied parameters
_build_usage();

sub run {
    my ($arg_href) = @_;

    say STDERR
q{Please choose a pipeline to start generation of the installation script};

    return;
}

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

    return;
}

1;
