package MIP::Cli::Mip;

use Carp;
use open qw{ :encoding(UTF-8) :std };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use MooseX::App qw{ Color BashCompletion };

## MIPs lib/
use MIP::Cli::Mip;

our $VERSION = 0.01;

## Define, check and get Cli supplied parameters
_build_usage();

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

    option(
        q{config_file} => (
            cmd_aliases => [qw{ config }],
            cmd_flag    => q{config_file},
            documentation =>
              q{YAML config file for analysis parameters (defaults to "")},
            is       => q{rw},
            isa      => q{Str},
            required => 0,
        )
    );

    option(
        q{version} => (
            cmd_flag      => q{v},
            documentation => q{Show version},
            is            => q{rw},
            isa           => q{Bool},
            required      => 0,
        )
    );
    return;
}

1;
