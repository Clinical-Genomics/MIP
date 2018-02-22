package MIP::Cli::Mip;

use Carp;
use open qw{ :encoding(UTF-8) :std };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use MooseX::App qw{ Config Color BashCompletion };

## MIPs lib/
use MIP::Cli::Mip;

our $VERSION = 0.01;

option(
    q{version} => (
        cmd_flag      => q{v},
        documentation => q{Show version},
        is            => q{rw},
        isa           => q{Bool},
        required      => 0,
    )
);

1;
