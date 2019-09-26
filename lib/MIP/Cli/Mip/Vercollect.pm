package MIP::Cli::Mip::Vercollect;

use 5.026;
use Carp;
use Cwd;
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Moose::Util::TypeConstraints;
use MooseX::App::Command;
use MooseX::Types::Moose qw{ ArrayRef Bool HashRef Int Str };

## MIPs lib/
use MIP::Constants qw{ $NEWLINE };

our $VERSION = 1.00;

command_short_description(q{MIP vercollect command});

command_long_description(q{Entry point for collecting MIP execuatable versions});

command_usage(
    q{vercollect <options> -si [sample_info.yaml] -r [regexp.yaml] -o [outfile]});

## Define, check and get Cli supplied parameters
_build_usage();

sub run {

    ## Input from Cli
    my ($arg_href) = @_;

    ## Remove Moose::App extra variable
    delete $arg_href->{extra_argv};

    # Flatten argument(s)
    my $infile_path = $arg_href->{infile_path};
    my $log_file    = $arg_href->{log_file};
    my $outfile     = $arg_href->{outfile};

    use MIP::Log::MIP_log4perl qw{ initiate_logger };
    use MIP::Main::Vercollect qw{ mip_vercollect };

    ## Creates log object
    my $log = initiate_logger(
        {
            file_path => $log_file,
            log_name  => uc q{mip_vercollect},
        }
    );

    mip_vercollect(
        {
            infile_path => $infile_path,
            log         => $log,
            outfile     => $outfile,
        }
    );
    return;
}

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

    my ($arg_href) = @_;

    option(
        q{log_file} => (
            cmd_aliases   => [qw{ l }],
            default       => catfile( cwd(), q{vercollect.log} ),
            documentation => q{Log file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{outfile} => (
            cmd_aliases   => [qw{ o }],
            cmd_tags      => [q{YAML}],
            default       => q{qcmetrics.yaml},
            documentation => q{Data file output},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{infile_path} => (
            cmd_aliases   => [qw{ i inf }],
            cmd_tags      => [q{YAML}],
            documentation => q{Binary file for executables used in the analysis},
            is            => q{rw},
            isa           => Str,
            required      => 1,
        )
    );

    return;
}

1;
