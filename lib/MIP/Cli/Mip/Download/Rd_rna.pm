package MIP::Cli::Mip::Download::Rd_rna;

use 5.026;
use Carp;
use open qw{ :encoding(UTF-8) :std };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };
use Params::Check qw{ check allow last_error };

## CPANM
use autodie qw{ :all };
use MooseX::App::Command;
use Moose::Util::TypeConstraints;
use MooseX::Types::Moose qw{ ArrayRef Bool HashRef Int Str };

## MIPs lib
use MIP::Definition qw{ get_parameter_from_definition_files };
use MIP::Main::Download qw{ mip_download };
use MIP::Script::Utils qw{ print_parameter_defaults };

our $VERSION = 1.04;

extends(qw{ MIP::Cli::Mip::Download });

command_short_description(q{Generate bash or sbatch for download of references});
command_long_description(
q{Generates a download script(s), which is used for downloading reference(s) for the rare disease RNA flavor of the Mutation Identification Pipeline (MIP).}
);
command_usage(q{mip <download> <rd_rna> [options]});

## Define, check and get Cli supplied parameters
_build_usage();

sub run {

    my ($arg_href) = @_;

    ## Remove Moose::App extra variable
    delete $arg_href->{extra_argv};

    ## Input from Cli
    my %active_parameter = %{$arg_href};

    ## %parameter holds all defined parameters for MIP download rd_rna
    ## CLI commands inheritance level
    my %parameter =
      get_parameter_from_definition_files( { level => q{download_rd_rna}, } );

    ## Print parameters from config file and exit
    print_parameter_defaults(
        {
            parameter_href          => \%parameter,
            print_parameter_default => $arg_href->{print_parameter_default},
        }
    );

    ## Start generating the installation script
    mip_download(
        {
            active_parameter_href => \%active_parameter,

            parameter_href => \%parameter,
        }
    );
    return;
}

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

## Special case:Set download_pipeline_type. Cannot be changed from cmd or config
    has(
        q{download_pipeline_type} => (
            default => q{rd_rna},
            is      => q{rw},
            isa     => Str,
        )
    );
    return;
}

1;
