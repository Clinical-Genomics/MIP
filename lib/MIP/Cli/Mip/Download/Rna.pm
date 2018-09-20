package MIP::Cli::Mip::Download::Rna;

use 5.026;
use Carp;
use Cwd qw{ abs_path };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use List::Util qw{ any };
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
use MooseX::Types::Moose qw{ Str Int HashRef Bool ArrayRef };
use Readonly;

## MIPs lib
use MIP::File::Format::Yaml qw{ load_yaml };
use MIP::Main::Download qw{ mip_download };
use MIP::Script::Utils qw{ print_parameter_defaults };

our $VERSION = 0.01;

extends(qw{ MIP::Cli::Mip::Download });

command_short_description(q{Generate mip.sh for download of references});
command_long_description(
q{Generates a download script (download_reference.sh), which is used for downloading reference(s) for the rna flavor of the Mutation Identification Pipeline (MIP).}
);
command_usage(q{mip <download> <rna> [options]});

## Define, check and get Cli supplied parameters
_build_usage();

sub run {
    my ($arg_href) = @_;

    ## Remove Moose::App extra variable
    delete $arg_href->{extra_argv};

    ## Load default parameters from config file
    my $download_parameters_path = abs_path( $arg_href->config_file );
    my %parameter                = load_yaml(
        {
            yaml_file => $download_parameters_path,
        }
    );

    ## Print parameters from config file and exit
    if ( $arg_href->{print_parameter_default} ) {
        print_parameter_defaults(
            {
                parameter_href => \%parameter,
            }
        );
    }

    ## Merge arrays and overwrite flat values in config YAML with command line
    @parameter{ keys %{$arg_href} } = values %{$arg_href};

    ## Start generating the installation script
    mip_download(
        {
            parameter_href => \%parameter,
        }
    );
    return;
}

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

    option(
        q{config_file} => (
            cmd_aliases => [qw{ config c }],
            documentation =>
              q{File with configuration parameters in YAML format},
            is      => q{rw},
            isa     => Str,
            default => catfile(
                dirname($Bin),
                qw{ MIP definitions download_rna_parameters.yaml }
            ),
        )
    );

    option(
        q{cmd_reference} => (
            cmd_aliases   => [qw{ ref }],
            cmd_flag      => q{reference},
            documentation => q{References to download},
            is            => q{rw},
            isa           => HashRef,
        ),
    );

    option(
        q{reference_genome_versions} => (
            cmd_aliases   => [qw{ rg }],
            cmd_flag      => q{reference_genome_versions},
            cmd_tags      => [q{Default: GRCh37, hg38}],
            documentation => q{Reference genomes to download},
            is            => q{rw},
            isa           => ArrayRef [ enum( [qw{ GRCh37 hg38 }] ), ],
        ),
    );

    return;
}

1;
