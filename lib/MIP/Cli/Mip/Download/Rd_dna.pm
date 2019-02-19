package MIP::Cli::Mip::Download::Rd_dna;

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

our $VERSION = 1.02;

extends(qw{ MIP::Cli::Mip::Download });

command_short_description(q{Generate bash or sbatch for download of references});
command_long_description(
q{Generates download script(s), which is used for downloading reference(s) for the rare disease DNA flavor of the Mutation Identification Pipeline (MIP).}
);
command_usage(q{mip <download> <rd_dna> [options]});

## Define, check and get Cli supplied parameters
_build_usage();

sub run {

    my ($arg_href) = @_;

    ## Remove Moose::App extra variable
    delete $arg_href->{extra_argv};

    ## Input from Cli
    my %active_parameter = %{$arg_href};

    use MIP::File::Format::Parameter qw{ parse_definition_file  };

    ## Mip analyse rd_dna parameters
    ## CLI commands inheritance
    my @definition_files = (
        catfile( $Bin, qw{ definitions mip_parameters.yaml } ),
        catfile( $Bin, qw{ definitions download_parameters.yaml } ),
        catfile( $Bin, qw{ definitions download_rd_dna_parameters.yaml } ),
    );

    ## Non mandatory parameter definition keys to check
    my $non_mandatory_parameter_keys_path =
      catfile( $Bin, qw{ definitions non_mandatory_parameter_keys.yaml } );

    ## Mandatory parameter definition keys to check
    my $mandatory_parameter_keys_path =
      catfile( $Bin, qw{ definitions mandatory_parameter_keys.yaml } );

    ## %parameter holds all defined parameters for MIP
    ## mip download rd_dna parameters
    my %parameter;

  DEFINITION_FILE:
    foreach my $definition_file (@definition_files) {

        %parameter = (
            %parameter,
            parse_definition_file(
                {
                    define_parameters_path        => $definition_file,
                    mandatory_parameter_keys_path => $mandatory_parameter_keys_path,
                    non_mandatory_parameter_keys_path =>
                      $non_mandatory_parameter_keys_path,
                }
            ),
        );
    }

    ## Print parameters from config file and exit
    if ( $arg_href->{print_parameter_default} ) {
        print_parameter_defaults(
            {
                parameter_href => \%parameter,
            }
        );
    }

    ## Start generating the installation script
    mip_download(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
        }
    );
    return;
}

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

    option(
        q{reference} => (
            cmd_aliases   => [qw{ ref }],
            cmd_flag      => q{reference},
            documentation => q{References to download},
            is            => q{rw},
            isa           => HashRef,
        ),
    );

    return;
}

1;
