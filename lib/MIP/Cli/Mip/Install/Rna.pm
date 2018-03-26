package MIP::Cli::Mip::Install::Rna;

use 5.022;
use Carp;
use Cwd qw{ abs_path };
use File::Spec::Functions qw{ catdir catfile };
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
use MooseX::Types::Structured qw{ Dict Optional };
use Readonly;

## MIPs lib
use MIP::File::Format::Yaml qw{ load_yaml };
use MIP::Main::Install qw{ mip_install };
use MIP::Script::Utils qw{ nest_hash print_install_defaults };

our $VERSION = '0.0.4';

extends(qw{ MIP::Cli::Mip::Install });

command_short_description(q{Generate mip.sh for installation});
command_long_description(
q{Generates an installation script (mip.sh), which is used for installation of the RNA flavor of the Mutation Identification Pipeline (MIP).}
);
command_usage(q{mip <install> <rare_disease> [options]});

## Constants
Readonly my $SPACE => q{ };

## Define, check and get Cli supplied parameters
_build_usage();

sub run {
    my ($arg_href) = @_;

    ## Remove Moose::App extra variable
    delete $arg_href->{extra_argv};

    ## Load default parameters from config file
    my $install_parameters_path = abs_path( $arg_href->config_file );
    my %parameter               = load_yaml(
        {
            yaml_file => $install_parameters_path
        }
    );

    ## Print parameters from config file and exit
    if ( $arg_href->print_parameter_default ) {
        print_install_defaults(
            {
                parameter_href => \%parameter,
            }
        );
    }

    ## Merge arrays and overwrite flat values in config YAML with command line
    @parameter{ keys %{$arg_href} } = values %{$arg_href};

    ## Nest the command line parameters and overwrite the default
    nest_hash( { cmd_href => \%parameter } );

    ## Start generating the installation script
    mip_install(
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
        q{bioconda_programs} => (
            cmd_aliases   => [qw{ bc }],
            cmd_flag      => q{bioconda},
            documentation => q{Set bioconda version of programs},
            is            => q{rw},
            isa           => Dict [
                cufflinks => Optional [Str],
                fastqc    => Optional [Str],
                htslib    => Optional [Str],
                picard    => Optional [Str],
                salmon    => Optional [Str],
                samtools  => Optional [Str],
                star      => Optional [Str],
            ],
            required => 0,
        ),
    );

    option(
        q{conda_packages} => (
            cmd_aliases   => [qw{ cpa }],
            cmd_flag      => q{conda_packages},
            cmd_tags      => [q{Default: pip, python=2.7}],
            documentation => q{Base conda packages that are always installed},
            is            => q{rw},
            isa           => HashRef,
            required      => 0,
        ),
    );

    option(
        q{pip} => (
            cmd_aliases   => [qw{ pip }],
            cmd_flag      => q{pip_programs},
            documentation => q{Set the version of programs installed via pip},
            is            => q{rw},
            isa           => HashRef,
            required      => 0,
        ),
    );

    option(
        q{conda_packages:python} => (
            cmd_aliases   => [qw{ pyv }],
            cmd_flag      => q{python_version},
            cmd_tags      => [q{Default: 2.7}],
            documentation => q{Python version to install},
            is            => q{rw},
            isa           => Str,
            required      => 0,
        ),
    );

    option(
        q{reference_dir} => (
            cmd_aliases   => [qw{ rd }],
            cmd_flag      => q{reference_dir},
            cmd_tags      => [q{Default: ""}],
            documentation => q{Install references to this dir},
            is            => q{rw},
            isa           => Str,
            required      => 0,
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
            required      => 0,
        ),
    );

    option(
        q{select_programs} => (
            cmd_aliases   => [qw{ sp }],
            cmd_flag      => q{select_programs},
            documentation => q{Install only selected programs},
            is            => q{rw},
            isa           => ArrayRef [
                enum(
                    [
                        qw{ cufflinks fastqc htslib mip_scripts picard salmon
                          samtools star star_fusion }
                    ]
                ),
            ],
            required => 0,
        ),
    );
    option(
        q{shell_install} => (
            cmd_aliases => [qw{ si }],
            cmd_flag    => q{shell_install},
            documentation =>
              q{Install supplied programs via shell instead of via conda},
            is       => q{rw},
            isa      => ArrayRef [ enum( [qw{ picard }] ), ],
            required => 0,
        ),
    );

    option(
        q{shell_programs} => (
            cmd_aliases   => [qw{ shell }],
            cmd_flag      => q{shell_programs},
            documentation => q{Set shell version of programs},
            is            => q{rw},
            isa           => Dict [
                picard      => Optional [Str],
                star_fusion => Optional [Str],
            ],
            required => 0,
        ),
    );

    option(
        q{skip_programs} => (
            cmd_aliases   => [qw{ skip }],
            cmd_flag      => q{skip_programs},
            documentation => q{Disable installation of supplied programs},
            is            => q{rw},
            isa           => ArrayRef [
                enum(
                    [
                        qw{ cufflinks fastqc htslib mip_scripts picard salmon
                          samtools star star_fusion }
                    ]
                ),
            ],
            required => 0,
        ),
    );

    return;
}

1;
