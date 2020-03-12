package MIP::Cli::Mip::Install::Rd_rna;

use 5.026;
use Carp;
use File::Spec::Functions qw{ catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use MooseX::App::Command;
use Moose::Util::TypeConstraints;
use MooseX::Types::Moose qw{ Str Int HashRef Bool ArrayRef };
use MooseX::Types::Structured qw{ Dict Optional };

## MIPs lib
use MIP::Definition qw{ get_parameter_from_definition_files };
use MIP::Main::Install qw{ mip_install };
use MIP::Script::Utils qw{ print_parameter_defaults };

our $VERSION = 2.18;

extends(qw{ MIP::Cli::Mip::Install });

command_short_description(q{Generate mip.sh for installation});
command_long_description(
q{Generates an installation script (mip.sh), which is used for installation of the rare disease RNA flavor of the Mutation Identification Pipeline (MIP).}
);
command_usage(q{mip <install> <rd_rna> [options]});

## Define, check and get Cli supplied parameters
_build_usage();

sub run {
    my ($arg_href) = @_;

    ## Remove Moose::App extra variable
    delete $arg_href->{extra_argv};

    ## Input from Cli
    my %active_parameter = %{$arg_href};

    ## %parameter holds all defined parameters for MIP install rd_rna
    ## CLI commands inheritance level
    my %parameter =
      get_parameter_from_definition_files( { level => q{install_rd_rna}, } );

    ## If no config from cmd
    if ( not $active_parameter{config_file} ) {

        ## Use default
        $active_parameter{config_file} =
          catfile( $Bin, qw{ templates mip_install_rd_rna_config_-1.0-.yaml } );
    }

    ## Print parameters from config file and exit
    print_parameter_defaults(
        {
            parameter_href          => \%parameter,
            print_parameter_default => $arg_href->{print_parameter_default},
        }
    );

    ## Start generating the installation script
    mip_install(
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
        q{config_file} => (
            cmd_aliases   => [qw{ config c }],
            documentation => q{File with configuration parameters in YAML format},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{environment_name} => (
            cmd_aliases   => [qw{ envn }],
            cmd_flag      => q{environment_name},
            cmd_tags      => [q{Default: mip_rd_rna }],
            documentation => q{Set environment name},
            is            => q{rw},
            isa           => Str,
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
                        qw{ arriba blobfish bootstrapann fastqc fusion-filter gatk4
                          gffcompare gtf2bed htslib mip_scripts multiqc picard preseq python
                          rseqc salmon sambamba star star-fusion stringtie trim-galore ucsc
                          utilities vep }
                    ]
                ),
            ],
            required => 0,
        ),
    );
    option(
        q{shell_install} => (
            cmd_aliases   => [qw{ si }],
            cmd_flag      => q{shell_install},
            documentation => q{Install supplied programs via shell instead of via conda},
            is            => q{rw},
            isa           => ArrayRef [ enum( [qw{ picard }] ), ],
            required      => 0,
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
                        qw{ arriba blobfish bootstrapann fastqc fusion-filter gatk4
                          gffcompare gtf2bed htslib mip_scripts multiqc picard preseq python
                          rseqc salmon sambamba star star-fusion stringtie trim-galore ucsc
                          utilities vep }
                    ]
                ),
            ],
            required => 0,
        ),
    );

    return;
}

1;
