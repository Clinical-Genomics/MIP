package MIP::Cli::Mip::Install::Rd_dna;

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
use MIP::Get::Parameter qw{ get_install_parameter_attribute };
use MIP::Main::Install qw{ mip_install };
use MIP::Script::Utils qw{ print_parameter_defaults };

our $VERSION = 2.21;

extends(qw{ MIP::Cli::Mip::Install });

command_short_description(q{Generate mip.sh for installation});
command_long_description(
q{Generates an installation script (mip.sh), which is used for installation of the Mutation Identification Pipeline (MIP) for rare diseases DNA.}
);
command_usage(q{mip <install> <rd_dna> [options]});

## Define, check and get Cli supplied parameters
_build_usage();

sub run {
    my ($arg_href) = @_;

    ## Remove Moose::App extra variable
    delete $arg_href->{extra_argv};

    ## Input from Cli
    my %active_parameter = %{$arg_href};

    ## %parameter holds all defined parameters for MIP install rd_dna
    ## CLI commands inheritance level
    my %parameter =
      get_parameter_from_definition_files( { level => q{install_rd_dna}, } );

    ## If no config from cmd
    if ( not $active_parameter{config_file} ) {

        ## Use default
        $active_parameter{config_file} =
          catfile( $Bin, qw{ templates mip_install_rd_dna_config_-1.0-.yaml } );
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
            cmd_tags      => [q{Default: mip_rd_dna}],
            documentation => q{Set environment name},
            is            => q{rw},
            isa           => Str,
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
                        qw{ bedtools bwa bwakit cadd chanjo chromograph
                          cnvnator delly expansionhunter fastqc
                          gatk gatk4 genmod gcc htslib libxml2 libxslt
                          manta mip_scripts multiqc numpy peddy picard pip
                          plink python rhocall rtg-tools sambamba
                          smncopynumbercaller scikit-learn stranger svdb tiddit ucsc upd
                          utilities varg variant_integrity vcf2cytosure vcfanno vep vt }
                    ]
                ),
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
                        qw{ bedtools bwa bwakit cadd chanjo chromograph
                          cnvnator delly expansionhunter fastqc
                          gatk gatk4 genmod gcc htslib libxml2 libxslt
                          manta mip_scripts multiqc numpy peddy picard pigz pip
                          plink python rhocall rtg-tools sambamba
                          smncopynumbercaller scikit-learn stranger svdb tiddit ucsc upd
                          utilities varg variant_integrity vcf2cytosure vcfanno vep vt }
                    ]
                ),
            ],
            required => 0,
        ),
    );

    return;
}

1;
