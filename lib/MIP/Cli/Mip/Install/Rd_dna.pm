package MIP::Cli::Mip::Install::Rd_dna;

use 5.026;
use Carp;
use Cwd qw{ abs_path };
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
use MooseX::Types::Structured qw{ Dict Optional };
use Readonly;

## MIPs lib
use MIP::File::Format::Parameter qw{ parse_definition_file  };
use MIP::File::Format::Yaml qw{ load_yaml };
use MIP::Get::Parameter qw{ get_install_parameter_attribute };
use MIP::Main::Install qw{ mip_install };
use MIP::Script::Utils qw{ print_parameter_defaults };

our $VERSION = 2.11;

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

    ## Mip analyse rd_dna parameters
    ## CLI commands inheritance
    my @definition_files = (
        catfile( $Bin, qw{ definitions mip_parameters.yaml } ),
        catfile( $Bin, qw{ definitions install_parameters.yaml } ),
        catfile( $Bin, qw{ definitions install_rd_dna_parameters.yaml } ),
    );

    ## Non mandatory parameter definition keys to check
    my $non_mandatory_parameter_keys_path =
      catfile( $Bin, qw{ definitions non_mandatory_parameter_keys.yaml } );

    ## Mandatory parameter definition keys to check
    my $mandatory_parameter_keys_path =
      catfile( $Bin, qw{ definitions mandatory_parameter_keys.yaml } );

    ## %parameter holds all defined parameters for MIP
    ## mip install rd_dna parameters
    my %parameter;

    ## If no config from cmd
    if ( not $active_parameter{config_file} ) {

        ## Use default
        $active_parameter{config_file} =
          catfile( $Bin, qw{ templates mip_install_rd_dna_config_-1.0-.yaml } );
    }

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
            cmd_aliases => [qw{ envn }],
            cmd_flag    => q{environment_name},
            cmd_tags    => [
q{Default: mip7_rd-dna mip7_rd-dna_perl5 mip7_rd-dna_py3 mip7_rd-dna_tiddit}
            ],
            documentation => q{Set environment names},
            is            => q{rw},
            isa           => Dict [
                emip   => Optional [Str],
                eperl5 => Optional [Str],
            ],
            required => 0,
        ),
    );

    option(
        q{installations} => (
            cmd_aliases   => [qw{ install }],
            cmd_flag      => q{installations},
            cmd_tags      => [q{Default: emip eperl5}],
            documentation => q{Environments to install},
            is            => q{rw},
            isa           => ArrayRef [ enum( [qw{ emip eperl5 }] ), ],
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
                        qw{ bcftools bedtools bwa bwakit cadd chanjo chromograph
                          cnvnator cramtools cutadapt delly expansionhunter fastqc
                          gatk gatk4 genmod gcc htslib libxml2 libxslt
                          manta mip_scripts multiqc numpy peddy picard pip
                          plink python rhocall rtg-tools sambamba samtools
                          scikit-learn stranger svdb tiddit upd varg
                          variant_integrity vcf2cytosure vcfanno vep vt }
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
            isa           => ArrayRef [ enum( [qw{ picard plink2 vt }] ), ],
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
                        qw{ bcftools bedtools bwa bwakit cadd chanjo chromograph
                          cnvnator cramtools cutadapt delly expansionhunter fastqc
                          gatk gatk4 genmod gcc htslib libxml2 libxslt
                          manta mip_scripts multiqc numpy peddy picard pip
                          plink python rhocall rtg-tools sambamba samtools
                          scikit-learn stranger svdb tiddit upd varg
                          variant_integrity vcf2cytosure vcfanno vep vt }
                    ]
                ),
            ],
            required => 0,
        ),
    );

    return;
}

1;
