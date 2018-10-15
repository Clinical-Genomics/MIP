package MIP::Cli::Mip::Install::Rd_dna;

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
use MooseX::Types::Structured qw{ Dict Optional };
use Readonly;

## MIPs lib
use MIP::File::Format::Yaml qw{ load_yaml };
use MIP::Main::Install qw{ mip_install };
use MIP::Script::Utils
  qw{ nest_hash print_parameter_defaults update_program_versions};

our $VERSION = 1.00;

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

    ## Load default parameters from config file
    my $install_parameters_path = abs_path( $arg_href->config_file );
    my %parameter               = load_yaml(
        {
            yaml_file => $install_parameters_path
        }
    );

    ## Print parameters from config file and exit
    if ( $arg_href->{print_parameter_default} ) {

        ## Set default for vep cache dir
        if ( $parameter{shell}{vep} ) {
            $parameter{shell}{vep}{vep_cache_dir} = catdir( qw{ PATH TO CONDA },
                q{ensembl-tools-release-} . $parameter{shell}{vep}{version},
                q{cache} );
        }

        print_parameter_defaults(
            {
                parameter_href => \%parameter,
            }
        );
    }

    ## Merge arrays and overwrite flat values in config YAML with command line
    @parameter{ keys %{$arg_href} } = values %{$arg_href};

    ## Add all environments to installation if full installation was selected
    if ( any { $_ eq q{full} } @{ $parameter{installations} } ) {
        @{ $parameter{installations} } =
          qw{ emip ecnvnator edelly efreebayes epeddy epy3 etiddit evep };
    }

    ## Make sure that the cnvnator environment is installed last
    if ( any { $_ eq q{ecnvnator} } @{ $parameter{installations} } ) {
        @{ $parameter{installations} } =
          grep { !m/ecnvnator/xms } @{ $parameter{installations} };
        push @{ $parameter{installations} }, q{ecnvnator};
    }

    ## Nest the command line parameters and overwrite the default
    nest_hash( { cmd_href => \%parameter } );

    ## Update the program versions with the user input
    update_program_versions(
        {
            parameter_href => \%parameter,
        }
    );

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
        q{environment_name} => (
            cmd_aliases   => [qw{ envn }],
            cmd_flag      => q{environment_name},
            documentation => q{Set environment names},
            is            => q{rw},
            isa           => Dict [
                emip       => Optional [Str],
                ecnvnator  => Optional [Str],
                edelly     => Optional [Str],
                efreebayes => Optional [Str],
                epeddy     => Optional [Str],
                epy3       => Optional [Str],
                etiddit    => Optional [Str],
                evep       => Optional [Str],
            ],
            required => 0,
        ),
    );

    option(
        q{config_file} => (
            cmd_aliases => [qw{ config c }],
            documentation =>
              q{File with configuration parameters in YAML format},
            is      => q{rw},
            isa     => Str,
            default => catfile(
                dirname($Bin),
                qw{ MIP definitions install_rd_dna_parameters.yaml }
            ),
        )
    );

    option(
        q{installations} => (
            cmd_aliases   => [qw{ install }],
            cmd_flag      => q{installations},
            cmd_tags      => [q{Default: emip, epeddy, epy3, evep}],
            documentation => q{Environments to install},
            is            => q{rw},
            isa           => ArrayRef [
                enum(
                    [
                        qw{ emip ecnvnator edelly efreebayes epeddy epy3 etiddit evep full }
                    ]
                ),
            ],
            required => 0,
        ),
    );

    option(
        q{program_versions} => (
            cmd_aliases   => [qw{ pv }],
            cmd_flag      => q{program_versions},
            documentation => q{Set program versions},
            is            => q{rw},
            isa           => Dict [
                bcftools          => Optional [Str],
                bedtools          => Optional [Str],
                bwa               => Optional [Str],
                bwakit            => Optional [Str],
                chanjo            => Optional [Str],
                cmake             => Optional [Str],
                cnvnator          => Optional [Str],
                cramtools         => Optional [Str],
                cutadapt          => Optional [Str],
                delly             => Optional [Str],
                expansionhunter   => Optional [Str],
                fastqc            => Optional [Str],
                freebayes         => Optional [Str],
                gatk              => Optional [Str],
                gatk4             => Optional [Str],
                gcc               => Optional [Str],
                genmod            => Optional [Str],
                htslib            => Optional [Str],
                libxml2           => Optional [Str],
                libxslt           => Optional [Str],
                manta             => Optional [Str],
                multiqc           => Optional [Str],
                numpy             => Optional [Str],
                peddy             => Optional [Str],
                picard            => Optional [Str],
                pip               => Optional [Str],
                plink2            => Optional [Str],
                python            => Optional [Str],
                q{rtg-tools}      => Optional [Str],
                q{scikit-learn}   => Optional [Str],
                rhocall           => Optional [Str],
                sambamba          => Optional [Str],
                samtools          => Optional [Str],
                snpeff            => Optional [Str],
                snpeff            => Optional [Str],
                snpsift           => Optional [Str],
                svdb              => Optional [Str],
                tiddit            => Optional [Str],
                variant_integrity => Optional [Str],
                vcf2cytosure      => Optional [Str],
                vcfanno           => Optional [Str],
                vep               => Optional [Str],
                vt                => Optional [Str],
            ],
            required => 0,
        ),
    );

    option(
        q{shell:cnvnator:cnvnator_root_binary} => (
            cmd_aliases => [qw{ cnvnr }],
            cmd_flag    => q{cnvnator_root_binary},
            cmd_tags =>
              [q{Default: root_v6.06.00.Linux-slc6-x86_64-gcc4.8.tar.gz}],
            documentation => q{Set the cnvnator root binary},
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
                        qw{ bcftools bedtools bwa bwakit chanjo cmake cnvnator
                          cramtools cutadapt delly expansionhunter fastqc
                          freebayes gatk gatk4 genmod gcc htslib libxml2 libxslt
                          manta mip_scripts multiqc numpy peddy picard pip
                          plink python rhocall rtg-tools sambamba samtools
                          scikit-learn snpeff snpsift svdb tiddit
                          variant_integrity vcf2cytosure vcfanno vep vt }
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
            is  => q{rw},
            isa => ArrayRef [
                enum( [qw{ bedtools picard plink2 sambamba snpeff vt }] ),
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
                        qw{ bcftools bedtools bwa bwakit chanjo cmake cnvnator
                          cramtools cutadapt delly expansionhunter fastqc
                          freebayes gatk gatk4 genmod gcc htslib libxml2 libxslt
                          manta mip_scripts multiqc numpy peddy picard pip
                          plink python rhocall rtg-tools sambamba samtools
                          scikit-learn snpeff snpsift svdb tiddit
                          variant_integrity vcf2cytosure vcfanno vep vt }
                    ]
                ),
            ],
            required => 0,
        ),
    );

    option(
        q{shell:snpeff:snpeff_genome_versions} => (
            cmd_aliases   => [qw{ snpg }],
            cmd_flag      => q{snpeff_genome_versions},
            cmd_tags      => [q{Default: GRCh37.75, GRCh38.86}],
            documentation => q{Set the SnpEff genome versions},
            is            => q{rw},
            isa           => ArrayRef,
            required      => 0,
        ),
    );

    option(
        q{shell:vep:vep_auto_flag} => (
            cmd_aliases   => [qw{ vepf }],
            cmd_flag      => q{vep_auto_flag},
            cmd_tags      => [q{Default: acfp}],
            documentation => q{Set the vep auto installer flags},
            is            => q{rw},
            isa           => Str,
            required      => 0,
        ),
    );

    option(
        q{shell:vep:vep_assemblies} => (
            cmd_aliases   => [qw{ vepa }],
            cmd_flag      => q{vep_assemblies},
            cmd_tags      => [q{Default: GRCh37, hg38}],
            documentation => q{Select the assembly version},
            is            => q{rw},
            isa           => ArrayRef [ enum( [qw{ GRCh37 hg38 }] ), ],
            required      => 0,
        ),
    );

    option(
        q{shell:vep:vep_cache_dir} => (
            cmd_aliases => [qw{ vepc }],
            cmd_flag    => q{vep_cache_dir},
            cmd_tags    => [
q{Default: [path_to_conda_env]/ensembl-tools-release-[vep_version]/cache}
            ],
            documentation => q{Specify the cache directory to use},
            is            => q{rw},
            isa           => Str,
            required      => 0,
        ),
    );

    option(
        q{shell:vep:vep_plugins} => (
            cmd_aliases   => [qw{ vepp }],
            cmd_flag      => q{vep_plugins},
            cmd_tags      => [q{Default: MaxEntScan, LoFtool}],
            documentation => q{Select the vep plugins to install},
            is            => q{rw},
            isa           => ArrayRef [ enum( [qw{ MaxEntScan Loftool }] ), ],
            required      => 0,
        ),
    );

    return;
}

1;
