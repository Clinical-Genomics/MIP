package MIP::Cli::Mip::Install;

use 5.026;
use Carp;
use File::Spec::Functions qw{ catfile };
use FindBin qw{ $Bin };
use Params::Check qw{ check allow last_error };
use open qw{ :encoding(UTF-8) :std };
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## CPANM
use autodie qw{ :all };
use MooseX::App::Command;
use Moose::Util::TypeConstraints;
use MooseX::Types::Moose qw{ ArrayRef Bool HashRef Int Str };
use MooseX::Types::Structured qw{ Dict Optional };

## MIPs lib/
use MIP::Definition qw{ get_parameter_from_definition_files };
use MIP::Get::Parameter qw{ get_install_parameter_attribute };
use MIP::Main::Install qw{ mip_install };

our $VERSION = 1.26;

extends(qw{ MIP::Cli::Mip });

command_short_description(q{MIP install command});

command_long_description(
q{Generates an installation script (mip.sh), which is used for installation of the Mutation Identification Pipeline (MIP) for rare diseases.}
);

command_usage(q{mip <install> [options]});

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
    my %parameter = get_parameter_from_definition_files( { level => q{install}, } );

    ## If no config from cmd
    if ( not $active_parameter{config_file} ) {

        ## Use default
        $active_parameter{config_file} =
          catfile( $Bin, qw{ templates mip_install_config.yaml } );
    }

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
            cmd_tags      => [q{Default: mip}],
            documentation => q{Set environment name},
            is            => q{rw},
            isa           => Str,
        ),
    );

    option(
        q{pipelines} => (
            cmd_aliases   => [qw{ p }],
            cmd_flag      => q{pipelines},
            documentation => q{Pipelines to install. Default: rd_dna, rd_rna },
            is            => q{rw},
            isa           => ArrayRef [ enum( [qw{ rd_dna rd_rna}] ), ],
            required      => 0,
        ),
    );

    option(
        q{quiet} => (
            cmd_aliases   => [qw{ q }],
            cmd_flag      => q{quiet},
            documentation => q{Limit output from programs},
            is            => q{rw},
            isa           => Bool,
            required      => 0,
        ),
    );

    option(
        q{reference_dir} => (
            cmd_aliases   => [qw{ rd }],
            cmd_tags      => [q{Default: ""}],
            documentation => q{Reference directory},
            is            => q{rw},
            isa           => Str,
        )
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
                        qw{ arriba bedtools blobfish bootstrapann bwa bwakit bwa-mem2 cadd chanjo
                          chromograph cnvnator cyrius deepvariant delly expansionhunter fastqc gatk gatk4
                          genmod gffcompare glnexus htslib manta mip mip_scripts multiqc peddy picard plink
                          preseq python rhocall rseqc rtg-tools salmon sambamba smncopynumbercaller star
                          star-fusion stranger stringtie svdb telomerecat tiddit trim-galore ucsc upd
                          utilities varg vcf2cytosure vcfanno vep vts }
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
                        qw{ arriba bedtools blobfish bootstrapann bwa bwakit bwa-mem2 cadd chanjo
                          chromograph cnvnator cyrius deepvariant delly expansionhunter fastqc gatk gatk4
                          genmod gffcompare glnexus htslib manta mip mip_scripts multiqc peddy picard plink
                          preseq python rhocall rseqc rtg-tools salmon sambamba smncopynumbercaller star
                          star-fusion stranger stringtie svdb telomerecat tiddit trim-galore ucsc upd
                          utilities varg vcf2cytosure vcfanno vep vts }
                    ]
                ),
            ],
            required => 0,
        ),
    );

    option(
        q{test_mode} => (
            documentation => q{Run MIP int test mode, i.e. not launching any child processes},
            is            => q{rw},
            isa           => Bool,
            required      => 0,
        ),
    );

    option(
        q{vep_assemblies} => (
            cmd_tags      => [q{Default: GRCh37, GRCh38}],
            cmd_flag      => q{vep_assemblies},
            documentation => q{VEP assemblies to download},
            is            => q{rw},
            isa           => ArrayRef,
            required      => 0,
        ),
    );

    option(
        q{vep_auto_flag} => (
            cmd_flag      => q{vep_auto_flag},
            cmd_tags      => [q{Default: cf}],
            documentation => q{VEP's --AUTO flag},
            is            => q{rw},
            isa           => Str,
            required      => 0,
        ),
    );

    option(
        q{vep_cache_dir} => (
            cmd_flag => q{vep_cache_dir},
            cmd_tags =>
              [q{Default: <reference_dir>/ensembl-tools-release-<version>/cache}],
            documentation => q{VEP's cache directory},
            is            => q{rw},
            isa           => Str,
            required      => 0,
        ),
    );

    option(
        q{vep_plugins} => (
            cmd_flag      => q{vep_plugins},
            documentation => q{VEP plugins to install},
            is            => q{rw},
            isa           => ArrayRef,
            required      => 0,
        ),
    );

    option(
        q{vep_species} => (
            cmd_tags      => [q{Default: homo_sapiens_merged}],
            cmd_flag      => q{vep_species},
            documentation => q{VEP species},
            is            => q{rw},
            isa           => ArrayRef,
            required      => 0,
        ),
    );

    return;
}

1;
