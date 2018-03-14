package MIP::Cli::Mip::Install::Rare_disease;

use Carp;
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
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
use MooseX::Types::Moose qw{ Str Int HashRef Num Bool ArrayRef };
use MooseX::Types::Structured qw{ Dict Optional };
use Readonly;

## MIPs lib
use MIP::File::Format::Yaml qw{ load_yaml };
use MIP::Main::Install qw{ mip_install };

our $version = 0.0.1;

extends(qw{ MIP::Cli::Mip::Install });

command_short_description(q{Generate mip.sh for installation});
command_long_description(
q{Generates an installation script (mip.sh), which is used for installation of the Mutation Identification Pipeline (MIP) for rare diseases.}
);
command_usage(q{mip <install> <rare_disease> [options]});

## Constants
Readonly my $SPACE => q{ };
Readonly my $COLON => q{:};

## Define, check and get Cli supplied parameters
_build_usage();

sub run {
    my ($arg_href) = @_;

    ## Remove Moose::App extra variable
    delete $arg_href->{extra_argv};

    ## Load default parameters from config file
    my $install_parameters_path = catfile( $Bin, $arg_href->config_file );
    my %parameter = load_yaml(
        {
            yaml_file => $install_parameters_path
        }
    );

    ## Print parameters from config file and exit
    if ( $arg_href->print_parameter_default ) {
        _print_defaults(
            {
                parameter_href => \%parameter,
            }
        );
    }

    ## Merge arrays and overwrite flat values in config YAML with command line
    @parameter{ keys %{$arg_href} } = values %{$arg_href};

    ## Nest the command line parameters and overwrite the default
    _nest_hash( { cmd_href => \%parameter } );

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
                bcftools        => Optional [Num],
                bedtools        => Optional [Num],
                bwa             => Optional [Num],
                bwakit          => Optional [Num],
                cmake           => Optional [Num],
                cramtools       => Optional [Str],
                cutadapt        => Optional [Num],
                delly           => Optional [Num],
                fastqc          => Optional [Num],
                freebayes       => Optional [Num],
                gatk            => Optional [Num],
                gcc             => Optional [Num],
                htslib          => Optional [Num],
                libxml2         => Optional [Num],
                libxslt         => Optional [Num],
                manta           => Optional [Num],
                numpy           => Optional [Num],
                peddy           => Optional [Num],
                picard          => Optional [Num],
                plink2          => Optional [Str],
                q{rtg-tools}    => Optional [Num],
                sambamba        => Optional [Num],
                samtools        => Optional [Num],
                q{scikit-learn} => Optional [Num],
                snpeff          => Optional [Num],
                snpsift         => Optional [Num],
                vcfanno         => Optional [Num],
                vt              => Optional [Num],
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
        q{conda_dir_path} => (
            cmd_aliases   => [qw{ cdp }],
            cmd_flag      => q{conda_dir_path},
            documentation => q{Path to conda_directory},
            is            => q{rw},
            isa           => Str,
            required      => 0,
        ),
    );

    option(
        q{conda_environment} => (
            cmd_aliases   => [qw{ env }],
            cmd_flag      => q{conda_environment},
            cmd_tags      => [q{Default: root/base environment}],
            documentation => q{Name of conda environmnet to use},
            is            => q{rw},
            isa           => Str,
            required      => 0,
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
        q{conda_update} => (
            cmd_aliases   => [qw{ cdu }],
            cmd_flag      => q{conda_update},
            documentation => q{Update conda},
            is            => q{rw},
            isa           => Bool,
            required      => 0,
        ),
    );

    option(
        q{disable_env_check} => (
            cmd_aliases   => [qw{ dec }],
            cmd_flag      => q{disable_env_check},
            documentation => q{Disable source environment check},
            is            => q{rw},
            isa           => Bool,
            required      => 0,
        ),
    );

    option(
        q{no_update} => (
            cmd_aliases   => [qw{ nup }],
            cmd_flag      => q{noupdate},
            documentation => q{Do not update existing shell programs},
            is            => q{rw},
            isa           => Bool,
            required      => 0,
        ),
    );

    option(
        q{pip} => (
            cmd_aliases   => [qw{ pip }],
            cmd_flag      => q{pip_programs},
            documentation => q{Set the version of programs installed via pip},
            is            => q{rw},
            isa           => Dict [
                chanjo            => Optional [Num],
                multiqc           => Optional [Num],
                genmod            => Optional [Num],
                variant_integrity => Optional [Num],
            ],
            required => 0,
        ),
    );

    option(
        q{prefer_shell} => (
            cmd_aliases => [qw{ psh }],
            cmd_flag    => q{prefer_shell},
            documentation =>
q{Shell will be used for overlapping shell and biconda installations},
            is       => q{rw},
            isa      => Bool,
            required => 0,
        ),
    );

    option(
        q{print_parameter_default} => (
            cmd_aliases   => [qw{ ppd }],
            cmd_flag      => q{print_parameter_default},
            documentation => q{print the default parameters},
            is            => q{rw},
            isa           => Bool,
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
            isa           => Num,
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
                        qw{ bcftools bedtools bwa bwakit cmake cnvnator
                          cramtools cutadapt delly fastqc freebayes gatk gcc
                          htslib libxml2 libxslt manta numpy peddy picard
                          plink rhocall rtg-tools sambamba samtools scikit-learn
                          snpeff snpsift svdb tiddit vcf2cytosure vcfanno vep vt
                          }
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
        q{shell_programs} => (
            cmd_aliases   => [qw{ shell }],
            cmd_flag      => q{shell_programs},
            documentation => q{Set shell version of programs},
            is            => q{rw},
            isa           => Dict [
                bedtools     => Optional [Num],
                cnvnator     => Optional [Num],
                plink2       => Optional [Num],
                picard       => Optional [Num],
                rhocall      => Optional [Num],
                sambamba     => Optional [Num],
                snpeff       => Optional [Str],
                svdb         => Optional [Num],
                tiddit       => Optional [Num],
                vcf2cytosure => Optional [Num],
                vep          => Optional [Num],
                vt           => Optional [Str],
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
                        qw{ bcftools bedtools bwa bwakit cmake cnvnator
                          cramtools cutadapt delly fastqc freebayes gatk gcc
                          htslib libxml2 libxslt manta numpy peddy picard
                          plink rhocall rtg-tools sambamba samtools scikit-learn
                          snpeff snpsift svdb tiddit vcf2cytosure vcfanno vep vt
                          }
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

    option(
        q{verbose} => (
            cmd_aliases   => [qw{ vb }],
            cmd_flag      => q{verbose},
            documentation => q{Turn on chatty output},
            is            => q{rw},
            isa           => Bool,
            required      => 0,
        ),
    );

    return;
}

sub _nest_hash {
## Function : If necessary, nests the command line hash to fit the structure used in mip_install
## Returns  :
## Arguments: $cmd_ref => Arguments from command line {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cmd_href;

    my $tmpl = {
        cmd_href => {
            defined  => 1,
            required => 1,
            store    => \$cmd_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Fix the shell programs into the hash
    if ( $cmd_href->{shell_programs} ) {
      PROGRAM:
        foreach my $program ( keys %{ $cmd_href->{shell_programs} } ) {

            $cmd_href->{shell}{$program}{version} =
              $cmd_href->{shell_programs}{$program};
            delete $cmd_href->{shell_programs}{$program};
        }
    }

    ## Nest the shell parameters
    my @colon_keys = grep { /:/ } keys %{$cmd_href};
  PARAMETER:
    foreach my $parameter (@colon_keys) {

        my $final_value = $cmd_href->{$parameter};
        _recursive_nesting(
            {
                array_to_shift_ref => [ ( split /:/, $parameter ) ],
                final_value => $final_value,
                hash_to_populate_href => $cmd_href,
            }
        );
        delete $cmd_href->{$parameter};
    }
    return;
}

sub _recursive_nesting {
## Function  : Recursive sub to nest values into a hash from an array
## Returns   :
## Arguments : $array_to_shift_ref    => Array of keys
##           : $final_value           => Value to be stored
##           : $hash_to_populate_href => Shift array values to this hash

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $hash_to_populate_href;
    my $array_to_shift_ref;
    my $final_value;

    my $tmpl = {
        final_value => {
            required => 1,
            store    => \$final_value,
        },
        hash_to_populate_href => {
            defined  => 1,
            required => 1,
            store    => \$hash_to_populate_href,
        },
        array_to_shift_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$array_to_shift_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Assign and remove the first value from the array
    my $value = shift @{$array_to_shift_ref};

    ## If the array is empty, give the last hash key the final value and return
    if ( scalar @{$array_to_shift_ref} == 0 ) {
        return $hash_to_populate_href->{$value} = $final_value;
    }

    ## Call same subroutine but increment the hash_ref to include the value as a key
    _recursive_nesting(
        {
            hash_to_populate_href => \%{ $hash_to_populate_href->{$value} },
            final_value           => $final_value,
            array_to_shift_ref    => $array_to_shift_ref,
        }
    );
}

sub _print_defaults {

## Function : Print all parameters and the default values
## Returns  :
## Arguments: $parameter_href => Holds all parameters {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set default for vep cache dir
    $parameter_href->{shell}{vep}{vep_cache_dir} = catdir( qw{ PATH TO CONDA },
        q{ensembl-tools-release-} . $parameter_href->{shell}{vep}{version},
        q{cache} );

    ## Looping over the parameter hash to extract keys and values
  KEY:
    foreach my $key ( keys %{$parameter_href} ) {
        ## If the first level value is not a hash or array ref
        if ( ref( $parameter_href->{$key} ) !~ / ARRAY | HASH /xms ) {
            print {*STDOUT} $key . $SPACE;
            ## Check if scalar exists and print
            if ( $parameter_href->{$key} ) {
                say {*STDOUT} $parameter_href->{$key};
            }
            ## Boolean value
            else {
                say {*STDOUT} q{0};
            }
        }
        ## If the first level value is a hash ref
        elsif ( ref( $parameter_href->{$key} ) =~ /HASH/xms ) {
            ## Loop over the next set of hash keys
          PROGRAM:
            foreach my $program ( keys %{ $parameter_href->{$key} } ) {
                ## If the value is a hash ref
                if ( ref( $parameter_href->{$key}{$program} ) =~ /HASH/xms ) {
                    ## Loop over the next set of hash keys
                  NESTED_PARAM:
                    foreach my $nested_param (
                        keys %{ $parameter_href->{$key}{$program} } )
                    {
                        ## Print the key
                        print {*STDOUT} $key
                          . $SPACE
                          . $program
                          . $SPACE
                          . $nested_param
                          . $COLON
                          . $SPACE;
                        ## If the value is an array ref
                        if (
                            ref(
                                $parameter_href->{$key}{$program}{$nested_param}
                            ) =~ /ARRAY/xms
                          )
                        {
                            ## Print array
                            say {*STDOUT} join $SPACE,
                              @{ $parameter_href->{$key}{$program}
                                  {$nested_param} };
                        }
                        else {
                            ## Otherwise print the hash value
                            say {*STDOUT}
                              $parameter_href->{$key}{$program}{$nested_param};
                        }
                    }
                }
                ## Print values
                else {
                    ## Don't print value if it is undef
                    if ( not $parameter_href->{$key}{$program} ) {
                        say {*STDOUT} $key . $SPACE . $program;
                    }
                    else {
                        ## Print hash value
                        say {*STDOUT} $key
                          . $SPACE
                          . $program
                          . $COLON
                          . $SPACE
                          . $parameter_href->{$key}{$program};
                    }
                }
            }
        }
        ## Check for ref to array and print
        elsif ( ref( $parameter_href->{$key} ) =~ /ARRAY/xms ) {
            say {*STDOUT} $key . $COLON . $SPACE . join $SPACE,
              @{ $parameter_href->{$key} };
        }
    }
    exit 0;
}

1;

