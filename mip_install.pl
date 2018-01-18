#!/usr/bin/env perl

use Carp;
use charnames qw{ :full :short };
use Cwd;
use Cwd qw{ abs_path };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname basename fileparse };
use File::Spec::Functions qw{ catfile catdir devnull };
use FindBin qw{ $Bin };
use Getopt::Long;
use IO::Handle;
use List::Util qw{ any };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use Time::Piece;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Array::Utils qw{ array_minus };
use Readonly;
use YAML;

## MIPs lib/
#Add MIPs internal lib
use lib catdir( $Bin, q{lib} );
use MIP::File::Format::Yaml qw{ load_yaml };
use MIP::Gnu::Coreutils qw{ gnu_rm };
use MIP::Language::Shell qw{ create_bash_file };
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Package_manager::Conda
  qw{ conda_source_activate conda_source_deactivate };
use MIP::Script::Utils qw{ help };

## Recipes
use MIP::Recipes::Install::Bedtools qw{ install_bedtools };
use MIP::Recipes::Install::Cnvnator qw{ install_cnvnator };
use MIP::Recipes::Install::Conda
  qw{ check_conda_installation setup_conda_env install_bioconda_packages };
use MIP::Recipes::Install::Mip_scripts qw{ install_mip_scripts };
use MIP::Recipes::Install::Picard qw{ install_picard };
use MIP::Recipes::Install::Pip qw{ install_pip_packages };
use MIP::Recipes::Install::Plink2 qw{ install_plink2 };
use MIP::Recipes::Install::Reference qw{ download_genome_references };
use MIP::Recipes::Install::Rhocall qw{ install_rhocall };
use MIP::Recipes::Install::Sambamba qw{ install_sambamba };
use MIP::Recipes::Install::SnpEff qw{ install_snpeff };
use MIP::Recipes::Install::Svdb qw{ install_svdb };
use MIP::Recipes::Install::Tiddit qw{ install_tiddit };
use MIP::Recipes::Install::Vcf2cytosure qw{ install_vcf2cytosure };
use MIP::Recipes::Install::Vep qw{ install_vep };
use MIP::Recipes::Install::Vt qw{ install_vt };

our $USAGE = build_usage( {} );

## Constants
Readonly my $COLON      => q{:};
Readonly my $COMMA      => q{,};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $TAB        => qq{\t};
Readonly my $UNDERSCORE => q{_};

### Set parameter default
my $config_file = catfile( $Bin, qw{ definitions install_parameters.yaml} );
my %parameter = load_yaml( { yaml_file => $config_file } );

our $VERSION = q{1.2.28};

GetOptions(
    q{see|bash_set_errexit}    => \$parameter{bash_set_errexit},
    q{snu|bash_set_nounset}    => \$parameter{bash_set_nounset},
    q{env|conda_environment:s} => \$parameter{conda_environment},
    q{cdp|conda_dir_path:s}    => \$parameter{conda_dir_path},
    q{cdu|conda_update}        => \$parameter{conda_update},
    q{bcv|bioconda=s}          => \%{ $parameter{bioconda} },
    q{pip|pip=s}               => \%{ $parameter{pip} },
    q{pyv|python_version=s}    => \$parameter{conda_packages}{python},

    # SHELL
    q{psh|prefer_shell}   => \$parameter{prefer_shell},
    q{si|shell_install:s} => \@{ $parameter{shell_install} },
    q{pic|picard:s}       => \$parameter{shell}{picard}{version},
    q{sbb|sambamba:s}     => \$parameter{shell}{sambamba}{version},
    q{bet|bedtools:s}     => \$parameter{shell}{bedtools}{version},
    q{vt|vt:s}            => \$parameter{shell}{vt}{version},
    q{plk|plink2:s}       => \$parameter{shell}{plink2}{version},
    q{snpg|snpeff_genome_versions:s{,}} => sub {
        @{ $parameter{shell}{snpeff}{snpeff_genome_versions} } =
          split /,/xms, $ARG[1];
    },
    q{v2cs|vcf2cytosure:s}          => \$parameter{shell}{v2cs}{version},
    q{vep|varianteffectpredictor:i} => \$parameter{shell}{vep}{version},
    q{vepf|vep_auto_flag:s}         => \$parameter{shell}{vep}{vep_auto_flag},
    q{vepc|vep_cache_dir:s}         => \$parameter{shell}{vep}{vep_cache_dir},
    q{vepa|vep_assemblies:s{,}}     => sub {
        @{ $parameter{shell}{vep}{vep_assemblies} } = split /,/xms, $ARG[1];
    },
    q{vepp|vep_plugins:s{,}} => sub {
        @{ $parameter{shell}{vep}{vep_plugins} } = split /,/xms, $ARG[1];
    },
    q{rhc|rhocall:s}       => \$parameter{shell}{rhocall}{version},
    q{rhcp|rhocall_path:s} => \$parameter{shell}{rhocall}{path},
    q{cnv|cnvnator:s}      => \$parameter{shell}{cnvnator}{version},
    q{cnvnr|cnvnator_root_binary:s} =>
      \$parameter{shell}{cnvnator}{cnvnator_root_binary},
    q{tid|tiddit:s} => \$parameter{shell}{tiddit}{version},
    q{svdb|svdb:s}  => \$parameter{shell}{svdb}{version},

    # Utility
    q{rd|reference_dir:s}                => \$parameter{reference_dir},
    q{rg|reference_genome_versions:s{,}} => sub {
        @{ $parameter{reference_genome_versions} } = split /,/xms, $ARG[1];
    },
    q{ppd|print_parameters_default} => sub {
        print_parameters(
            {
                parameter_href => \%parameter,
            }
        );
        exit;
    },
    q{nup|noupdate}        => \$parameter{noupdate},
    q{sp|select_program:s} => \@{ $parameter{select_program} },
    q{skip|skip_program:s} => \@{ $parameter{skip_program} },
    q{l|log:s}             => \$parameter{log_file},
    q{q|quiet}             => \$parameter{quiet},
    q{h|help}              => sub {
        say {*STDOUT} $USAGE;
        exit;
    },
    q{ver|version} => sub {
        say {*STDOUT} $NEWLINE
          . basename($PROGRAM_NAME)
          . $SPACE
          . $VERSION
          . $NEWLINE;
        exit;
    },
    q{v|verbose} => \$parameter{verbose},
  )
  or croak help(
    {
        USAGE     => $USAGE,
        exit_code => 1,
    }
  );

## Create default log name
if ( not $parameter{log_file} ) {
    ## Get local time
    my $date_time       = localtime;
    my $date_time_stamp = $date_time->datetime;

    $parameter{log_file} = catfile(
        q{mip_install} . $UNDERSCORE . $date_time_stamp . $DOT . q{log} );
}

## Initiate logger
my $log = initiate_logger(
    {
        file_path => $parameter{log_file},
        log_name  => q{mip_install},
    }
);
$log->info(
    q{Writing log messages to} . $COLON . $SPACE . $parameter{log_file} );

## Check the conda installation and get the conda path
$parameter{conda_prefix_path} = check_conda_installation(
    {
        conda_dir_path => $parameter{conda_dir_path},
        conda_env      => $parameter{conda_environment},
        quiet          => $parameter{quiet},
        verbose        => $parameter{verbose},
    }
);

if ( not $parameter{vep_cache_dir} ) {

    # Cache directory
    $parameter{shell}{vep}{vep_cache_dir} =
      catdir( $parameter{conda_prefix_path},
        q{ensembl-tools-release-} . $parameter{shell}{vep}{version}, q{cache} );
}

##########
###MAIN###
##########

# Create anonymous filehandle
my $FILEHANDLE = IO::Handle->new();

# Installation instruction file
my $file_name_path = catfile( cwd(), q{mip.sh} );

open $FILEHANDLE, q{>}, $file_name_path
  or $log->logcroak( q{Cannot write to}
      . $SPACE . q{'}
      . $file_name_path . q{'}
      . $SPACE
      . $COLON
      . $OS_ERROR
      . $NEWLINE );

## Create bash file for writing install instructions
create_bash_file(
    {
        file_name   => $file_name_path,
        FILEHANDLE  => $FILEHANDLE,
        remove_dir  => catfile( cwd(), $DOT . q{MIP} ),
        set_errexit => $parameter{bash_set_errexit},
        set_nounset => $parameter{bash_set_nounset},
        log         => $log
    }
);

$log->info( q{Writing install instructions to:} . $SPACE . $file_name_path );

## Process input parameters to get a correct combination of programs that are to be installed
%parameter = get_programs_for_installation(
    {
        log            => $log,
        parameter_href => \%parameter,
    }
);

## Seting up conda environment and installing default packages
setup_conda_env(
    {
        conda_packages_href => $parameter{conda_packages},
        conda_env           => $parameter{conda_environment},
        conda_env_path      => $parameter{conda_prefix_path},
        FILEHANDLE          => $FILEHANDLE,
        conda_update        => $parameter{conda_update},
        quiet               => $parameter{quiet},
        verbose             => $parameter{verbose},
    }
);

## Installing bioconda packages
install_bioconda_packages(
    {
        bioconda_packages_href => $parameter{bioconda},
        conda_env              => $parameter{conda_environment},
        conda_env_path         => $parameter{conda_prefix_path},
        snpeff_genome_versions_ref =>
          $parameter{shell}{snpeff}{snpeff_genome_versions},
        FILEHANDLE => $FILEHANDLE,
        quiet      => $parameter{quiet},
        verbose    => $parameter{verbose},
    }
);

## Install PIP packages
install_pip_packages(
    {
        pip_packages_href => $parameter{pip},
        quiet             => $parameter{quiet},
        conda_env         => $parameter{conda_environment},
        FILEHANDLE        => $FILEHANDLE,
    }
);

### Install shell programs
## Create dispatch table for shell installation subs
my %shell_subs = (
    picard       => \&install_picard,
    sambamba     => \&install_sambamba,
    bedtools     => \&install_bedtools,
    vt           => \&install_vt,
    snpeff       => \&install_snpeff,
    plink2       => \&install_plink2,
    rhocall      => \&install_rhocall,
    mip_scripts  => \&install_mip_scripts,
    vep          => \&install_vep,
    cnvnator     => \&install_cnvnator,
    tiddit       => \&install_tiddit,
    svdb         => \&install_svdb,
    vcf2cytosure => \&install_vcf2cytosure,
);

## Launch shell installation subroutines
SHELL_PROGRAM:
for my $shell_program ( @{ $parameter{shell_programs_to_install} } ) {
    $shell_subs{$shell_program}->(
        {
            program_parameters_href => $parameter{shell}{$shell_program},
            conda_prefix_path       => $parameter{conda_prefix_path},
            conda_environment       => $parameter{conda_environment},
            noupdate                => $parameter{noupdate},
            quiet                   => $parameter{quiet},
            verbose                 => $parameter{verbose},
            FILEHANDLE              => $FILEHANDLE,
        }
    );
}

## Download reference genome if requested
if ( $parameter{reference_dir} ) {
    download_genome_references(
        {
            reference_genome_versions_ref =>
              $parameter{reference_genome_versions},
            reference_dir_path => $parameter{reference_dir},
            conda_prefix_path  => $parameter{conda_prefix_path},
            conda_environment  => $parameter{conda_environment},
            quiet              => $parameter{quiet},
            verbose            => $parameter{verbose},
            FILEHANDLE         => $FILEHANDLE,
        }
    );
}

## Add final message to FILEHANDLE
display_final_message(
    {
        bioconda_programs_href => $parameter{bioconda},
        conda_env_name         => $parameter{conda_environment},
        conda_programs_href    => $parameter{conda_packages},
        FILEHANDLE             => $FILEHANDLE,
        pip_programs_href      => $parameter{pip},
        shell_programs_ref     => $parameter{shell_programs_to_install},
    }
);

$log->info(q{Finished writing installation instructions for MIP});

close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

#################
###SubRoutines###
#################

sub build_usage {

##Function : Build the USAGE instructions
##Returns  :
##Arguments: $script_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $script_name;

    my $tmpl = {
        script_name => {
            default     => basename($PROGRAM_NAME),
            strict_type => 1,
            store       => \$script_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $script_name [options]
    -bse/--bash_set_errexit         Set errexit in bash scripts (defaults to "0")
    -bsu/--bash_set_nounset         Set nounset in bash scripts (defaults to "0")
    -env/--conda_environment        Conda environment (Default: "")
    -cdp/--conda_dir_path           The conda directory path (Default: "HOME/miniconda")
    -cdu/--conda_update             Update conda before installing (Supply flag to enable)
    -bvc/--bioconda                 Set the module version of the programs that can be installed with bioconda (e.g. 'bwa=0.7.12')
    -pip/--pip                      Set the module version of the programs that can be installed with pip (e.g. 'genmod=3.7.2')
    -pyv/--python_version           Set the env python version (Default: "2.7")

    ## SHELL
    -psh/--prefer_shell             Shell will be used for overlapping shell and biconda installations (Supply flag to enable)
    -si/--shell_install             Install supplied programs via shell e.g. -si bedtools -si picard (Default: "")
    -pic/--picard                   Set the picard version (Default: "2.5.9"),
    -sbb/--sambamba                 Set the sambamba version (Default: "0.6.6")
    -bet/--bedtools                 Set the bedtools version (Default: "2.26.0")
    -vt/--vt                        Set the vt version (Default: "0.57")
    -plk/--plink                    Set the plink version (Default: "160224")
    -snpg/--snpeff_genome_versions  Set the snpEff genome version (Default: ["GRCh37.75", "GRCh38.82"])
    -vep/--varianteffectpredictor   Set the VEP version (Default: "90")
    -vepf/--vep_auto_flag           Set the VEP auto installer flags
    -vepc/--vep_cache_dir           Specify the cache directory to use (whole path;
                                        defaults to "[--conda_dir_path]/ensembl-tools-release-varianteffectpredictorVersion/cache")
    -vepa/--vep_assemblies          Select the assembly version (Default: ["GRCh37", "GRCh38"])
    -vepp/--vep_plugins             Supply VEP plugins (Default: "UpDownDistance, LoFtool, Lof")
    -rhc/--rhocall                  Set the rhocall version (Default: "0.4")
    -rhcp/--rhocall_path            Set the path to where to install rhocall (Defaults: "HOME/rhocall")
    -cnvn/--cnvnator                Set the cnvnator version (Default: 0.3.3)
    -cnvnr/--cnvnator_root_binary   Set the cnvnator root binary (Default: "root_v6.06.00.Linux-slc6-x86_64-gcc4.8.tar.gz")
    -tid/--tiddit                   Set the tiddit version (Default: "1.1.6")
    -svdb/--svdb                    Set the svdb version (Default: "1.0.6")
    -v2cs/--vcf2cytosure            Set the vcf2cytosure version (Default: "0.2.0")

    ## Utility
    -rd/--reference_dir             Reference(s) directory (Default: "")
    -rd/--reference_genome_versions Reference versions to download ((Default: ["GRCh37", "hg38"]))
    -ppd/--print_parameters_default Print the parameter defaults
    -nup/--noupdate                 Do not update already installed programs (Supply flag to enable)
    -skip/--skip_program            Exclude a default program from installation
    -sp/--select_program            Only install selected program(s) (e.g. -sp bedtools -sp TIDDIT)
    -l/--log                        File for writing log messages (Default: "mip_install_TIMESTAMP.log")
    -q/--quiet                      Quiet (Supply flag to enable; no output from individual program that has a quiet flag)
    -h/--help                       Display this help message
    -ver/--version                  Display version
    -v/--verbose                    Set verbosity
END_USAGE
}

sub print_parameters {

##Function : Print all parameters and the default values
##Returns  : ""
##Arguments: $parameter_href => Holds all parameters {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set default for vep cache dir
    $parameter{shell}{vep}{vep_cache_dir} = catdir( qw{ PATH TO CONDA },
        q{ensembl-tools-release-} . $parameter{shell}{vep}{version}, q{cache} );

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
    return;
}

sub get_programs_for_installation {

## Function : Procces the lists of programs that has been seleceted for installation and returns them
## Returns  : %{ $parameter_href }
## Arguments: $log            => Log
##          : $parameter_href => The entire parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;

    my $tmpl = {
        log => {
            required => 1,
            defined  => 1,
            store    => \$log,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Remove selected programs from installation
    if ( @{ $parameter_href->{skip_program} } ) {
      PROGRAM:
        foreach my $program ( @{ $parameter_href->{skip_program} } ) {
            delete $parameter_href->{shell}{$program};
            delete $parameter_href->{bioconda}{$program};
            delete $parameter_href->{pip}{$program};
        }
    }
    ## Exit if a python 3 env has ben specified for something else than Chanjo or Genmod
    _assure_python_3_compability(
        {
            log                => $log,
            py3_packages_ref   => $parameter_href->{py3_packages},
            python_version     => $parameter_href->{conda_packages}{python},
            select_program_ref => $parameter_href->{select_program}
        }
    );

    ## Chanjo requires sambamba to run and thus sambamba is added to the
    ## select program array if a Chanjo installation has been requested
    if ( any { $_ eq q{chanjo} } @{ $parameter_href->{select_program} } ) {
        push @{ $parameter_href->{select_program} }, q{sambamba};
    }

    ## vcf2cytosure requires libxml2 and libxslt to run and these are added to the
    ## select program array if a vcf2cytosure installation has been requested
    if ( any { $_ eq q{vcf2cytosure} } @{ $parameter_href->{select_program} } ) {
        push @{ $parameter_href->{select_program} }, qw{ libxml2 libxslt };
    }

    ## Remove all programs except those selected for installation
    if ( @{ $parameter_href->{select_program} } ) {
        my @programs = (
            keys %{ $parameter_href->{shell} },
            keys %{ $parameter_href->{bioconda} },
            keys %{ $parameter_href->{pip} }
        );
        my @programs_to_skip =
          array_minus( @programs, @{ $parameter_href->{select_program} } );

      PROGRAM:
        foreach my $program (@programs_to_skip) {
            delete $parameter_href->{shell}{$program};
            delete $parameter_href->{bioconda}{$program};
            delete $parameter_href->{pip}{$program};
        }
    }

    ## Some programs have conflicting dependencies and require seperate environments to function properly
    ## These are excluded from installation unless specified with the select_program flag
    my @conflicting_programs = qw{ cnvnator peddy };
  CONFLICTING_PROGRAM:
    foreach my $conflicting_program (@conflicting_programs) {
        if (
            not any { $_ eq $conflicting_program }
            @{ $parameter_href->{select_program} }
          )
        {
            delete $parameter_href->{shell}{$conflicting_program};
            delete $parameter_href->{bioconda}{$conflicting_program};
        }
    }

    ## Assure a gcc version of 4.8 in the case of a cnvnator installation
    if ( any { $_ eq q{cnvnator} } @{ $parameter_href->{select_program} } ) {
        $parameter_href->{conda_packages}{gcc} = q{4.8};
    }

    ## Exclude Chanjo, Genmod and Variant_integrity unless a python 3 env has been specified
    $parameter_href = _assure_python_2_compability(
        {
            log            => $log,
            parameter_href => $parameter_href,
        }
    );

    my @shell_programs_to_install = get_programs_for_shell_installation(
        {
            shell_programs_href        => $parameter_href->{shell},
            conda_programs_href        => $parameter_href->{bioconda},
            shell_install_programs_ref => $parameter_href->{shell_install},
            prefer_shell               => $parameter_href->{prefer_shell},
        }
    );

    ## Removing the bioconda packages that has been selected to be installed via SHELL
    delete @{ $parameter_href->{bioconda} }{@shell_programs_to_install};
    ## Special case for snpsift since it is installed together with SnpEff
    ## if shell installation of SnpEff has been requested.
    if ( any { $_ eq q{snpeff} } @shell_programs_to_install ) {
        delete $parameter_href->{bioconda}{snpsift};
    }
    $parameter_href->{shell_programs_to_install} = [@shell_programs_to_install];

    return %{$parameter_href};
}

sub _assure_python_3_compability {

## Function : Test if specified programs are to be installed in a python 3 environment
## Returns  :
## Arguments: $sub_log            => Log
##          : $py3_packages_ref   => Array with packages that requires python 3 {REF}
##          : $python_version     => The python version that are to be used for the environment
##          : $select_program_ref => Programs selected for installation by the user {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sub_log;
    my $py3_packages_ref;
    my $python_version;
    my $select_program_ref;

    my $tmpl = {
        log => {
            required => 1,
            defined  => 1,
            store    => \$sub_log,
        },
        py3_packages_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$py3_packages_ref,
        },
        python_version => {
            required => 1,
            defined  => 1,
            allow    => qr{
                         ^( 2 | 3 )    # Assert that the python major version starts with 2 or 3
                         \.            # Major version separator
                         ( \d+$        # Assert that the minor version is a digit
                         | \d+\.\d+$ ) # Case when minor and patch version has been supplied, allow only digits
                         }xms,
            store => \$python_version,
        },
        select_program_ref => {
            required => 1,
            default  => [],
            store    => \$select_program_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Array::Utils qw{ array_minus };

    ## Check if a python 3 environment has been specified and a python 2 program has been specified for installation
    if (
        $python_version =~ m/
        3\.\d+ |    # Python 3 release with minor version eg 3.6
        3\.\d+\.\d+ # Python 3 release with minor and patch e.g. 3.6.2
        /xms
        and array_minus( @{$select_program_ref}, @{$py3_packages_ref} )
      )
    {
        $sub_log->fatal(
q{A python 3 env has been specified. Please use a python 2 environment for all programs except:}
              . $NEWLINE
              . join $TAB,
            @{$py3_packages_ref}
        );
        exit 1;
    }
    return;
}

sub _assure_python_2_compability {

## Function : Exclude programs that are not compatible with python 2 and exit when a program with python 3 dependency has been selected.
## Returns  : $parameter_href
## Arguments: $log            => Log
##          : $parameter_href => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sub_log;
    my $parameter_href;

    my $tmpl = {
        log => {
            required => 1,
            defined  => 1,
            store    => \$sub_log,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Array::Utils qw{ intersect };

    if (
        $parameter_href->{conda_packages}{python} =~
        m/ 2\.\d+ |    # Python 2 release with minor version eg 2.7.14
           2\.\d+\.\d+ # Python 3 release with minor and patch e.g. 3.6.2
        /xms
      )
    {
        ## Delete python 3 packages if a python 2 env has been specified
        foreach my $py3_package ( @{ $parameter_href->{py3_packages} } ) {
            delete $parameter_href->{pip}{$py3_package};
        }

        ## Check if a python 3 package has been selected for installation in a python 2 environment
        if (
            intersect(
                @{ $parameter_href->{select_program} },
                @{ $parameter_href->{py3_packages} }
            )
          )
        {
            $sub_log->fatal(
                q{Please specify a python 3 environment for:}
                  . $NEWLINE
                  . join $TAB,
                @{ $parameter_href->{py3_packages} }
            );
            exit 1;
        }
    }

    return $parameter_href;
}

sub get_programs_for_shell_installation {

## Function  : Get the programs that are to be installed via SHELL
## Returns   : @shell_programs
## Arguments : $shell_programs_href       => Hash with shell programs {REF}
##           : $conda_programs_href       => Hash with conda progrmas {REF}
##           : $shell_install_programs_ref => Array with programs selected for shell installation {REF}
##           : $prefer_shell              => Path to conda environment

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $shell_programs_href;
    my $conda_programs_href;
    my $shell_install_programs_ref;
    my $prefer_shell;

    my $tmpl = {
        shell_programs_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$shell_programs_href,
        },
        conda_programs_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$conda_programs_href,
        },
        shell_install_programs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$shell_install_programs_ref,
        },
        prefer_shell => {
            required    => 1,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$prefer_shell
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Array::Utils qw{ intersect array_minus unique };

    my @shell_programs = keys %{$shell_programs_href};
    my @conda_programs = keys %{$conda_programs_href};

    if ($prefer_shell) {

        # Only get the selected programs otherwise leave the array unaltered
        if ( @{$shell_install_programs_ref} ) {

            # Get the intersect between the two arrays
            @shell_programs =
              intersect( @shell_programs, @{$shell_install_programs_ref} );
        }
    }
    elsif ( @{$shell_install_programs_ref} ) {

        # Assert that the selected program has shell install instructions.
        my @faulty_selects =
          array_minus( @{$shell_install_programs_ref}, @shell_programs );
        if (@faulty_selects) {
            $log->fatal(
                q{No shell installation instructions available for}
                  . $COLON
                  . join $SPACE,
                @faulty_selects
            );
            exit 1;
        }

        # Get elements in @shell_programs that are not part of the conda hash
        my @shell_only_programs =
          array_minus( @shell_programs, @conda_programs );

        # Add the selected program(s) and remove possible duplicates
        @shell_programs =
          unique( @shell_only_programs, @{$shell_install_programs_ref} );
    }
    else {
        # If no shell preferences only add programs lacking conda counterpart
        @shell_programs = array_minus( @shell_programs, @conda_programs );
    }

    return @shell_programs;
}

sub display_final_message {

## Function : Displays a final message to the user at the end of the installation process
## Returns  :
## Arguments: $bioconda_programs_href => Hash with bioconda programs {REF}
##          : $conda_env_name         => Name of conda environment
##          : $conda_programs_href    => Hash with conda programs {REF}
##          : pip_programs_href       => Hash with pip programs {REF}
##          : shell_programs_ref      => Array with shell programs {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bioconda_programs_href;
    my $conda_env_name;
    my $conda_programs_href;
    my $pip_programs_href;
    my $shell_programs_ref;

    my $tmpl = {
        bioconda_programs_href => {
            default => {},
            store   => \$bioconda_programs_href,
        },
        conda_env_name => {
            store => \$conda_env_name,
        },
        conda_programs_href => {
            default => {},
            store   => \$conda_programs_href,
        },
        pip_programs_href => {
            default => {},
            store   => \$pip_programs_href,
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE,
        },
        shell_programs_ref => {
            default => [],
            store   => \$shell_programs_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};
    use Array::Utils qw{ unique };

    ## Get the programs that mip has tried to install
    my @programs_to_install = unique(
        keys %{$bioconda_programs_href},
        keys %{$conda_programs_href},
        keys %{$pip_programs_href},
        @{$shell_programs_ref},
    );

    ## Set conda env variable to Root if no conda environment was specified
    if ( not $conda_env_name ) {
        $conda_env_name = q{Root environment};
    }

    say $FILEHANDLE
q{echo -e '\n##############################################################\n'};
    if ( any { $_ eq q{cnvnator} } @programs_to_install ) {
        say $FILEHANDLE
q{echo -e "\tMIP's installation script has attempted to install CNVnator"};
        say $FILEHANDLE q{echo -e "\tin the specified conda environment: }
          . $conda_env_name . q{\n"};
        say $FILEHANDLE
          q{echo -e "\tPlease exit the current session before continuing"};
    }
    else {
        say $FILEHANDLE q{echo -e "\tMIP's installation script has finished\n"};
        say $FILEHANDLE
          q{echo -e "\tMIP has attempted to install the following programs"};
        say $FILEHANDLE q{echo -e "\tin the specified conda environment: }
          . $conda_env_name . q{\n"};

        foreach my $program_to_install ( sort @programs_to_install ) {
            say $FILEHANDLE q{echo -e "\t"} . $program_to_install;
        }
    }
    say $FILEHANDLE
q{echo -e '\n##############################################################\n'};
    return;
}
