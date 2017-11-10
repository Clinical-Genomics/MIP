#!/usr/bin/env perl

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use Getopt::Long;
use Cwd;
use Cwd qw{ abs_path };
use FindBin qw{ $Bin };
use IO::Handle;
use File::Basename qw{ dirname basename fileparse };
use File::Spec::Functions qw{ catfile catdir devnull };
use Time::Piece;
use List::Util qw{ any };

## CPANM
use Readonly;

## MIPs lib/
#Add MIPs internal lib
use lib catdir( $Bin, q{lib} );
use MIP::Check::Path qw{ check_dir_path_exist };
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
  qw{ setup_conda_env install_bioconda_packages };
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
use MIP::Recipes::Install::Vep qw{ install_vep };
use MIP::Recipes::Install::Vt qw{ install_vt };

our $USAGE = build_usage( {} );

## Constants
Readonly my $DOT        => q{.};
Readonly my $COMMA      => q{,};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};
Readonly my $COLON      => q{:};

### Set parameter default
my %parameter;

## Bash
$parameter{bash_set_errexit} = 0;
$parameter{bash_set_nounset} = 0;

## Conda
$parameter{conda_dir_path} = [
    catdir( $ENV{HOME}, q{miniconda} ),
    catdir( $ENV{HOME}, q{miniconda2} ),
    catdir( $ENV{HOME}, q{miniconda3} )
];
$parameter{conda_packages}{python} = q{2.7};
$parameter{conda_packages}{pip}    = undef;

## Bioconda channel
$parameter{bioconda}{bwa}          = q{0.7.15};
$parameter{bioconda}{bwakit}       = q{0.7.12};
$parameter{bioconda}{fastqc}       = q{0.11.4};
$parameter{bioconda}{cramtools}    = q{3.0.b47};
$parameter{bioconda}{samtools}     = q{1.4.1};
$parameter{bioconda}{bcftools}     = q{1.4.1};
$parameter{bioconda}{snpeff}       = q{4.3.1};
$parameter{bioconda}{snpsift}      = q{4.3.1};
$parameter{bioconda}{picard}       = q{2.9.2};
$parameter{bioconda}{htslib}       = q{1.4.1};
$parameter{bioconda}{bedtools}     = q{2.26.0};
$parameter{bioconda}{vt}           = q{2015.11.10};
$parameter{bioconda}{sambamba}     = q{0.6.6};
$parameter{bioconda}{freebayes}    = q{1.1.0};
$parameter{bioconda}{delly}        = q{0.7.7};
$parameter{bioconda}{manta}        = q{1.1.0};
$parameter{bioconda}{multiqc}      = q{0.9.1a0};
$parameter{bioconda}{peddy}        = q{0.2.9};
$parameter{bioconda}{plink2}       = q{1.90b3.35};
$parameter{bioconda}{vcfanno}      = q{0.1.0};
$parameter{bioconda}{q{rtg-tools}} = q{3.8.4};
$parameter{bioconda}{q{gatk}}      = q{3.8};

# Required for CNVnator
$parameter{bioconda}{gcc}   = q{4.8.5};
$parameter{bioconda}{cmake} = q{3.3.1};

## Bioconda pathes
# For correct softlinking in share and bin in conda env
$parameter{bioconda_patches}{bioconda_bwakit_patch}  = q{-0};
$parameter{bioconda_patches}{bioconda_snpeff_patch}  = q{r-0};
$parameter{bioconda_patches}{bioconda_snpsift_patch} = q{r-0};
$parameter{bioconda_patches}{bioconda_picard_patch}  = q{-1};
$parameter{bioconda_patches}{bioconda_manta_patch}   = q{-0};

## PIP
$parameter{pip}{genmod}            = q{3.7.2};
$parameter{pip}{variant_integrity} = q{0.0.4};
$parameter{pip}{chanjo}            = q{4.2.0};

## Programs currently installable by SHELL
$parameter{shell}{mip_scripts}{version} = q{Your current MIP version};
$parameter{shell}{picard}{version}      = q{2.3.0};
$parameter{shell}{sambamba}{version}    = q{0.6.1};
$parameter{shell}{bedtools}{version}    = q{2.25.0};
$parameter{shell}{vt}{version}          = q{gitRepo};
$parameter{shell}{plink2}{version}      = q{171013};
$parameter{shell}{snpeff}{version}      = q{v4_3s};
$parameter{shell}{vep}{version}         = q{90};
$parameter{shell}{vep}{vep_auto_flag}   = q{alcfp};
$parameter{shell}{rhocall}{version}     = q{0.4};
$parameter{shell}{rhocall}{path}        = catdir( $ENV{HOME}, q{rhocall} );
$parameter{shell}{cnvnator}{version}    = q{0.3.3};
$parameter{shell}{cnvnator}{cnvnator_root_binary} =
  q{root_v6.06.00.Linux-slc6-x86_64-gcc4.8.tar.gz};
$parameter{shell}{tiddit}{version} = q{1.1.6};
$parameter{shell}{svdb}{version}   = q{1.0.7};

## Define default array parameters
$parameter{shell}{vep}{vep_assemblies} = [qw{ GRCh37 GRCh38 }];
$parameter{shell}{vep}{vep_plugins} =
  [qw{ UpDownDistance LoFtool MaxEntScan }];

# GRCh38.86 but check current on the snpEff sourceForge
$parameter{shell}{snpeff}{snpeff_genome_versions} =
  [qw{ GRCh37.75 GRCh38.86 }];
$parameter{reference_genome_versions} = [qw{ GRCh37 hg38 }];

our $VERSION = q{1.2.19};

GetOptions(
    q{see|bash_set_errexit}    => \$parameter{bash_set_errexit},
    q{snu|bash_set_nounset}    => \$parameter{bash_set_nounset},
    q{env|conda_environment:s} => \$parameter{conda_environment},
    q{cdp|conda_dir_path:s{,}} => sub {
        @{ $parameter{conda_dir_path} } = split /,/xms, $ARG[1];
    },
    q{cdu|conda_update}     => \$parameter{conda_update},
    q{bcv|bioconda=s}       => \%{ $parameter{bioconda} },
    q{pip|pip=s}            => \%{ $parameter{pip} },
    q{pyv|python_version=s} => \$parameter{python_version},

    # SHELL
    q{psh|prefer_shell}     => \$parameter{prefer_shell},
    q{sp|select_programs:s} => \@{ $parameter{select_programs} },
    q{pic|picard:s}         => \$parameter{shell}{picard}{version},
    q{sbb|sambamba:s}       => \$parameter{shell}{sambamba}{version},
    q{bet|bedtools:s}       => \$parameter{shell}{bedtools}{version},
    q{vt|vt:s}              => \$parameter{shell}{vt}{version},
    q{plk|plink2:s}         => \$parameter{shell}{plink2}{version},
    q{snpg|snpeff_genome_versions:s{,}} => sub {
        @{ $parameter{shell}{snpeff}{snpeff_genome_versions} } = 
          split /,/xms, $ARG[1];
    },
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
    q{nup|noupdate} => \$parameter{noupdate},
    q{l|log:s}      => \$parameter{log_file},
    q{q|quiet}      => \$parameter{quiet},
    q{h|help}       => sub {
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

## Establish path to conda
my @conda_dir_paths = check_dir_path_exist(
    {
        dir_paths_ref => $parameter{conda_dir_path},
    }
);
my $conda_dir_path = $conda_dir_paths[0];
if ( not defined $conda_dir_path ) {
    $log->error(
        q{Could not find miniconda directory in}
          . $COLON
          . $SPACE
          . join $SPACE,
        @{ $parameter{conda_dir_path} }
    );
    exit 1;
}

## Update default parameter dependent on other parameters
if (   ( exists $parameter{conda_environment} )
    && ( $parameter{conda_environment} ) )
{
    $parameter{conda_prefix_path} =
      catdir( $conda_dir_path, q{envs}, $parameter{conda_environment} );
}
else {
    $parameter{conda_prefix_path} = $conda_dir_path;
}

if ( not $parameter{vep_cache_dir} ) {

    # Cache directory
    $parameter{shell}{vep}{vep_cache_dir} =
      catdir( $parameter{conda_prefix_path},
        q{ensembl-tools-release-} . $parameter{shell}{vep}{version}, q{cache} );
}

##########
###MAIN###
##########

say STDERR q{This is vep_autoflag: } . $parameter{shell}{vep}{vep_auto_flag};
say STDERR q{This is vep assemblies: } . join $SPACE,
  @{ $parameter{shell}{vep}{vep_assemblies} };
say STDERR q{This is vep plugins: } . join $SPACE,
  @{ $parameter{shell}{vep}{vep_plugins} };

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

## Check whether the user wants to do a shell installation of some or all overlapping bioconda packages
my @shell_programs_to_install = get_programs_for_shell_installation(
    {
        shell_programs_href       => $parameter{shell},
        conda_programs_href       => $parameter{bioconda},
        shell_select_programs_ref => $parameter{select_programs},
        prefer_shell              => $parameter{prefer_shell},
    }
);

## Removing the bioconda packages that has been selected to be installed via SHELL
delete @{ $parameter{bioconda} }{@shell_programs_to_install};
## Special case for snpsift since it is installed together with SnpEff
## if shell installation of SnpEff has been requested.
if ( any { $_ eq q{snpeff} } @shell_programs_to_install ) {
    delete $parameter{bioconda}{snpsift};
}

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
        bioconda_patches_href  => $parameter{bioconda_patches},
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
    picard      => \&install_picard,
    sambamba    => \&install_sambamba,
    bedtools    => \&install_bedtools,
    vt          => \&install_vt,
    snpeff      => \&install_snpeff,
    plink2      => \&install_plink2,
    rhocall     => \&install_rhocall,
    mip_scripts => \&install_mip_scripts,
    vep         => \&install_vep,
    cnvnator    => \&install_cnvnator,
    tiddit      => \&install_tiddit,
    svdb        => \&install_svdb,
);

## Launch shell installation subroutines
SHELL_PROGRAM:
for my $shell_program (@shell_programs_to_install) {
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
    -sp/--select_programs           Install supplied programs via shell e.g. -sp bedtools -sp picard (Default: "")
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
    
    ## Utility
    -rd/--reference_dir             Reference(s) directory (Default: "")
    -rd/--reference_genome_versions Reference versions to download ((Default: ["GRCh37", "hg38"]))
    -ppd/--print_parameters_default Print the parameter defaults
    -nup/--noupdate                 Do not update already installed programs (Supply flag to enable)
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

sub get_programs_for_shell_installation {

## Function  : Get the programs that are to be installed via SHELL
## Returns   : @shell_programs
## Arguments : $shell_programs_href       => Hash with shell programs {REF}
##           : $conda_programs_href       => Hash with conda progrmas {REF}
##           : $shell_select_programs_ref => Array with programs selected for shell installation {REF}
##           : $prefer_shell              => Path to conda environment

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $shell_programs_href;
    my $conda_programs_href;
    my $shell_select_programs_ref;
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
        shell_select_programs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$shell_select_programs_ref,
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
        if ( @{$shell_select_programs_ref} ) {

            # Get the intersect between the two arrays
            @shell_programs =
              intersect( @shell_programs, @{$shell_select_programs_ref} );
        }
    }
    elsif ( @{$shell_select_programs_ref} ) {

        # Assert that the selected program has shell install instructions.
        my @faulty_selects =
          array_minus( @{$shell_select_programs_ref}, @shell_programs );
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
          unique( @shell_only_programs, @{$shell_select_programs_ref} );
    }
    else {
        # If no shell preferences only add programs lacking conda counterpart
        @shell_programs = array_minus( @shell_programs, @conda_programs );
    }

    return @shell_programs;
}
