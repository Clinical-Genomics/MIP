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
use FindBin qw{ $Bin };    #Find directory of script
use IO::Handle;
use File::Basename qw{ dirname basename fileparse };
use File::Spec::Functions qw{ catfile catdir devnull };
use Readonly;
use Time::Piece;

## MIPs lib/
use lib catdir( $Bin, q{lib} );    #Add MIPs internal lib
use MIP::Language::Shell qw{ create_bash_file };
use MIP::Program::Download::Wget qw{ wget };
use MIP::Gnu::Bash qw{ gnu_cd };
use MIP::Gnu::Coreutils qw{ gnu_cp gnu_rm gnu_mv gnu_mkdir gnu_ln gnu_chmod };
use MIP::Package_manager::Conda
  qw{ conda_source_activate conda_source_deactivate };
use MIP::Script::Utils qw{ help set_default_array_parameters };
use MIP::Check::Path qw{ check_dir_path_exist };
use MIP::Package_manager::Pip qw{ pip_install };
use MIP::Log::MIP_log4perl qw{ initiate_logger };

## Recipes
use MIP::Recipes::Install::Conda
  qw{ setup_conda_env install_bioconda_packages finish_bioconda_package_install };
use MIP::Recipes::Install::Vep qw{ install_varianteffectpredictor };
use MIP::Recipes::Install::Pip qw{ install_pip_packages };

our $USAGE = build_usage( {} );

## Constants
Readonly my $DOT        => q{.};
Readonly my $COMMA      => q{,};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

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
$parameter{bioconda}{fastqc}       = q{0.11.5};
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
$parameter{mip_scripts}            = 'Your current MIP version';
$parameter{picardtools}            = '2.3.0';
$parameter{sambamba}               = '0.6.1';
$parameter{bedtools}               = '2.25.0';
$parameter{vt}                     = 'gitRepo';
$parameter{plink2}                 = '160316';
$parameter{snpeff}                 = 'v4_2';
$parameter{varianteffectpredictor} = '90';
$parameter{vep_auto_flag}          = q{alcfp};
$parameter{rhocall}                = '0.4';
$parameter{rhocall_path}           = catdir( $ENV{HOME}, 'rhocall' );
$parameter{cnvnator}               = '0.3.3';
$parameter{cnvnator_root_binary} =
  'root_v6.06.00.Linux-slc6-x86_64-gcc4.8.tar.gz';
$parameter{tiddit} = q{1.1.6};
$parameter{svdb}   = q{1.0.7};

## Define default parameters
my %array_parameter;

$array_parameter{vep_assemblies}{default} = [qw(GRCh37 GRCh38)];
$array_parameter{vep_plugins}{default} =
  [qw{ UpDownDistance LoFtool MaxEntScan }];

# GRCh38.86 but check current on the snpEff sourceForge
$array_parameter{snpeff_genome_versions}{default} =
  [qw(GRCh37.75 GRCh38.86)];
$array_parameter{reference_genome_versions}{default} = [qw(GRCh37 hg38)];

my $VERSION = q{1.2.12};

###User Options
GetOptions(
    q{see|bash_set_errexit}    => \$parameter{bash_set_errexit},
    q{snu|bash_set_nounset}    => \$parameter{bash_set_nounset},
    q{env|conda_environment:s} => \$parameter{conda_environment},
    q{cdp|conda_dir_path:s}    => \@{ $parameter{conda_dir_path} },
    q{cdu|conda_update}        => \$parameter{conda_update},
    q{bcv|bioconda=s}          => \%{ $parameter{bioconda} },
    q{pip|pip=s}               => \%{ $parameter{pip} },
    q{pyv|python_version=s}    => \$parameter{python_version},
    q{pic|picardtools:s}       => \$parameter{picardtools},
    q{sbb|sambamba:s}          => \$parameter{sambamba},
    q{bet|bedtools:s}          => \$parameter{bedtools},
    q{vt|vt:s}                 => \$parameter{vt},
    q{plk|plink2:s}            => \$parameter{plink2},
    q{snpg|snpeff_genome_versions:s} =>
      \@{ $parameter{snpeff_genome_versions} },
    q{vep|varianteffectpredictor:s} => \$parameter{varianteffectpredictor},
    q{vepai|vep_auto_flag:s}        => \$parameter{vep_auto_flag},
    q{vepc|vep_cache_dir:s}         => \$parameter{vep_cache_dir},
    q{vepa|vep_assemblies:s}        => \@{ $parameter{vep_assemblies} },
    q{vepp|vep_plugins:s}           => \@{ $parameter{vep_plugins} },
    q{rhc|rhocall:s}                => \$parameter{rhocall},
    q{rhcp|rhocall_path:s}          => \$parameter{rhocall_path},
    q{cnv|cnvnator:s}               => \$parameter{cnvnator},
    q{cnvnr|cnvnator_root_binary:s} => \$parameter{cnvnator_root_binary},
    q{tid|tiddit:s}                 => \$parameter{tiddit},
    q{svdb|svdb:s}                  => \$parameter{svdb},
    q{psh|prefer_shell}             => \$parameter{prefer_shell},

    # Display parameter defaults
    q{ppd|print_parameters_default} => sub {
        print_parameters(
            {
                parameter_href       => \%parameter,
                array_parameter_href => \%array_parameter,
            }
        );
        exit;
    },
    q{nup|noupdate}         => \$parameter{noupdate},
    q{sp|select_programs:s} => \@{ $parameter{select_programs} },
    q{rd|reference_dir:s}   => \$parameter{reference_dir},
    q{l|log:s}              => \$parameter{log_file},
    q{rg|reference_genome_versions:s} =>
      \@{ $parameter{reference_genome_versions} },
    q{q|quiet} => \$parameter{quiet},

    #Display help text
    q{h|help} => sub {
        print STDOUT $USAGE, "\n";
        exit;
    },

    #Display version number
    q{ver|version} => sub {
        print STDOUT "\n" . basename($PROGRAM_NAME) . q{ } . $VERSION, "\n\n";
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

$log->info( q{Writing log messages to:} . $SPACE . $parameter{log_file} );

## Establish path to conda
my @conda_dir_paths = check_dir_path_exist(
    {
        dir_paths_ref => $parameter{conda_dir_path},
    }
);
my $conda_dir_path = $conda_dir_paths[0];
if ( not defined $conda_dir_path ) {
    $log->error(
        q{Could not find miniconda directory in:} . $SPACE . join $SPACE,
        @{ $parameter{conda_dir_path} } );
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

if ( !$parameter{vep_cache_dir} ) {

    # Cache directory
    $parameter{vep_cache_dir} = catdir( $parameter{conda_prefix_path},
        'ensembl-tools-release-' . $parameter{varianteffectpredictor},
        'cache' );
}

## Set default for array parameters
set_default_array_parameters(
    {
        parameter_href       => \%parameter,
        array_parameter_href => \%array_parameter,
    }
);

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
      . $SPACE . q{:}
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
        FILEHANDLE             => $FILEHANDLE,
    }
);

## Custom solutions for BWA, SnpEff and Manta
## Copying files, downloading necessary databases and make files executable
finish_bioconda_package_install(
    {
        bioconda_packages_href     => $parameter{bioconda},
        bioconda_patches_href      => $parameter{bioconda_patches},
        conda_env                  => $parameter{conda_environment},
        conda_env_path             => $parameter{conda_prefix_path},
        FILEHANDLE                 => $FILEHANDLE,
        snpeff_genome_versions_ref => $parameter{snpeff_genome_versions},
        verbose                    => $parameter{verbose},
        quiet                      => $parameter{quiet},
    }
);

install_pip_packages(
    {
        pip_packages_href => $parameter{pip},
        quiet             => $parameter{quiet},
        conda_env         => $parameter{conda_environment},
        FILEHANDLE        => $FILEHANDLE,
    }
);

if ( $parameter{prefer_shell} ) {

    if ( @{ $parameter{select_programs} } ) {

        if ( ( grep { $_ eq 'picardtools' } @{ $parameter{select_programs} } ) )
        {    #If element is part of array

            picardtools(
                {
                    parameter_href => \%parameter,
                    FILEHANDLE     => $FILEHANDLE,
                }
            );
        }
        if ( ( grep { $_ eq 'sambamba' } @{ $parameter{select_programs} } ) )
        {    #If element is part of array

            sambamba(
                {
                    parameter_href => \%parameter,
                    FILEHANDLE     => $FILEHANDLE,
                }
            );
        }
        if ( ( grep { $_ eq 'bedtools' } @{ $parameter{select_programs} } ) )
        {    #If element is part of array

            bedtools(
                {
                    parameter_href => \%parameter,
                    FILEHANDLE     => $FILEHANDLE,
                }
            );
        }
        if ( ( grep { $_ eq 'vt' } @{ $parameter{select_programs} } ) )
        {    #If element is part of array

            vt(
                {
                    parameter_href => \%parameter,
                    FILEHANDLE     => $FILEHANDLE,
                }
            );
        }
        if ( ( grep { $_ eq 'snpeff' } @{ $parameter{select_programs} } ) )
        {    #If element is part of array

            snpeff(
                {
                    parameter_href => \%parameter,
                    FILEHANDLE     => $FILEHANDLE,
                }
            );
        }
        if ( ( grep { $_ eq 'plink2' } @{ $parameter{select_programs} } ) )
        {    #If element is part of array

            plink2(
                {
                    parameter_href => \%parameter,
                    FILEHANDLE     => $FILEHANDLE,
                }
            );
        }
        if ( ( grep { $_ eq 'rhocall' } @{ $parameter{select_programs} } ) )
        {    #If element is part of array

            rhocall(
                {
                    parameter_href => \%parameter,
                    FILEHANDLE     => $FILEHANDLE,
                }
            );
        }
    }
    else {

        picardtools(
            {
                parameter_href => \%parameter,
                FILEHANDLE     => $FILEHANDLE,
            }
        );

        sambamba(
            {
                parameter_href => \%parameter,
                FILEHANDLE     => $FILEHANDLE,
            }
        );

        bedtools(
            {
                parameter_href => \%parameter,
                FILEHANDLE     => $FILEHANDLE,
            }
        );

        vt(
            {
                parameter_href => \%parameter,
                FILEHANDLE     => $FILEHANDLE,
            }
        );

        snpeff(
            {
                parameter_href => \%parameter,
                FILEHANDLE     => $FILEHANDLE,
            }
        );

        plink2(
            {
                parameter_href => \%parameter,
                FILEHANDLE     => $FILEHANDLE,
            }
        );

        rhocall(
            {
                parameter_href => \%parameter,
                FILEHANDLE     => $FILEHANDLE,
            }
        );
    }
}

if ( @{ $parameter{select_programs} } ) {

    if ( ( grep { $_ eq 'mip_scripts' } @{ $parameter{select_programs} } ) )
    {    #If element is part of array

        mip_scripts(
            {
                parameter_href => \%parameter,
                FILEHANDLE     => $FILEHANDLE,
            }
        );
    }
    if (
        (
            grep { $_ eq q{varianteffectpredictor} }
            @{ $parameter{select_programs} }
        )
      )
    {

        install_varianteffectpredictor(
            {
                plugins_ref       => \@{ $parameter{vep_plugins} },
                assemblies_ref    => \@{ $parameter{vep_assemblies} },
                conda_prefix_path => $parameter{conda_prefix_path},
                conda_environment => $parameter{conda_environment},
                noupdate          => $parameter{noupdate},
                vep_version       => $parameter{varianteffectpredictor},
                auto              => $parameter{vep_auto_flag},
                cache_directory   => $parameter{vep_cache_dir},
                quiet             => $parameter{quiet},
                verbose           => $parameter{verbose},
                FILEHANDLE        => $FILEHANDLE,
            }
        );
    }
    if ( ( grep { $_ eq 'cnvnator' } @{ $parameter{select_programs} } ) )
    {    #If element is part of array

        cnvnator(
            {
                parameter_href => \%parameter,
                FILEHANDLE     => $FILEHANDLE,
            }
        );
    }
    if ( ( grep { $_ eq 'tiddit' } @{ $parameter{select_programs} } ) )
    {    #If element is part of array

        tiddit(
            {
                parameter_href => \%parameter,
                FILEHANDLE     => $FILEHANDLE,
            }
        );
    }
    if ( ( grep { $_ eq 'svdb' } @{ $parameter{select_programs} } ) )
    {    #If element is part of array

        svdb(
            {
                parameter_href => \%parameter,
                FILEHANDLE     => $FILEHANDLE,
            }
        );
    }
}
else {

    mip_scripts(
        {
            parameter_href => \%parameter,
            FILEHANDLE     => $FILEHANDLE,
        }
    );

    install_varianteffectpredictor(
        {
            plugins_ref       => \@{ $parameter{vep_plugins} },
            assemblies_ref    => \@{ $parameter{vep_assemblies} },
            conda_prefix_path => $parameter{conda_prefix_path},
            conda_environment => $parameter{conda_environment},
            noupdate          => $parameter{noupdate},
            vep_version       => $parameter{varianteffectpredictor},
            auto              => $parameter{vep_auto_flag},
            cache_directory   => $parameter{vep_cache_dir},
            quiet             => $parameter{quiet},
            verbose           => $parameter{verbose},
            FILEHANDLE        => $FILEHANDLE,
        }
    );

    rhocall(
        {
            parameter_href => \%parameter,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    cnvnator(
        {
            parameter_href => \%parameter,
            FILEHANDLE     => $FILEHANDLE,
        }
    );

    tiddit(
        {
            parameter_href => \%parameter,
            FILEHANDLE     => $FILEHANDLE,
        }
    );

    svdb(
        {
            parameter_href => \%parameter,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
}

if ( exists( $parameter{reference_dir} ) && ( $parameter{reference_dir} ) ) {

    references(
        {
            parameter_href => \%parameter,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
}

close($FILEHANDLE);

#close($LOGFILEHANDLE);

#################
###SubRoutines###
#################

sub build_usage {

##build_usage

##Function : Build the USAGE instructions
##Returns  : ""
##Arguments: $script_name
##         : $script_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $script_name;

    my $tmpl = {
        script_name => {
            default     => basename($0),
            strict_type => 1,
            store       => \$script_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $script_name [options]
    -env/--conda_environment        Conda environment (Default: "")
    -cdp/--conda_dir_path           The conda directory path (Default: "HOME/miniconda")
    -cdu/--conda_update             Update conda before installing (Supply flag to enable)
    -bvc/--bioconda                 Set the module version of the programs that can be installed with bioconda (e.g. 'bwa=0.7.12')
    -pip/--pip                      Set the module version of the programs that can be installed with pip (e.g. 'genmod=3.7.2')
    -pyv/--python_version           Set the env python version (Default: "2.7")

    ## SHELL
    -pic/--picardtools              Set the picardtools version (Default: "2.5.9"),
    -sbb/--sambamba                 Set sthe 
    sambamba version (Default: "0.6.6")
    -bet/--bedtools                 Set the bedtools version (Default: "2.26.0")
    -vt/--vt                        Set the vt version (Default: "0.57")
    -plk/--plink                    Set the plink version (Default: "160224")
    -snpg/--snpeff_genome_versions  Set the snpEff genome version (Default: ["GRCh37.75", "GRCh38.82"])
    -vep/--varianteffectpredictor   Set the VEP version (Default: "90")
    -vepa/--vep_auto_flag           Set the VEP auto installer flags
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
    -psh/--prefer_shell             Shell will be used for overlapping shell and biconda installations (Supply flag to enable)
    -ppd/--print_parameters_default Print the parameter defaults
    -nup/--noupdate                 Do not update already installed programs (Supply flag to enable)
    -sp/--select_programs           Install supplied programs e.g. -sp perl -sp bedtools (Default: "")
    -rd/--reference_dir             Reference(s) directory (Default: "")
    -rd/--reference_genome_versions Reference versions to download ((Default: ["GRCh37", "hg38"]))
    -l/--log                        File for writing log messages (Default: "mip_install_TIMESTAMP.log")
    -q/--quiet                      Quiet (Supply flag to enable; no output from individual program that has a quiet flag)
    -h/--help                       Display this help message
    -ver/--version                  Display version
    -v/--verbose                    Set verbosity
END_USAGE
}

sub print_parameters {

##print_parameters

##Function : Print all parameters and the default values
##Returns  : ""
##Arguments: $parameter_href, $array_parameter_href
##         : $parameter_href => Holds all parameters {REF}
##         : $array_parameter_href => Hold the array parameter defaults as {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $array_parameter_href;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        array_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$array_parameter_href
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    use MIP::Script::Utils qw{ set_default_array_parameters };

    ## Set default for array parameters
    set_default_array_parameters(
        {
            parameter_href       => $parameter_href,
            array_parameter_href => \%array_parameter,
        }
    );

    foreach my $key ( keys %{$parameter_href} ) {

        if ( ref( $parameter_href->{$key} ) !~ /ARRAY|HASH/ ) {

            print STDOUT $key . q{ };
            if ( $parameter_href->{$key} ) {

                print $parameter_href->{$key}, "\n";
            }
            else {    ##Boolean value

                print '0', "\n";
            }
        }
        elsif ( ref( $parameter_href->{$key} ) =~ /HASH/ ) {

            foreach my $program ( keys %{ $parameter_href->{$key} } ) {

                print STDOUT $key . q{ } . $program . q{: }
                  . $parameter_href->{$key}{$program}, "\n";
            }
        }
        elsif ( ref( $parameter_href->{$key} ) =~ /ARRAY/ ) {

            print STDOUT $key . q{: }
              . join( " ", @{ $parameter_href->{$key} } ), "\n";
        }
    }
    return;
}

sub picardtools {

##picardtools

##Function : Install picardtools
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (
        check_conda_bin_file_exists(
            {
                parameter_href => $parameter_href,
                program_name   => 'picard.jar',
            }
        )
      )
    {    # Assumes that picard.jar is there as well then

        return;
    }

    ## Install picard
    print $FILEHANDLE '### Install Picard', "\n";

    ## Create the temporary install directory
    create_install_dir( { FILEHANDLE => $FILEHANDLE, } );

    ## Download
    print $FILEHANDLE '## Download Picard', "\n";
    wget(
        {
            url => 'https://github.com/broadinstitute/picard/releases/download/'
              . $parameter_href->{picardtools}
              . '/picard-tools-'
              . $parameter_href->{picardtools} . '.zip',
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $parameter_href->{quiet},
            verbose      => $parameter_href->{verbose},
            outfile_path => 'picard-tools-'
              . $parameter_href->{picardtools} . '.zip',
        }
    );
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE '## Extract', "\n";
    print $FILEHANDLE 'unzip picard-tools-'
      . $parameter_href->{picardtools} . '.zip';
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    if (
        -d catdir(
            $parameter_href->{conda_prefix_path}, 'share',
            'picard-tools-' . $parameter_href->{picardtools}
        )
      )
    {

        gnu_rm(
            {
                infile_path => catdir(
                    $parameter_href->{conda_prefix_path}, 'share',
                    'picard-tools-' . $parameter_href->{picardtools}
                ),
                force      => 1,
                recursive  => 1,
                FILEHANDLE => $FILEHANDLE,
            }
        );
        print $FILEHANDLE "\n\n";
    }

    print $FILEHANDLE '## Make available from conda environment', "\n";
    gnu_mv(
        {
            infile_path => 'picard-tools-' . $parameter_href->{picardtools},
            outfile_path =>
              catdir( $parameter_href->{conda_prefix_path}, 'share' ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Specifying target and link paths
    my $target_path = catfile(
        $parameter_href->{conda_prefix_path},              q{share},
        q{picard-tools-} . $parameter_href->{picardtools}, q{picard.jar}
    );
    my $link_path =
      catfile( $parameter_href->{conda_prefix_path}, q{picard.jar} );
    gnu_ln(
        {
            FILEHANDLE  => $FILEHANDLE,
            target_path => $target_path,
            link_path   => $link_path,
            symbolic    => 1,
            force       => 1,
        }
    );
    print $FILEHANDLE $NEWLINE;

    ## Remove the temporary install directory
    remove_install_dir(
        {
            FILEHANDLE => $FILEHANDLE,
            pwd        => $pwd,
        }
    );
    return;
}

sub sambamba {

##sambamba

##Function : Install sambamba
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (
        check_conda_bin_file_exists(
            {
                parameter_href  => $parameter_href,
                program_name    => 'sambamba',
                program_version => $parameter_href->{sambamba},
            }
        )
      )
    {
        return;
    }

    ## Install sambamba
    print $FILEHANDLE '### Install sambamba', "\n";

    ## Create the temporary install directory
    create_install_dir( { FILEHANDLE => $FILEHANDLE, } );

    ## Download
    print $FILEHANDLE '## Download sambamba release', "\n";
    wget(
        {
            url => 'https://github.com/lomereiter/sambamba/releases/download/v'
              . q?$parameter_href->{sambamba}?
              . '/sambamba_v'
              . q?$parameter_href->{sambamba}?
              . '_linux.tar.bz2',
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $parameter_href->{quiet},
            verbose      => $parameter_href->{verbose},
            outfile_path => 'sambamba_v'
              . $parameter_href->{sambamba}
              . '_linux.tar.bz2',
        }
    );
    print $FILEHANDLE "\n\n";

    ## Decompress
    print $FILEHANDLE '## Decompress sambamba file', "\n";
    print $FILEHANDLE 'bzip2 ';
    print $FILEHANDLE '-f ';    #Force
    print $FILEHANDLE '-d ';    #Decompress
    print $FILEHANDLE 'sambamba_v'
      . $parameter_href->{sambamba}
      . '_linux.tar.bz2';
    print $FILEHANDLE "\n\n";

    ## Extract files
    print $FILEHANDLE '## Extract files', "\n";
    print $FILEHANDLE 'tar xvf sambamba_v'
      . $parameter_href->{sambamba}
      . '_linux.tar';
    print $FILEHANDLE "\n\n";

    ## Make executable
    print $FILEHANDLE '## Make executable', "\n";
    print $FILEHANDLE 'chmod 755 ';
    print $FILEHANDLE 'sambamba_v' . $parameter_href->{sambamba};
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE '## Make available from conda environment', "\n";
    gnu_mv(
        {
            infile_path => 'sambamba_v' . $parameter_href->{sambamba},
            outfile_path =>
              catdir( $parameter_href->{conda_prefix_path}, 'bin' ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Specifying target and link paths
    my $target_path = q{sambamba_v} . $parameter_href->{bioconda}{sambamba};
    my $link_path =
      catfile( $parameter_href->{conda_prefix_path}, q{sambamba} );
    gnu_ln(
        {
            FILEHANDLE  => $FILEHANDLE,
            target_path => $target_path,
            link_path   => $link_path,
            symbolic    => 1,
            force       => 1,
        }
    );
    print $FILEHANDLE $NEWLINE;

    ## Remove the temporary install directory
    remove_install_dir(
        {
            FILEHANDLE => $FILEHANDLE,
            pwd        => $pwd,
        }
    );
    return;
}

sub bedtools {

##bedtools

##Function : Install bedtools
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (
        check_conda_bin_file_exists(
            {
                parameter_href => $parameter_href,
                program_name   => 'bedtools',
            }
        )
      )
    {

        return;
    }

    my $bedtools_main_version = substr( $parameter_href->{bedtools}, 0, 1 );

    ## Install bedtools
    print $FILEHANDLE '### Install bedtools', "\n";

    ## Create the temporary install directory
    create_install_dir( { FILEHANDLE => $FILEHANDLE, } );

    ## Download
    print $FILEHANDLE '## Download bedtools', "\n";
    wget(
        {
            url => 'https://github.com/arq5x/bedtools'
              . $bedtools_main_version
              . '/releases/download/v'
              . $parameter_href->{bedtools}
              . '/bedtools-'
              . $parameter_href->{bedtools}
              . '.tar.gz',
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $parameter_href->{quiet},
            verbose      => $parameter_href->{verbose},
            outfile_path => 'bedtools-'
              . $parameter_href->{bedtools}
              . '.tar.gz',
        }
    );
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE '## Extract', "\n";
    print $FILEHANDLE 'tar xvf bedtools-'
      . $parameter_href->{bedtools}
      . '.tar.gz';
    print $FILEHANDLE "\n\n";

    ## Move to bedtools directory
    print $FILEHANDLE '## Move to bedtools directory', "\n";
    gnu_cd(
        {
            directory_path => 'bedtools' . $bedtools_main_version,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE 'make';
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE '## Make available from conda environment', "\n";
    gnu_mv(
        {
            infile_path => catfile(qw(. bin * )),
            outfile_path =>
              catdir( $parameter_href->{conda_prefix_path}, 'bin' ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Remove the temporary install directory
    remove_install_dir(
        {
            FILEHANDLE => $FILEHANDLE,
            pwd        => $pwd,
        }
    );
    return;
}

sub vt {

##vt

##Function : Install vt
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (
        check_conda_bin_file_exists(
            {
                parameter_href => $parameter_href,
                program_name   => 'vt',
            }
        )
      )
    {

        return;
    }

    ## Install vt
    print $FILEHANDLE '### Install VT', "\n";

    ## Create the temporary install directory
    create_install_dir( { FILEHANDLE => $FILEHANDLE, } );

    ## Download
    print $FILEHANDLE '## Download VT', "\n";

    print $FILEHANDLE 'git clone https://github.com/atks/vt.git ';
    print $FILEHANDLE "\n\n";

    ## Move to vt directory
    print $FILEHANDLE '## Move to vt directory', "\n";
    gnu_cd(
        {
            directory_path => 'vt',
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Configure
    print $FILEHANDLE '## Configure', "\n";
    print $FILEHANDLE 'make',         "\n";
    print $FILEHANDLE 'make test';
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE '## Make available from conda environment', "\n";
    gnu_mv(
        {
            infile_path => 'vt',
            outfile_path =>
              catdir( $parameter_href->{conda_prefix_path}, 'bin' ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Remove the temporary install directory
    remove_install_dir(
        {
            FILEHANDLE => $FILEHANDLE,
            pwd        => $pwd,
        }
    );
    return;
}

sub plink2 {

##plink2

##Function : Install plink2
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (
        check_conda_bin_file_exists(
            {
                parameter_href => $parameter_href,
                program_name   => 'plink',
            }
        )
      )
    {

        return;
    }

    ## Install Plink
    print $FILEHANDLE '### Install Plink', "\n";

    ## Create the temporary install directory
    create_install_dir( { FILEHANDLE => $FILEHANDLE, } );

    ## Download
    print $FILEHANDLE '## Download Plink', "\n";
    wget(
        {
            url => 'https://www.cog-genomics.org/static/bin/plink'
              . $parameter_href->{plink2}
              . '/plink_linux_x86_64.zip',
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $parameter_href->{quiet},
            verbose      => $parameter_href->{verbose},
            outfile_path => 'plink-'
              . $parameter_href->{plink2}
              . '-x86_64.zip',
        }
    );
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE '## Extract', "\n";
    print $FILEHANDLE 'unzip plink-'
      . $parameter_href->{plink2}
      . '-x86_64.zip';
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE '## Make available from conda environment', "\n";
    gnu_mv(
        {
            infile_path => 'plink',
            outfile_path =>
              catdir( $parameter_href->{conda_prefix_path}, qw(bin plink2) ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Remove the temporary install directory
    remove_install_dir(
        {
            FILEHANDLE => $FILEHANDLE,
            pwd        => $pwd,
        }
    );
    return;
}

sub snpeff {

##snpeff

##Function : Install snpeff
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (
        check_conda_bin_file_exists(
            {
                parameter_href => $parameter_href,
                program_name   => 'snpEff.jar',
            }
        )
      )
    {    # Assumes that SnpSift.jar is there as well then

        return;
    }

    ## Install snpeff
    print $FILEHANDLE '### Install snpeff', "\n";

    ## Create the temporary install directory
    create_install_dir( { FILEHANDLE => $FILEHANDLE, } );

    ## Download
    print $FILEHANDLE '## Download snpeff', "\n";
    wget(
        {
            url => 'http://sourceforge.net/projects/snpeff/files/snpEff_'
              . $parameter_href->{snpeff}
              . '_core.zip/download',
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $parameter_href->{quiet},
            verbose      => $parameter_href->{verbose},
            outfile_path => 'snpEff_' . $parameter_href->{snpeff} . '_core.zip',
        }
    );
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE '## Extract', "\n";
    print $FILEHANDLE 'unzip snpEff_' . $parameter_href->{snpeff} . '_core.zip';
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    if (
        -d catdir(
            $parameter_href->{conda_prefix_path}, 'share',
            'snpEff',                             $parameter_href->{snpeff}
        )
      )
    {

        gnu_rm(
            {
                infile_path => catdir(
                    $parameter_href->{conda_prefix_path},
                    'share', 'snpEff', $parameter_href->{snpeff}
                ),
                force      => 1,
                recursive  => 1,
                FILEHANDLE => $FILEHANDLE,
            }
        );
        print $FILEHANDLE "\n\n";
    }

    print $FILEHANDLE '## Make available from conda environment', "\n";
    gnu_mkdir(
        {
            indirectory_path => catdir(
                $parameter_href->{conda_prefix_path}, 'share',
                'snpEff.' . $parameter_href->{snpeff}
            ),
            parents    => 1,
            FILEHANDLE => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    gnu_mv(
        {
            infile_path  => catfile(qw(snpEff *.jar)),
            outfile_path => catdir(
                $parameter_href->{conda_prefix_path}, 'share',
                'snpEff.' . $parameter_href->{snpeff}
            ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    gnu_mv(
        {
            infile_path  => catfile(qw(snpEff snpEff.config)),
            outfile_path => catdir(
                $parameter_href->{conda_prefix_path}, 'share',
                'snpEff.' . $parameter_href->{snpeff}
            ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Define binaries
    my @snpeff_binaries = qw(snpEff.jar SnpSift.jar snpEff.config);

    foreach my $binary (@snpeff_binaries) {
        ## Specifying target and link path
        my $target_path = catfile(
            $parameter_href->{conda_prefix_path},   q{share},
            q{snpEff.} . $parameter_href->{snpeff}, $binary
        );
        my $link_path =
          catfile( $parameter_href->{conda_prefix_path}, $binary );
        gnu_ln(
            {
                FILEHANDLE  => $FILEHANDLE,
                target_path => $target_path,
                link_path   => $link_path,
                symbolic    => 1,
                gorce       => 1,
            }
        );
        print $FILEHANDLE $NEWLINE;
    }

    foreach
      my $genome_version ( @{ $parameter_href->{snpeff_genome_versions} } )
    {
        ## Check and if required add the vertebrate mitochondrial codon table to snpeff config
        check_mt_codon_table(
            {
                parameter_href => $parameter_href,
                FILEHANDLE     => $FILEHANDLE,
                share_dir      => catdir(
                    $parameter_href->{conda_prefix_path}, 'share',
                    'snpEff.' . $parameter_href->{snpeff}
                ),
                config_file        => 'snpEff.config',
                genome_version_ref => \$genome_version,
            }
        );

        unless (
            -d catdir(
                $parameter_href->{conda_prefix_path},  'share',
                'snpEff.' . $parameter_href->{snpeff}, 'data',
                $genome_version
            )
          )
        {

            ## Write instructions to download snpeff database.
            ## This is done by install script to avoid race conditin when doing first analysis run in MIP
            snpeff_download(
                {
                    parameter_href     => $parameter_href,
                    FILEHANDLE         => $FILEHANDLE,
                    genome_version_ref => \$genome_version,
                }
            );
        }
    }

    ## Remove the temporary install directory
    remove_install_dir(
        {
            FILEHANDLE => $FILEHANDLE,
            pwd        => $pwd,
        }
    );
    return;
}

sub cnvnator {

##cnvnator

##Function : Install cnvnator
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (
        check_conda_bin_file_exists(
            {
                parameter_href => $parameter_href,
                program_name   => 'cnvnator',
            }
        )
      )
    {

        return;
    }

    my $miniconda_bin_dir =
      catdir( $parameter_href->{conda_prefix_path}, 'root' );

    if ( -d $miniconda_bin_dir ) {

        $log->info( q{Found root in miniconda directory:}
              . $SPACE
              . $miniconda_bin_dir );

        if ( $parameter_href->{noupdate} ) {

            $log->info(q{Skipping writting installation process for Root});
            return;
        }
        else {

            ## Removing Root
            print $FILEHANDLE '### Removing Root', "\n";
            gnu_rm(
                {
                    infile_path => $miniconda_bin_dir,
                    force       => 1,
                    recursive   => 1,
                    FILEHANDLE  => $FILEHANDLE,
                }
            );
            print $FILEHANDLE "\n\n";
        }
    }
    else {

        $log->info(q{Writting install instructions for Root});
    }

    ## Install Root
    print $FILEHANDLE '### Install cnvnator/Root', "\n";

    ## Move to miniconda environment
    gnu_cd(
        {
            directory_path => catdir( $parameter_href->{conda_prefix_path} ),
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Download
    print $FILEHANDLE '## Download Root', "\n";
    wget(
        {
            url => 'https://root.cern.ch/download/'
              . $parameter_href->{cnvnator_root_binary}
            , #root_v5.34.34.Linux-slc6-x86_64-gcc4.4.tar.gz", #Currently hardcoded
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $parameter_href->{quiet},
            verbose      => $parameter_href->{verbose},
            outfile_path => $parameter_href->{cnvnator_root_binary},
        }
    );
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE '## Extract', "\n";
    print $FILEHANDLE 'tar xvf '
      . $parameter_href->{cnvnator_root_binary} . q{ };
    print $FILEHANDLE "\n\n";

    my $binary_regexp =
      catdir( $parameter_href->{conda_prefix_path} . qw(\\root\\bin) );
    unless ( $ENV{PATH} =~ /$binary_regexp/ ) {

        ## Export path
        print $FILEHANDLE "## Export path\n";
        print $FILEHANDLE q{echo 'source }
          . catdir( $parameter_href->{conda_prefix_path},
            qw(root bin thisroot.sh) )
          . q{' >> ~/.bashrc};
        print $FILEHANDLE "\n\n";

        ## Use newly installed root
        print $FILEHANDLE 'source '
          . catdir( $parameter_href->{conda_prefix_path},
            qw(root bin thisroot.sh) );
        print $FILEHANDLE "\n\n";
    }

    ## Go back to subroutine origin
    print $FILEHANDLE '## Moving back to original working directory', "\n";
    gnu_cd(
        {
            directory_path => $pwd,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Install CNVNator
    print $FILEHANDLE '### Install cnvnator', "\n";

    ## Only activate conda environment if supplied by user
    if ( $parameter_href->{conda_environment} ) {
        ## Activate conda environment
        say $FILEHANDLE q{## Activate conda environment};
        conda_source_activate(
            {
                FILEHANDLE => $FILEHANDLE,
                env_name   => $parameter_href->{conda_environment},
            }
        );
        say $FILEHANDLE $NEWLINE;
    }

    ## Create the temporary install directory
    create_install_dir( { FILEHANDLE => $FILEHANDLE, } );

    ## Download
    print $FILEHANDLE '## Download CNVNator', "\n";
    wget(
        {
            url => 'https://github.com/abyzovlab/CNVnator/releases/download/v'
              . $parameter_href->{cnvnator}
              . '/CNVnator_v'
              . $parameter_href->{cnvnator} . '.zip',
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $parameter_href->{quiet},
            verbose      => $parameter_href->{verbose},
            outfile_path => 'CNVnator_v' . $parameter_href->{cnvnator} . '.zip',
        }
    );
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE '## Extract', "\n";
    print $FILEHANDLE 'unzip CNVnator_v' . $parameter_href->{cnvnator} . '.zip';
    print $FILEHANDLE "\n\n";

    ## Move to CNVnator directory
    print $FILEHANDLE '## Move to CNVnator directory', "\n";
    gnu_cd(
        {
            directory_path => catdir(
                'CNVnator_v' . $parameter_href->{cnvnator}, 'src',
                'samtools'
            ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Configure
    print $FILEHANDLE '## Configure CNVnator samTools specific version', "\n";
    print $FILEHANDLE './configure --without-curses',                    "\n";
    print $FILEHANDLE 'make ';
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE '## Move to CNVnator directory', "\n";
    gnu_cd(
        {
            directory_path => '..',
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n";

    print $FILEHANDLE 'make OMP=no';
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE '## Make available from conda environment', "\n";
    gnu_mv(
        {
            infile_path => 'cnvnator',
            outfile_path =>
              catdir( $parameter_href->{conda_prefix_path}, 'bin' ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE '## Make available from conda environment', "\n";
    gnu_cd(
        {
            directory_path => '..',
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n";
    gnu_mv(
        {
            infile_path => 'cnvnator2VCF.pl',
            outfile_path =>
              catdir( $parameter_href->{conda_prefix_path}, 'bin' ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE '## Make executable from conda environment', "\n";
    print $FILEHANDLE 'chmod +x ';
    print $FILEHANDLE catfile( $parameter_href->{conda_prefix_path},
        qw(bin cnvnator2VCF.pl) );
    print $FILEHANDLE "\n\n";

    ## Remove the temporary install directory
    remove_install_dir(
        {
            FILEHANDLE => $FILEHANDLE,
            pwd        => $pwd,
        }
    );

    ## Deactivate conda environment if conda_environment exists
    if ( $parameter_href->{conda_environment} ) {
        say $FILEHANDLE q{## Deactivate conda environment};
        conda_source_deactivate(
            {
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say $FILEHANDLE $NEWLINE;
    }

    return;
}

sub tiddit {

##tiddit

##Function : Install tiddit
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE, $quiet, $verbose
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (
        check_conda_bin_file_exists(
            {
                parameter_href => $parameter_href,
                program_name   => 'TIDDIT',
            }
        )
      )
    {

        return;
    }

    ## Install tiddit
    print $FILEHANDLE '### Install tiddit', "\n";

    ## Only activate conda environment if supplied by user
    if ( $parameter_href->{conda_environment} ) {
        ## Activate conda environment
        say $FILEHANDLE q{## Activate conda environment};
        conda_source_activate(
            {
                FILEHANDLE => $FILEHANDLE,
                env_name   => $parameter_href->{conda_environment},
            }
        );
        say $FILEHANDLE $NEWLINE;
    }

    ## Move to miniconda environment
    gnu_cd(
        {
            directory_path => catdir( $parameter_href->{conda_prefix_path} ),
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Download
    print $FILEHANDLE '## Download Tiddit', "\n";
    wget(
        {
            url => 'https://github.com/J35P312/TIDDIT/archive/'
              . $parameter_href->{tiddit} . '.zip',
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $parameter_href->{quiet},
            verbose      => $parameter_href->{verbose},
            outfile_path => 'TIDDIT-' . $parameter_href->{tiddit} . '.zip',
        }
    );

    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE '## Extract', "\n";
    gnu_rm(
        {
            infile_path => 'TIDDIT-' . $parameter_href->{tiddit},
            force       => 1,
            recursive   => 1,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";
    print $FILEHANDLE 'unzip TIDDIT-' . $parameter_href->{tiddit} . '.zip',
      "\n\n";

    ## Move to Tiddit directory
    print $FILEHANDLE '## Move to Tiddit directory', "\n";
    gnu_cd(
        {
            directory_path => 'TIDDIT-' . $parameter_href->{tiddit},
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    gnu_mkdir(
        {
            indirectory_path => 'build',
            parents          => 1,
            FILEHANDLE       => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    gnu_cd(
        {
            directory_path => 'build',
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Configure
    print $FILEHANDLE '## Configure', "\n";
    print $FILEHANDLE 'cmake .. ',    "\n\n";

    print $FILEHANDLE 'make', "\n\n";

    gnu_cd(
        {
            directory_path => catdir(qw(.. bin)),
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";
    print $FILEHANDLE 'chmod a+x TIDDIT', "\n\n";

    ## Make available from conda environment
    my $cwd = cwd();
    print $FILEHANDLE '## Make available from conda environment', "\n";

    ## Specifying target and link_path
    my $target_path = catfile(
        $parameter_href->{conda_prefix_path},
        q{TIDDIT-} . $parameter_href->{tiddit},
        qw{ bin TIDDIT }
    );
    my $link_path =
      catfile( $parameter_href->{conda_prefix_path}, qw{ bin TIDDIT } );
    gnu_ln(
        {
            FILEHANDLE  => $FILEHANDLE,
            target_path => $target_path,
            link_path   => $link_path,
            symbolic    => 1,
            force       => 1,
        }
    );

    print $FILEHANDLE $NEWLINE;

    ## Clean-up
    gnu_rm(
        {
            infile_path => catdir(
                $parameter_href->{conda_prefix_path},
                'Tiddit-' . $parameter_href->{tiddit} . '.zip'
            ),
            force      => 1,
            FILEHANDLE => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Go back to subroutine origin
    print $FILEHANDLE '## Moving back to original working directory', "\n";
    gnu_cd(
        {
            directory_path => $pwd,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Deactivate conda environment if conda_environment exists
    if ( $parameter_href->{conda_environment} ) {
        say $FILEHANDLE q{## Deactivate conda environment};
        conda_source_deactivate(
            {
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say $FILEHANDLE $NEWLINE;
    }

    return;
}

sub svdb {

##svdb

##Function : Install svdb
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE, $quiet, $verbose
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (
        check_conda_bin_file_exists(
            {
                parameter_href => $parameter_href,
                program_name   => 'svdb',
            }
        )
      )
    {

        return;
    }

    ## Install svdb
    print $FILEHANDLE '### Install svdb', "\n";

    ## Only activate conda environment if supplied by user
    if ( $parameter_href->{conda_environment} ) {
        ## Activate conda environment
        say $FILEHANDLE q{## Activate conda environment};
        conda_source_activate(
            {
                FILEHANDLE => $FILEHANDLE,
                env_name   => $parameter_href->{conda_environment},
            }
        );
        say $FILEHANDLE $NEWLINE;
    }

    ## Move to miniconda environment
    gnu_cd(
        {
            directory_path => catdir( $parameter_href->{conda_prefix_path} ),
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Download
    print $FILEHANDLE '## Download Svdb', "\n";
    wget(
        {
            url => 'https://github.com/J35P312/SVDB/archive/SVDB-'
              . $parameter_href->{svdb} . '.zip',
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $parameter_href->{quiet},
            verbose      => $parameter_href->{verbose},
            outfile_path => 'SVDB-' . $parameter_href->{svdb} . '.zip',
        }
    );

    print $FILEHANDLE "\n\n";

    ## Clean-up
    print $FILEHANDLE '## Clean-up', "\n";
    gnu_rm(
        {
            infile_path => 'SVDB-' . $parameter_href->{svdb},
            force       => 1,
            recursive   => 1,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE '## Extract', "\n";
    print $FILEHANDLE 'unzip SVDB-' . $parameter_href->{svdb} . '.zip', "\n\n";

    ## Move to Svdb directory
    print $FILEHANDLE '## Move to Svdb directory', "\n";
    gnu_cd(
        {
            directory_path => 'SVDB-SVDB-' . $parameter_href->{svdb},
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Pip install the downloaded SVDB package
    say {$FILEHANDLE} q{## Install};
    pip_install(
        {
            packages_ref => [$DOT],
            quiet        => $parameter_href->{quiet},
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Go back to subroutine origin
    print $FILEHANDLE '## Moving back to original working directory', "\n";
    gnu_cd(
        {
            directory_path => $pwd,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Deactivate conda environment if conda_environment exists
    if ( $parameter_href->{conda_environment} ) {
        say $FILEHANDLE q{## Deactivate conda environment};
        conda_source_deactivate(
            {
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say $FILEHANDLE $NEWLINE;
    }

    return;
}

sub mip_scripts {

##mip_scripts

##Function : Install mip_scripts
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Define MIP scripts and yaml files
    my @mip_scripts = (
        'calculate_af.pl', 'download_reference.pl',
        'mip_install.pl',  'max_af.pl',
        'mip.pl',          'qccollect.pl',
        'vcfparser.pl',    'covplots_genome.R'
    );

    my %mip_sub_scripts;
    $mip_sub_scripts{definitions} =
      [qw(define_download_references.yaml define_parameters.yaml)];
    $mip_sub_scripts{t} = [qw(mip_install.t mip.t run_tests.t mip_analysis.t)];
    $mip_sub_scripts{templates} = [qw(mip_config.yaml)];

    my @mip_directories = qw(lib t);

    ## Check if the binary of the program being installed already exists
    if (
        check_conda_bin_file_exists(
            {
                parameter_href => $parameter_href,
                program_name   => 'mip.pl',
            }
        )
      )
    {    #Proxy for all

        return;
    }

    ## Install mip_scripts
    print $FILEHANDLE '### Install mip_scripts', "\n";

    ## Create directories
    print $FILEHANDLE '## Create directories', "\n";
    foreach my $directory ( keys %mip_sub_scripts ) {

        gnu_mkdir(
            {
                indirectory_path => catdir(
                    $parameter_href->{conda_prefix_path},
                    'bin', $directory
                ),
                parents    => 1,
                FILEHANDLE => $FILEHANDLE,
            }
        );
        print $FILEHANDLE "\n\n";
    }
    ## Copy directory to conda env
    print $FILEHANDLE '## Copy directory to conda env', "\n\n";
    foreach my $directory (@mip_directories) {

        gnu_cp(
            {
                FILEHANDLE  => $FILEHANDLE,
                recursive   => 1,
                force       => 1,
                infile_path => catdir( $Bin, $directory ),
                outfile_path =>
                  catdir( $parameter_href->{conda_prefix_path}, 'bin' ),
            }
        );
        print $FILEHANDLE "\n\n";
    }

    ## Copy mip scripts and sub scripts to conda env and make executable
    print $FILEHANDLE
'## Copy mip scripts and subdirectory scripts to conda env and make executable',
      "\n";
    foreach my $script (@mip_scripts) {

        my $script_no_ending = fileparse( $script, qr/\.[^.]*/ );

        gnu_cp(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => catfile( $Bin, $script ),
                outfile_path => catdir(
                    $parameter_href->{conda_prefix_path}, 'bin',
                    $script_no_ending
                ),
            }
        );
        print $FILEHANDLE "\n";

        print $FILEHANDLE 'chmod a+x '
          . catfile( $parameter_href->{conda_prefix_path},
            'bin', $script_no_ending );
        print $FILEHANDLE "\n\n";
    }

    foreach my $directory ( keys %mip_sub_scripts ) {

        foreach my $script ( @{ $mip_sub_scripts{$directory} } ) {

            gnu_cp(
                {
                    FILEHANDLE   => $FILEHANDLE,
                    infile_path  => catfile( $Bin, $directory, $script ),
                    outfile_path => catdir(
                        $parameter_href->{conda_prefix_path}, 'bin',
                        $directory
                    ),
                }
            );
            print $FILEHANDLE "\n";

            print $FILEHANDLE 'chmod a+x '
              . catfile( $parameter_href->{conda_prefix_path},
                'bin', $directory, $script );
            print $FILEHANDLE "\n\n";
        }
    }
    return;
}

sub rhocall {

##rhocall

##Function : Install rhocall
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (
        check_conda_bin_file_exists(
            {
                parameter_href => $parameter_href,
                program_name   => 'rhocall',
            }
        )
      )
    {

        return;
    }

    ## Only activate conda environment if supplied by user
    if ( $parameter_href->{conda_environment} ) {
        ## Activate conda environment
        say $FILEHANDLE q{## Activate conda environment};
        conda_source_activate(
            {
                FILEHANDLE => $FILEHANDLE,
                env_name   => $parameter_href->{conda_environment},
            }
        );
        say $FILEHANDLE $NEWLINE;
    }

    ## Install rhocall
    print $FILEHANDLE '### Install rhocall', "\n";

    ## Create the temporary install directory
    create_install_dir(
        {
            FILEHANDLE        => $FILEHANDLE,
            install_directory => $parameter_href->{rhocall_path},
        }
    );

    ## Downloads files
    print $FILEHANDLE '## Download rhocall', "\n";
    wget(
        {
            url => 'https://github.com/dnil/rhocall/archive/'
              . $parameter_href->{rhocall} . '.zip',
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $parameter_href->{quiet},
            verbose      => $parameter_href->{verbose},
            outfile_path => 'rhocall-' . $parameter_href->{rhocall} . '.zip',
        }
    );
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE '## Extract', "\n";
    print $FILEHANDLE 'unzip rhocall-' . $parameter_href->{rhocall} . '.zip';
    print $FILEHANDLE "\n\n";

    ## Move to rhocall directory
    print $FILEHANDLE '## Move to rhocall directory', "\n";
    gnu_cd(
        {
            directory_path => 'rhocall-' . $parameter_href->{rhocall},
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Configure
    say {$FILEHANDLE} q{## Configure};
    pip_install(
        {
            packages_ref => [qw{ numpy Cython }],
            quiet        => $parameter_href->{quiet},
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    print {$FILEHANDLE} $NEWLINE;
    pip_install(
        {
            requirement => q{requirements.txt},
            quiet       => $parameter_href->{quiet},
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    print {$FILEHANDLE} $NEWLINE;
    pip_install(
        {
            editable   => $DOT,
            quiet      => $parameter_href->{quiet},
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Go back to subroutine origin
    print $FILEHANDLE '## Moving back to original working directory', "\n";
    gnu_cd(
        {
            directory_path => $pwd,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Deactivate conda environment if conda_environment exists
    if ( $parameter_href->{conda_environment} ) {
        say $FILEHANDLE q{## Deactivate conda environment};
        conda_source_deactivate(
            {
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say $FILEHANDLE $NEWLINE;
    }

    return;
}

sub remove_install_dir {

##remove_install_dir

##Function : Remove the temporary install directory
##Returns  : ""
##Arguments: $FILEHANDLE, $pwd, $install_directory
##         : $FILEHANDLE        => FILEHANDLE to write to
##         : $pwd               => The original working directory
##         : $install_directory => Temporary installation directory

    my ($arg_href) = @_;

    ## Default(s)
    my $install_directory;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $pwd;

    my $tmpl = {
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
        pwd =>
          { required => 1, defined => 1, strict_type => 1, store => \$pwd },
        install_directory => {
            default     => '.MIP',
            allow       => qr/^\.\S+$/,
            strict_type => 1,
            store       => \$install_directory
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    ## Go back to subroutine origin
    print $FILEHANDLE '## Moving back to original working directory', "\n";
    gnu_cd(
        {
            directory_path => $pwd,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Clean up
    print $FILEHANDLE '## Clean up', "\n";
    gnu_rm(
        {
            infile_path => $install_directory,
            force       => 1,
            recursive   => 1,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    return;
}

sub create_install_dir {

##create_install_dir

##Function : Create the temporary install directory
##Returns  : ""
##Arguments: $FILEHANDLE, $install_directory
##         : $FILEHANDLE        => FILEHANDLE to write to
##         : $install_directory => Temporary installation directory

    my ($arg_href) = @_;

    ## Default(s)
    my $install_directory;

    ## Flatten argument(s)
    my $FILEHANDLE;

    my $tmpl = {
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
        install_directory => {
            default     => '.MIP',
            strict_type => 1,
            store       => \$install_directory
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    ## Create temp install directory
    print $FILEHANDLE '## Create temp install directory', "\n";
    gnu_mkdir(
        {
            indirectory_path => $install_directory,
            parents          => 1,
            FILEHANDLE       => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    gnu_cd(
        {
            directory_path => $install_directory,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    return;
}

sub check_conda_bin_file_exists {

##check_conda_bin_file_exists

##Function : Check if the binary of the program being installed already exists
##Returns  : ""
##Arguments: $parameter_href, $program_name, $program_version
##         : $parameter_href  => Holds all parameters
##         : $program_name    => Program name
##         : $program_version => Program version

    my ($arg_href) = @_;

    ## Default(s)
    my $program_version;

    ## Flatten argument(s)
    my $parameter_href;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name
        },
        program_version => {
            default     => 0,
            strict_type => 1,
            store       => \$program_version
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $miniconda_bin_file;

    if ($program_version) {

        $miniconda_bin_file = catfile( $parameter_href->{conda_prefix_path},
            'bin', $program_name . $program_version );
    }
    else {

        $miniconda_bin_file =
          catfile( $parameter_href->{conda_prefix_path}, 'bin', $program_name );
    }

    if ( -f $miniconda_bin_file ) {

        if ($program_version) {

            $log->info( q{Found}
                  . $SPACE
                  . $program_name
                  . $SPACE
                  . q{version}
                  . $SPACE
                  . $program_version
                  . $SPACE
                  . q{in miniconda directory:}
                  . $SPACE
                  . catdir( $parameter_href->{conda_prefix_path}, q{bin} ) );

            if ( $parameter_href->{noupdate} ) {

                $log->info( q{Skipping writting installation process for}
                      . $SPACE
                      . $program_name
                      . $SPACE
                      . $program_version );
                return 1;
            }
            $log->info(
                q{Writting install instructions for} . $SPACE . $program_name );
        }
        else {

            $log->info( q{Found}
                  . $SPACE
                  . $program_name
                  . $SPACE
                  . q{in miniconda directory:}
                  . $SPACE
                  . catdir( $parameter_href->{conda_prefix_path}, q{bin} ) );

            if ( $parameter_href->{noupdate} ) {

                $log->info( q{Skipping writting installation process for}
                      . $SPACE
                      . $program_name );
                return 1;
            }
            $log->info(
                q{Writting install instructions for} . $SPACE . $program_name );
        }
        return 0;
    }
    else {

        $log->info(
            q{Writting install instructions for} . $SPACE . $program_name );
        return 0;
    }
    return;
}

sub references {

##references

##Function : Install references
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Only activate conda environment if supplied by user
    if ( $parameter_href->{conda_environment} ) {
        ## Activate conda environment
        say $FILEHANDLE q{## Activate conda environment};
        conda_source_activate(
            {
                FILEHANDLE => $FILEHANDLE,
                env_name   => $parameter_href->{conda_environment},
            }
        );
        say $FILEHANDLE $NEWLINE;
    }

    $log->info(q{Writting install instructions for references});

    print $FILEHANDLE 'download_reference ';
    print $FILEHANDLE '--reference_dir '
      . $parameter_href->{reference_dir} . q{ };

    print $FILEHANDLE '--reference_genome_versions '
      . join( ' --reference_genome_versions ',
        @{ $parameter_href->{reference_genome_versions} } )
      . q{ };
    print $FILEHANDLE "\n\n";

    ##Launch bash
    print $FILEHANDLE 'bash download_reference.sh', "\n\n";

    ##Cleanup
    gnu_rm(
        {
            infile_path => 'download_reference.sh',
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Deactivate conda environment if conda_environment exists
    if ( $parameter_href->{conda_environment} ) {
        say $FILEHANDLE q{## Deactivate conda environment};
        conda_source_deactivate(
            {
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say $FILEHANDLE $NEWLINE;
    }

    return;
}

