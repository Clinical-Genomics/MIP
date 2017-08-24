#!/usr/bin/env perl

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use English qw(-no_match_vars);
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

use Getopt::Long;
use Cwd;
use Cwd qw(abs_path);
use FindBin qw($Bin);                 #Find directory of script
use IO::Handle;
use File::Basename qw(dirname basename fileparse);
use File::Spec::Functions qw(catfile catdir devnull);
use Readonly;

## MIPs lib/
use lib catdir( $Bin, 'lib' );        #Add MIPs internal lib
use MIP::Language::Shell qw(create_bash_file);
use Program::Download::Wget qw(wget);
use MIP::Gnu::Bash qw(gnu_cd);
use MIP::Gnu::Coreutils qw(gnu_cp gnu_rm gnu_mv gnu_mkdir gnu_link );
use MIP::PacketManager::Conda qw{ conda_create conda_source_activate conda_source_deactivate conda_update conda_check };
use Script::Utils qw(help set_default_array_parameters);


our $USAGE = build_usage( {} );

### Set parameter default

my %parameter;

##Bash
$parameter{bash_set_errexit} = 0;
$parameter{bash_set_nounset} = 0;

## Conda
$parameter{conda_dir_path} = catdir( $ENV{HOME}, 'miniconda' );
$parameter{python_version} = '2.7';

## Bioconda channel
$parameter{bioconda}{bwa}       = '0.7.15';
$parameter{bioconda}{bwakit}    = '0.7.12';
$parameter{bioconda}{fastqc}    = '0.11.5';
$parameter{bioconda}{cramtools} = '3.0.b47';
$parameter{bioconda}{samtools}  = '1.4.1';
$parameter{bioconda}{bcftools}  = '1.4.1';
$parameter{bioconda}{snpeff}    = '4.3.1';
$parameter{bioconda}{snpsift}   = '4.3.1';
$parameter{bioconda}{picard}    = '2.9.2';
$parameter{bioconda}{htslib}    = '1.4.1';
$parameter{bioconda}{bedtools}  = '2.26.0';
$parameter{bioconda}{vt}        = '2015.11.10';
$parameter{bioconda}{sambamba}  = '0.6.6';
$parameter{bioconda}{freebayes} = '1.1.0';
$parameter{bioconda}{delly}     = '0.7.7';
$parameter{bioconda}{manta}     = '1.1.0';
$parameter{bioconda}{multiqc}   = '0.9.1a0';
$parameter{bioconda}{peddy}     = '0.2.9';
$parameter{bioconda}{plink2}    = '1.90b3.35';
$parameter{bioconda}{vcfanno}   = '0.1.0';

# Required for CNVnator
$parameter{bioconda}{gcc}   = '4.8.5';
$parameter{bioconda}{cmake} = '3.3.1';

## Bioconda pathes
# For correct softlinking in share and bin in conda env
$parameter{bioconda_bwakit_patch}  = '-0';
$parameter{bioconda_snpeff_patch}  = 'p-1';
$parameter{bioconda_snpsift_patch} = 'p-0';
$parameter{bioconda_picard_patch}  = '-1';
$parameter{bioconda_manta_patch}   = '-0';

## Perl Modules
$parameter{perl_version} = '5.18.2';

## PIP
$parameter{pip}{genmod}            = '3.7.1';
$parameter{pip}{variant_integrity} = '0.0.4';
$parameter{pip}{chanjo}            = '4.0.0';

## Programs currently installable by SHELL
$parameter{mip_scripts}            = 'Your current MIP version';
$parameter{picardtools}            = '2.3.0';
$parameter{sambamba}               = '0.6.1';
$parameter{bedtools}               = '2.25.0';
$parameter{vt}                     = 'gitRepo';
$parameter{plink2}                 = '160316';
$parameter{snpeff}                 = 'v4_2';
$parameter{varianteffectpredictor} = '88.8';
$parameter{vep_auto_flag}          = 'alcf';
$parameter{rhocall}                = '0.4';
$parameter{rhocall_path}           = catdir( $ENV{HOME}, 'rhocall' );
$parameter{cnvnator}               = '0.3.3';
$parameter{cnvnator_root_binary} =
  'root_v6.06.00.Linux-slc6-x86_64-gcc4.8.tar.gz';
$parameter{tiddit} = '1.1.5';
$parameter{svdb}   = '1.0.6';

## Define default parameters
my %array_parameter;

$array_parameter{vep_assemblies}{default} = [qw(GRCh37 GRCh38)];
$array_parameter{vep_plugins}{default}    = [qw(UpDownDistance LoFtool Lof)];

# GRCh38.86 but check current on the snpEff sourceForge
$array_parameter{snpeff_genome_versions}{default} =
  [qw(GRCh37.75 GRCh38.86)];
$array_parameter{reference_genome_versions}{default} = [qw(GRCh37 hg38)];
$array_parameter{perl_modules}{default}              = [
    'Modern::Perl',              # MIP
    'IPC::System::Simple',       # MIP
    'Path::Iterator::Rule',      # MIP
    'YAML',                      # MIP
    'Log::Log4perl',             # MIP
    'List::Util',                # MIP
    'List::MoreUtils',           # MIP
    'Readonly',                  # MIP
    'Scalar::Util::Numeric',     # MIP
    'Set::IntervalTree',         # MIP/vcfParser.pl
    'Net::SSLeay',               # VEP
    'LWP::Simple',               # VEP
    'LWP::Protocol::https',      # VEP
    'PerlIO::gzip',              # VEP
    'IO::Uncompress::Gunzip',    # VEP
    'HTML::Lint',                # VEP
    'Archive::Zip',              # VEP
    'Archive::Extract',          # VEP
    'DBI',                       # VEP
    'JSON',                      # VEP
    'DBD::mysql',                # VEP
    'CGI',                       # VEP
    'Sereal::Encoder',           # VEP
    'Sereal::Decoder',           # VEP
    'Bio::Root::Version',        # VEP
    'Module::Build',             # VEP
    'File::Copy::Recursive',     # VEP
];

my $VERSION = '1.2.7';

###User Options
GetOptions(
    'see|bash_set_errexit'          => \$parameter{bash_set_errexit},
    'snu|bash_set_nounset'          => \$parameter{bash_set_nounset},
    'env|conda_environment:s'       => \$parameter{conda_environment},
    'cdp|conda_dir_path:s'          => \$parameter{conda_dir_path},
    'cdu|conda_update'              => \$parameter{conda_update},
    'bcv|bioconda=s'                => \%{ $parameter{bioconda} },
    'pip|pip=s'                     => \%{ $parameter{pip} },
    'pyv|python_version=s'          => \$parameter{python_version},
    'pev|perl_version=s'            => \$parameter{perl_version},
    'pei|perl_install'              => \$parameter{perl_install},
    'pevs|perl_skip_test'           => \$parameter{perl_skip_test},
    'pm|perl_modules:s'             => \@{ $parameter{perl_modules} },
    'pmf|perl_modules_force'        => \$parameter{perl_modules_force},
    'pic|picardtools:s'             => \$parameter{picardtools},
    'sbb|sambamba:s'                => \$parameter{sambamba},
    'bet|bedtools:s'                => \$parameter{bedtools},
    'vt|vt:s'                       => \$parameter{vt},
    'plk|plink2:s'                  => \$parameter{plink2},
    'snpg|snpeff_genome_versions:s' => \@{ $parameter{snpeff_genome_versions} },
    'vep|varianteffectpredictor:s'  => \$parameter{varianteffectpredictor},
    'vepai|vep_auto_flag:s'         => \$parameter{vep_auto_flag},
    'vepc|vep_cache_dir:s'          => \$parameter{vep_cache_dir},
    'vepa|vep_assemblies:s'         => \@{ $parameter{vep_assemblies} },
    'vepp|vep_plugins:s'            => \@{ $parameter{vep_plugins} },
    'rhc|rhocall:s'                 => \$parameter{rhocall},
    'rhcp|rhocall_path:s'           => \$parameter{rhocall_path},
    'cnv|cnvnator:s'                => \$parameter{cnvnator},
    'cnvnr|cnvnator_root_binary:s'  => \$parameter{cnvnator_root_binary},
    'tid|tiddit:s'                  => \$parameter{tiddit},
    'svdb|svdb:s'                   => \$parameter{svdb},
    'psh|prefer_shell'              => \$parameter{prefer_shell},
    'ppd|print_parameters_default'  => sub {
        print_parameters(
            {
                parameter_href       => \%parameter,
                array_parameter_href => \%array_parameter,
            }
        );
        exit;
    },    # Display parameter defaults
    'nup|noupdate'         => \$parameter{noupdate},
    'sp|select_programs:s' => \@{ $parameter{select_programs} },
    'rd|reference_dir:s'   => \$parameter{reference_dir},
    'rg|reference_genome_versions:s' =>
      \@{ $parameter{reference_genome_versions} },
    'q|quiet' => \$parameter{quiet},
    'h|help'  => sub {
        print STDOUT $USAGE, "\n";
        exit;
    },    #Display help text
    'ver|version' => sub {
        print STDOUT "\n" . basename($PROGRAM_NAME) . q{ } . $VERSION, "\n\n";
        exit;
    },    #Display version number
    'v|verbose' => \$parameter{verbose},
  )
  or croak Script::Utils::help(
    {
        USAGE     => $USAGE,
        exit_code => 1,
    }
  );

## Update default parameter dependent on other parameters
if (   ( exists $parameter{conda_environment} )
    && ( $parameter{conda_environment} ) )
{

    $parameter{conda_prefix_path} = catdir( $parameter{conda_dir_path},
        'envs', $parameter{conda_environment} );
}
else {

    $parameter{conda_prefix_path} = $parameter{conda_dir_path};
}

if ( !$parameter{vep_cache_dir} ) {

    # Cache directory
    $parameter{vep_cache_dir} = catdir( $parameter{conda_prefix_path},
        'ensembl-tools-release-' . $parameter{varianteffectpredictor},
        'cache' );
}

## Set default for array parameters
Script::Utils::set_default_array_parameters(
    {
        parameter_href       => \%parameter,
        array_parameter_href => \%array_parameter,
    }
);

## Constants
Readonly my $SPACE => q{ };
Readonly my $NEWLINE => qq{\n};

##########
###MAIN###
##########

# Create anonymous filehandle
my $FILEHANDLE = IO::Handle->new();

# Installation instruction file
my $file_name_path = catfile( cwd(), 'mip.sh' );

open $FILEHANDLE, '>', $file_name_path
  or
  croak( q{Cannot write to '} . $file_name_path . q{' :} . $OS_ERROR . "\n" );

## Create bash file for writing install instructions
create_bash_file(
    {
        file_name   => $file_name_path,
        FILEHANDLE  => $FILEHANDLE,
        remove_dir  => catfile( cwd(), '.MIP' ),
        set_errexit => $parameter{bash_set_errexit},
        set_nounset => $parameter{bash_set_nounset},
    }
);

print STDOUT q{Will write install instructions to '} . $file_name_path, "'\n";

## Check existance of conda environment
conda_check(
    {
        conda_dir_path  => $parameter{conda_dir_path},
    }
);

## Optionally update conda
if ( $parameter{conda_update} ) {
    say $FILEHANDLE q{### Updating Conda};
    conda_update(
        {
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say $FILEHANDLE $NEWLINE;
}


if ( exists( $parameter{conda_environment} ) ) {

    ## Check conda environment
    if ( !-d catdir( $parameter{conda_prefix_path} ) ) {

        ## Create conda environment and install pip
        say $FILEHANDLE q{## Creating conda environment: } 
          . $parameter{conda_environment} 
          . q{and install packages}; 
        conda_create(
            {
                env_name => $parameter{conda_environment},
                python_version => $parameter{python_version},
                packages_ref => [ qw{pip} ],
                FILEHANDLE     => $FILEHANDLE,
            }
        );
        say $FILEHANDLE $NEWLINE;
    }
}

## Install modules into conda environment using channel Bioconda
install_bioconda_modules(
    {
        parameter_href => \%parameter,
        FILEHANDLE     => $FILEHANDLE,
    }
);

if ( @{ $parameter{select_programs} } ) {

    if ( ( grep { $_ eq 'perl' } @{ $parameter{select_programs} } ) )
    {    #If element is part of array

        perl(
            {
                parameter_href => \%parameter,
                FILEHANDLE     => $FILEHANDLE,
            }
        );
    }
}
else {

    perl(
        {
            parameter_href => \%parameter,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
}

pip_install(
    {
        parameter_href => \%parameter,
        FILEHANDLE     => $FILEHANDLE,
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
            grep { $_ eq 'varianteffectpredictor' }
            @{ $parameter{select_programs} }
        )
      )
    {    #If element is part of array

        varianteffectpredictor(
            {
                parameter_href => \%parameter,
                FILEHANDLE     => $FILEHANDLE,
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

    varianteffectpredictor(
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

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    return <<"END_USAGE";
 $script_name [options]
    -env/--conda_environment conda environment (Default: "")
    -cdp/--conda_dir_path The conda directory path (Default: "HOME/miniconda")
    -cdu/--conda_update Update conda before installing (Supply flag to enable)
    -bvc/--bioconda Set the module version of the programs that can be installed with bioconda (e.g. 'bwa=0.7.12')
    -pip/--pip Set the module version of the programs that can be installed with pip (e.g. 'genmod=3.7.01')
    -pyv/--python_version Set the env python version (Default: "2.7")

    ## SHELL
    -pei/--perl_install Install perl (Supply flag to enable)
    -pev/--perl_version Set the perl version (defaults: "5.18.2")
    -pevs/--perl_skip_test Skip "tests" in perl installation
    -pm/--perl_modules Set the perl modules to be installed via cpanm (Default: ["Modern::Perl", "List::Util", "IPC::System::Simple", "Path::Iterator::Rule", "YAML", "Log::Log4perl", "Set::IntervalTree", "Net::SSLeay",P, "LWP::Simple", "LWP::Protocol::https", "Archive::Zip", "Archive::Extract", "DBI","JSON", "DBD::mysql", "CGI", "Sereal::Encoder", "Sereal::Decoder", "Bio::Root::Version", "Module::Build"])
    -pmf/--perl_modules_force Force installation of perl modules
    -pic/--picardtools Set the picardtools version (Default: "2.5.9"),
    -sbb/--sambamba Set the sambamba version (Default: "0.6.6")
    -bet/--bedtools Set the bedtools version (Default: "2.26.0")
    -vt/--vt Set the vt version (Default: "0.57")
    -plk/--plink  Set the plink version (Default: "160224")
    -snpg/--snpeff_genome_versions Set the snpEff genome version (Default: ["GRCh37.75", "GRCh38.82"])
    -vep/--varianteffectpredictor Set the VEP version (Default: "88")
    -vepa/--vep_auto_flag Set the VEP auto installer flags
    -vepc/--vep_cache_dir Specify the cache directory to use (whole path; defaults to "[--conda_dir_path]/ensembl-tools-release-varianteffectpredictorVersion/cache")
    -vepa/--vep_assemblies Select the assembly version (Default: ["GRCh37", "GRCh38"])
    -vepp/--vep_plugins Supply VEP plugins (Default: "UpDownDistance, LoFtool, Lof")
    -rhc/--rhocall Set the rhocall version (Default: "0.4")
    -rhcp/--rhocall_path Set the path to where to install rhocall (Defaults: "HOME/rhocall")
    -cnvn/--cnvnator Set the cnvnator version (Default: 0.3.3)
    -cnvnr/--cnvnator_root_binary Set the cnvnator root binary (Default: "root_v6.06.00.Linux-slc6-x86_64-gcc4.8.tar.gz")
    -tid/--tiddit Set the tiddit version (Default: "1.1.5")
    -svdb/--svdb Set the svdb version (Default: "1.0.6")

    ## Utility
    -psh/--prefer_shell Shell will be used for overlapping shell and biconda installations (Supply flag to enable)
    -ppd/--print_parameters_default Print the parameter defaults
    -nup/--noupdate Do not update already installed programs (Supply flag to enable)
    -sp/--select_programs Install supplied programs e.g. -sp perl -sp bedtools (Default: "")
    -rd/--reference_dir Reference(s) directory (Default: "")
    -rd/--reference_genome_versions Reference versions to download ((Default: ["GRCh37", "hg38"]))
    -q/--quiet Quiet (Supply flag to enable; no output from individual program that has a quiet flag)
    -h/--help Display this help message
    -ver/--version Display version
    -v/--verbose Set verbosity
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

    ## Set default for array parameters
    Script::Utils::set_default_array_parameters(
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

sub install_bioconda_modules {

##install_bioconda_modules

##Function : Install modules into conda environment using channel Bioconda
##Returns  : ""
##Arguments: $parameter_href
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

    if ( exists( $parameter_href->{conda_environment} )
        && ( $parameter_href->{conda_environment} ) )
    {

        ## Install into conda environment using bioconda channel
        print $FILEHANDLE '### Installing into conda environment: '
          . $parameter_href->{conda_environment}, "\n";
    }
    else {

        print $FILEHANDLE '### Installing into conda main environment', "\n";
    }
    print $FILEHANDLE 'conda install ';

    if ( $parameter_href->{quiet} ) {

        print $FILEHANDLE '--quiet ';    #Do not display progress bar
    }

    if ( exists( $parameter_href->{conda_environment} )
        && ( $parameter_href->{conda_environment} ) )
    {

        print $FILEHANDLE '-n ' . $parameter_href->{conda_environment} . q{ };
    }
    print $FILEHANDLE '-y ';
    print $FILEHANDLE '-c bioconda ';

    ## Install all bioconda packages
    foreach my $program ( keys %{ $parameter_href->{bioconda} } ) {

        print $FILEHANDLE $program . '='
          . $parameter_href->{bioconda}{$program} . q{ };
    }

    print $FILEHANDLE "\n\n";

    ## Custom
    foreach my $program ( keys %{ $parameter_href->{bioconda} } ) {

        if ( $program eq 'bwakit' ) {

            ## Define binaries
            my @bwakit_binaries = (
                'k8',                'seqtk',
                'bwa-postalt.js',    'run-HLA',
                'typeHLA.sh',        'fermi2',
                'fermi2.pl',         'ropebwt2',
                'typeHLA-selctg.js', 'typeHLA.js'
            );

            foreach my $binary (@bwakit_binaries) {
                # Specifying target and link paths
                my $target_path = catfile(
                    $parameter_href->{conda_prefix_path}, q{share},
                    q{bwakit-} . $parameter_href->{bioconda}{bwakit}
                    . $parameter_href->{bioconda_bwakit_patch}, $binary
                );
                my $link_path = catfile(
                    $parameter_href->{conda_prefix_path}, q{bin}, $binary
                );
                gnu_link(
                    {
                        FILEHANDLE  => $FILEHANDLE,
                        target_path => $target_path,
                        link_path   => $link_path,
                        symbolic    => 1,
                        force       => 1,
                    }
                );
                print $FILEHANDLE $NEWLINE;
            }
            print $FILEHANDLE $NEWLINE;

            gnu_cp(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    recursive   => 1,
                    force       => 1,
                    infile_path => catdir(
                        $parameter_href->{conda_prefix_path},
                        'share',
                        'bwakit-'
                          . $parameter_href->{bioconda}{bwakit}
                          . $parameter_href->{bioconda_bwakit_patch},
                        'resource-human-HLA'
                    ),
                    outfile_path =>
                      catdir( $parameter_href->{conda_prefix_path}, 'bin' ),
                }
            );
            print $FILEHANDLE "\n\n";
        }

        if ( $program eq 'picard' ) {
            # Specifying target and link paths
            my $target_path = catfile(
                $parameter_href->{conda_prefix_path}, q{share}, q{picard-} 
                . $parameter_href->{bioconda}{picard} 
                . $parameter_href->{bioconda_picard_patch}, q{picard.jar}
            );
            my $link_path = catfile(
                $parameter_href->{conda_prefix_path}, q{picard.jar}
            );
            gnu_link(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    target_path => $target_path,
                    link_path   => $link_path,
                    symbolic    => 1,
                    force       => 1,
                }
            );
            print $FILEHANDLE $NEWLINE;
        }

        if ( $program eq 'snpeff' ) {
            ## Define binaries
            my @snpeff_binaries = qw(snpEff.jar snpEff.config);
            
            foreach my $binary (@snpeff_binaries) {
                # Specifying target and link paths
                my $target_path = catfile(
                    $parameter_href->{conda_prefix_path}, q{share}, 
                    q{snpeff-} . $parameter_href->{bioconda}{snpeff} 
                    . $parameter_href->{bioconda_snpeff_patch}, $binary
                );
                my $link_path = catfile(
                    $parameter_href->{conda_prefix_path}, $binary
                );
                gnu_link(
                    {
                        FILEHANDLE  => $FILEHANDLE,
                        target_path => $target_path,
                        link_path   => $link_path,
                        symbolic    => 1,
                        force       => 1,
                    }
                );
                print $FILEHANDLE $NEWLINE;
            }
            print $FILEHANDLE $NEWLINE;

            foreach my $genome_version (
                @{ $parameter_href->{snpeff_genome_versions} } )
            {

                ## Check and if required add the vertebrate mitochondrial codon table to snpeff config
                check_mt_codon_table(
                    {
                        parameter_href => $parameter_href,
                        FILEHANDLE     => $FILEHANDLE,
                        share_dir      => catdir(
                            $parameter_href->{conda_prefix_path},
                            'share',
                            'snpeff-'
                              . $parameter_href->{bioconda}{snpeff}
                              . $parameter_href->{bioconda_snpeff_patch}
                        ),
                        config_file        => 'snpEff.config',
                        genome_version_ref => \$genome_version,
                    }
                );

                unless (
                    -d catdir(
                        $parameter_href->{conda_prefix_path},
                        'share',
                        'snpeff-'
                          . $parameter_href->{bioconda}{snpeff}
                          . $parameter_href->{bioconda_snpeff_patch},
                        'data',
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
        }

        if ( $program eq 'snpsift' ) {
            ## Define binaries
            my @snpsift_binaries = qw(SnpSift.jar);
            
            foreach my $binary (@snpsift_binaries) {
                ## Specifying target and link paths
                my $target_path = catfile(
                    $parameter_href->{conda_prefix_path}, q{share}, 
                    q{snpsift-} . $parameter_href->{bioconda}{snpsift}
                    . $parameter_href->{bioconda_snpsift_patch}, $binary
                );
                my $link_path = catfile(
                    $parameter_href->{conda_prefix_path}, $binary
                );
                gnu_link(
                    {
                        FILEHANDLE  => $FILEHANDLE,
                        target_path => $target_path,
                        link_path   => $link_path,
                        symbolic    => 1,
                        force       => 1,
                    }
                );
                print $FILEHANDLE $NEWLINE;
            }
            print $FILEHANDLE $NEWLINE;
        }

        if ( $program eq 'manta' ) {
            ## Define binaries
            my @manta_binaries = qw(configManta.py configManta.py.ini);
            
            foreach my $binary (@manta_binaries) {
                ## Specifying target and link paths
                my $target_path = catfile(
                    $parameter_href->{conda_prefix_path}, q{share}, q{manta-}
                    . $parameter_href->{bioconda}{manta}
                    . $parameter_href->{bioconda_manta_patch}, q{bin}, $binary
                );
                my $link_path = catfile(
                    $parameter_href->{conda_prefix_path}, $binary
                );
                gnu_link(
                    {
                        FILEHANDLE  => $FILEHANDLE,
                        target_path => $target_path,
                        link_path   => $link_path,
                        symbolic    => 1,
                        force       => 1,
                    }
                );
                print $FILEHANDLE $NEWLINE;
            }
            print $FILEHANDLE $NEWLINE;

            ## Make file executable
            enable_executable(
                {
                    parameter_href => $parameter_href,
                    FILEHANDLE     => $FILEHANDLE,
                    binary         => q?configManta.py?,
                }
            );
        }
    }
    return;
}

sub perl {

##perl

##Function : Installs perl
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

    if ( $ENV{PATH} =~ /perl-$parameter_href->{perl_version}/ ) {

        if ( $parameter_href->{noupdate} ) {

            print STDERR 'Found perl-'
              . $parameter_href->{perl_version}
              . ' in your path', "\n";
            print STDERR 'Skipping writting installation for perl-'
              . $parameter_href->{perl_version}, "\n";
        }
        else {

            if ( $parameter_href->{perl_install} ) {

                ## Removing specific Perl version
                print $FILEHANDLE '### Removing specific perl version', "\n";
                gnu_rm(
                    {
                        infile_path => '$HOME/perl-'
                          . $parameter_href->{perl_version},
                        force      => 1,
                        recursive  => 1,
                        FILEHANDLE => $FILEHANDLE,
                    }
                );
                print $FILEHANDLE "\n\n";

                install_perl_cpnam(
                    {
                        parameter_href => $parameter_href,
                        FILEHANDLE     => $FILEHANDLE,
                    }
                );
            }

            perl_modules(
                {
                    parameter_href => $parameter_href,
                    FILEHANDLE     => $FILEHANDLE,
                }
            );
        }
    }
    else {

        if ( $parameter_href->{perl_install} ) {

            install_perl_cpnam(
                {
                    parameter_href => $parameter_href,
                    FILEHANDLE     => $FILEHANDLE,
                    path           => 1,
                }
            );
        }

        perl_modules(
            {
                parameter_href => $parameter_href,
                FILEHANDLE     => $FILEHANDLE,
            }
        );
    }
    return;
}

sub install_perl_cpnam {

##install_perl_cpnam

##Function : Install perl CPANM
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to
##         : $path           => Export path if provided {Optional}

    my ($arg_href) = @_;

    ## Default(s)
    my $path;

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
        path       => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$path
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $pwd = cwd();

    print STDERR 'Writting install instructions for perl and Cpanm', "\n";

    ## Install specific perl version
    print $FILEHANDLE '### Install specific perl version', "\n";

    ## Move to Home
    print $FILEHANDLE '## Move to $HOME', "\n";
    gnu_cd(
        {
            directory_path => q?$HOME?,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Download
    print $FILEHANDLE '## Download perl', "\n";
    Program::Download::Wget::wget(
        {
            url => 'http://www.cpan.org/src/5.0/perl-'
              . $parameter_href->{perl_version}
              . '.tar.gz',
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $parameter_href->{quiet},
            verbose      => $parameter_href->{verbose},
            outfile_path => 'perl-'
              . $parameter_href->{perl_version}
              . '.tar.gz',
        }
    );
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE '## Extract', "\n";
    print $FILEHANDLE "tar xzf perl-"
      . $parameter_href->{perl_version}
      . ".tar.gz";
    print $FILEHANDLE "\n\n";

    ## Move to perl directory
    print $FILEHANDLE '## Move to perl directory', "\n";
    gnu_cd(
        {
            directory_path => 'perl-' . $parameter_href->{perl_version},
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Configure
    print $FILEHANDLE '## Configure', "\n";
    print $FILEHANDLE './Configure -des -Dprefix=$HOME/perl-'
      . $parameter_href->{perl_version}, "\n";
    print $FILEHANDLE 'make', "\n";

    if ( !$parameter_href->{perl_skip_test} ) {

        print $FILEHANDLE 'make test', "\n";
    }
    print $FILEHANDLE 'make install', "\n\n";

    if ($path) {

        ## Export path
        print $FILEHANDLE '## Export path', "\n";
        print $FILEHANDLE q{echo 'export PATH=$HOME/perl-}
          . $parameter_href->{perl_version}
          . q{/:$PATH' >> ~/.bashrc};
        print $FILEHANDLE "\n\n";
        print $FILEHANDLE 'export PATH=$HOME/perl-'
          . $parameter_href->{perl_version}
          . '/:$PATH';    #Use newly installed perl
        print $FILEHANDLE "\n\n";
    }

    ## Remove tar file
    print $FILEHANDLE '## Remove tar file', "\n";
    gnu_cd( { FILEHANDLE => $FILEHANDLE, } );

    print $FILEHANDLE '&& ';

    gnu_rm(
        {
            infile_path => 'perl-'
              . $parameter_href->{perl_version}
              . '.tar.gz',
            FILEHANDLE => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Move to back
    print $FILEHANDLE '## Move to original working directory', "\n";
    gnu_cd(
        {
            directory_path => $pwd,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE q{echo 'eval `perl -I ~/perl-}
      . $parameter_href->{perl_version}
      . q{/lib/perl5/ -Mlocal::lib=~/perl-}
      . $parameter_href->{perl_version}
      . q{/`' >> ~/.bash_profile };    #Add at start-up
    print $FILEHANDLE "\n\n";
    print $FILEHANDLE
      q{echo 'export PERL_UNICODE=SAD' >> ~/.bash_profile };    #Add at start-up
    print $FILEHANDLE "\n\n";

    ## Install perl modules via cpanm
    print $FILEHANDLE '## Install cpanm', "\n";
    Program::Download::Wget::wget(
        {
            url          => 'http://cpanmin.us',
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $parameter_href->{quiet},
            verbose      => $parameter_href->{verbose},
            outfile_path => '-',
        }
    );
    print $FILEHANDLE q{ | perl - -l $HOME/perl-}
      . $parameter_href->{perl_version}
      . q{/bin App::cpanminus --local-lib=~/perl-}
      . $parameter_href->{perl_version}
      . q{/ local::lib };
    print $FILEHANDLE "\n\n";

    ## Use newly installed perl
    print $FILEHANDLE q{eval `perl -I ~/perl-}
      . $parameter_href->{perl_version}
      . q{/lib/perl5/ -Mlocal::lib=~/perl-}
      . $parameter_href->{perl_version} . q{/` };
    print $FILEHANDLE "\n\n";

    ## Use newly installed perl
    print $FILEHANDLE q{PERL5LIB=~/perl-}
      . $parameter_href->{perl_version}
      . q{/lib/perl5};
    print $FILEHANDLE "\n\n";

    return;
}

sub perl_modules {

##perl_modules

##Function : Install perl modules via cpanm
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

    ## Install perl modules via cpanm
    print $FILEHANDLE '## Install perl modules via cpanm', "\n";
    print $FILEHANDLE 'cpanm ';

    if ( $parameter_href->{perl_modules_force} ) {

        print $FILEHANDLE '--force ';
    }
    print $FILEHANDLE join( q{ }, @{ $parameter_href->{perl_modules} } ) . q{ };
    print $FILEHANDLE "\n\n";

    return;
}

sub pip_install {

##pip_install

##Function : Writes install instructions for pip packages
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

    print STDERR 'Writting install instructions for pip packages', "\n";

    ## Install PIP packages in conda environment
    if ( exists( $parameter_href->{conda_environment} )
        && ( $parameter_href->{conda_environment} ) )
    {

        print $FILEHANDLE '### Install PIP packages in conda environment: '
          . $parameter_href->{conda_environment}, "\n";
    }
    else {

        print $FILEHANDLE '### Install PIP packages in conda main environment',
          "\n";
    }
   
    ## Only activate conda environment if supplied by user 
    if ( exists $parameter_href->{conda_environment} ) {
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

    ## Install PIP packages
    print $FILEHANDLE '## Install PIP packages', "\n";
    print $FILEHANDLE 'pip install ';

    if ( $parameter_href->{quiet} ) {

        print $FILEHANDLE '--quiet ';    #Do not display progress bar
    }

    ## Install all PIP packages
    foreach my $program ( keys %{ $parameter_href->{pip} } ) {

        print $FILEHANDLE $program . '=='
          . $parameter_href->{pip}{$program} . q{ };
    }
    print $FILEHANDLE "\n\n";
    
    ## Deactivate conda environment if conda_environment exists
    if ( exists $parameter_href->{conda_environment} ) {
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
    Program::Download::Wget::wget(
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
        $parameter_href->{conda_prefix_path}, q{share}, q{picard-tools-} 
        . $parameter_href->{picardtools}, q{picard.jar}
    );
    my $link_path = catfile(
        $parameter_href->{conda_prefix_path}, q{picard.jar}
    );
    gnu_link(
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
    Program::Download::Wget::wget(
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
            outfile_path => catdir( 
                $parameter_href->{conda_prefix_path}, 'bin' 
            ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Specifying target and link paths 
    my $target_path = q{sambamba_v} . $parameter_href->{bioconda}{sambamba};
    my $link_path = catfile(
        $parameter_href->{conda_prefix_path}, q{sambamba}
    );
    gnu_link(
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
    Program::Download::Wget::wget(
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
    Program::Download::Wget::wget(
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
    Program::Download::Wget::wget(
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
            $parameter_href->{conda_prefix_path}, q{share},  q{snpEff.} 
            . $parameter_href->{snpeff}, $binary
        );
        my $link_path = catfile(
            $parameter_href->{conda_prefix_path}, $binary
        );
        gnu_link(
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
      my $genome_version ( @{ $parameter_href->{snpeff_genome_versions} } ) {
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

sub varianteffectpredictor {

##varianteffectpredictor

##Function : Install varianteffectpredictor
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

    my $miniconda_bin_dir = catdir( $parameter_href->{conda_prefix_path},
        'ensembl-tools-release-' . $parameter_href->{varianteffectpredictor} );

    if ( -d $miniconda_bin_dir ) {

        print STDERR 'Found varianteffectpredictor in miniconda directory: '
          . $miniconda_bin_dir, "\n";

        if ( $parameter_href->{noupdate} ) {

            print STDERR
'Skipping writting installation process for varianteffectpredictor',
              "\n";
            return;
        }
        else {

            ## Removing varianteffectpredictor
            print $FILEHANDLE '### Removing varianteffectpredictor', "\n";
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

        print STDERR 'Writting install instructions for varianteffectpredictor',
          "\n";
    }

    ## Install VEP
    print $FILEHANDLE '### Install varianteffectpredictor', "\n";
    
    ## Only activate conda environment if supplied by user
    if ( exists $parameter_href->{conda_environment} ) {
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

    ##Make sure that the cache directory exists
    gnu_mkdir(
        {
            indirectory_path => $parameter_href->{vep_cache_dir},
            parents          => 1,
            FILEHANDLE       => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Move to miniconda environment
    gnu_cd(
        {
            directory_path => catdir( $parameter_href->{conda_prefix_path} ),
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Download
    print $FILEHANDLE '## Download VEP', "\n";
    Program::Download::Wget::wget(
        {
            url => 'https://github.com/Ensembl/ensembl-tools/archive/release/'
              . $parameter_href->{varianteffectpredictor} . '.zip',
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $parameter_href->{quiet},
            verbose      => $parameter_href->{verbose},
            outfile_path => 'VariantEffectPredictor-'
              . $parameter_href->{varianteffectpredictor} . '.zip',
        }
    );
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE '## Extract', "\n";
    print $FILEHANDLE 'unzip VariantEffectPredictor-'
      . $parameter_href->{varianteffectpredictor} . '.zip';
    print $FILEHANDLE "\n\n";

    ## Move to VariantEffectPredictor directory
    print $FILEHANDLE '## Move to VariantEffectPredictor directory', "\n";

    ## Needed since VEP only uses major version to name the zip file
    my $vep_install_version = $parameter_href->{varianteffectpredictor};
    if ( $vep_install_version =~ /(\d+)/ ) {

        $vep_install_version = $1;
    }
    gnu_cd(
        {
            directory_path => catdir(
                'ensembl-tools-release-' . $vep_install_version, 'scripts',
                'variant_effect_predictor'
            ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Install VEP
    print $FILEHANDLE '## Install VEP', "\n";
    print $FILEHANDLE 'perl INSTALL.pl ';
    print $FILEHANDLE '--AUTO '
      . $parameter_href->{vep_auto_flag}
      ;    #a (API), l (FAIDX/htslib), c (cache), f (FASTA)

    if (   ( exists( $parameter_href->{vep_plugins} ) )
        && ( @{ $parameter_href->{vep_plugins} } ) )
    {

        print $FILEHANDLE 'p ';    #p (plugins)
        print $FILEHANDLE '-g '
          . join( ',', @{ $parameter_href->{vep_plugins} } )
          . q{ };                  #Plugins in comma sep string
    }
    else {

        print $FILEHANDLE q{ };    #Add whitespace if no plugins
    }

    print $FILEHANDLE '-c '
      . $parameter_href->{vep_cache_dir}
      . q{ };                      #Cache directory
    print $FILEHANDLE '-s homo_sapiens ';

    ## Only install first assembly version since VEP install cannot handle multiple versions at the same time
    print $FILEHANDLE '--ASSEMBLY '
      . $parameter_href->{vep_assemblies}[0] . q{ };
    print $FILEHANDLE "\n\n";

    if (   ( scalar( @{ $parameter_href->{vep_assemblies} } ) > 1 )
        && ( $parameter_href->{vep_auto_flag} =~ /c|f/ ) )
    {

        for (
            my $assembly_version = 1 ;
            $assembly_version <
            scalar( @{ $parameter_href->{vep_assemblies} } ) ;
            $assembly_version++
          )
        {

            print $FILEHANDLE
              '## Install additional VEP cache assembly version', "\n";
            print $FILEHANDLE 'perl INSTALL.pl ';
            print $FILEHANDLE '--AUTO cf '
              ;    #a (API), l (FAIDX/htslib), c (cache), f (FASTA), p (plugins)
            print $FILEHANDLE '-c '
              . $parameter_href->{vep_cache_dir}
              . q{ };    #Cache directory
            print $FILEHANDLE '-s homo_sapiens ';

            ## Only install first assembly version since VEP install cannot handle multiple versions at the same time
            print $FILEHANDLE '--ASSEMBLY '
              . $parameter_href->{vep_assemblies}[$assembly_version] . q{ };
            print $FILEHANDLE "\n\n";
        }
    }

    if ( exists( $parameter_href->{vep_plugins} ) ) {

        if ( grep { $_ eq 'LoFtool' } @{ $parameter_href->{vep_plugins} } ) {

            ##Add LofTool required text file
            print $FILEHANDLE '##Add LofTool required text file', "\n";
            Program::Download::Wget::wget(
                {
                    url =>
'https://raw.githubusercontent.com/Ensembl/VEP_plugins/master/LoFtool_scores.txt',
                    FILEHANDLE   => $FILEHANDLE,
                    quiet        => $parameter_href->{quiet},
                    verbose      => $parameter_href->{verbose},
                    outfile_path => q{$HOME/.vep/Plugins/LoFtool_scores.txt},
                }
            );
            print $FILEHANDLE "\n\n";
        }

        if ( grep { $_ eq 'Lof' } @{ $parameter_href->{vep_plugins} } ) {

            ## Add Lof required perl splice script
            print $FILEHANDLE '##Add Lof required perl splice script', "\n";
            Program::Download::Wget::wget(
                {
                    url =>
'https://raw.githubusercontent.com/konradjk/loftee/master/splice_module.pl',
                    FILEHANDLE   => $FILEHANDLE,
                    quiet        => $parameter_href->{quiet},
                    verbose      => $parameter_href->{verbose},
                    outfile_path => q{$HOME/.vep/Plugins/splice_module.pl},
                }
            );
            print $FILEHANDLE "\n\n";

            ## Add Lof optional human_ancestor_fa
            print $FILEHANDLE '##Add Lof optional human_ancestor_fa', "\n";
            Program::Download::Wget::wget(
                {
                    url =>
'https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz',
                    FILEHANDLE   => $FILEHANDLE,
                    quiet        => $parameter_href->{quiet},
                    verbose      => $parameter_href->{verbose},
                    outfile_path => catfile(
                        $parameter_href->{vep_cache_dir},
                        'human_ancestor.fa.gz'
                    ),
                }
            );
            print $FILEHANDLE "\n\n";

            ## Uncompress
            print $FILEHANDLE 'bgzip -d '
              . catfile( $parameter_href->{vep_cache_dir},
                'human_ancestor.fa.gz' )
              . q{ };
            print $FILEHANDLE "\n\n";

            ## Add Lof optional human_ancestor_fa index
            Program::Download::Wget::wget(
                {
                    url =>
'https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.fai',
                    FILEHANDLE   => $FILEHANDLE,
                    quiet        => $parameter_href->{quiet},
                    verbose      => $parameter_href->{verbose},
                    outfile_path => catfile(
                        $parameter_href->{vep_cache_dir},
                        'human_ancestor.fa.fai'
                    ),
                }
            );
            print $FILEHANDLE "\n\n";
        }
    }

    ## Clean up
    print $FILEHANDLE '## Clean up', "\n";
    gnu_rm(
        {
            infile_path => catdir(
                $parameter_href->{conda_prefix_path},
                'VariantEffectPredictor-'
                  . $parameter_href->{varianteffectpredictor} . '.zip'
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
    if ( exists $parameter_href->{conda_environment} ) {
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

        print STDERR 'Found Root in miniconda directory: ' . $miniconda_bin_dir,
          "\n";

        if ( $parameter_href->{noupdate} ) {

            print STDERR 'Skipping writting installation process for Root',
              "\n";
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

        print STDERR 'Writting install instructions for Root', "\n";
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
    Program::Download::Wget::wget(
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
    if ( exists $parameter_href->{conda_environment} ) {
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
    Program::Download::Wget::wget(
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
    if ( exists $parameter_href->{conda_environment} ) {
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
    if ( exists $parameter_href->{conda_environment} ) {
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
    Program::Download::Wget::wget(
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
        $parameter_href->{conda_prefix_path}, q{TIDDIT-} 
        . $parameter_href->{tiddit}, qw{ bin TIDDIT }
    );
    my $link_path = catfile(
        $parameter_href->{conda_prefix_path}, qw{ bin TIDDIT }
    );
    gnu_link(
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
    if ( exists $parameter_href->{conda_environment} ) {
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
    if ( exists $parameter_href->{conda_environment} ) {
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
    Program::Download::Wget::wget(
        {
            url => 'https://github.com/J35P312/SVDB/archive/'
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
            directory_path => 'SVDB-' . $parameter_href->{svdb},
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    ## Install
    print $FILEHANDLE '## Install', "\n";
    print $FILEHANDLE 'pip install . ';
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
    if ( exists $parameter_href->{conda_environment} ) {
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
    if ( exists $parameter_href->{conda_environment} ) {
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
    Program::Download::Wget::wget(
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
    print $FILEHANDLE '## Configure',                    "\n";
    print $FILEHANDLE 'pip install numpy Cython',        "\n";
    print $FILEHANDLE 'pip install -r requirements.txt', "\n";
    print $FILEHANDLE 'pip install -e .';
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
    if ( exists $parameter_href->{conda_environment} ) {
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

            print STDERR 'Found '
              . $program_name
              . ' version '
              . $program_version
              . ' in miniconda directory: '
              . catdir( $parameter_href->{conda_prefix_path}, 'bin' ), "\n";

            if ( $parameter_href->{noupdate} ) {

                print STDERR 'Skipping writting installation process for '
                  . $program_name . q{ }
                  . $program_version, "\n";
                return 1;
            }
            print STDERR 'Writting install instructions for ' . $program_name,
              "\n";
        }
        else {

            print STDERR 'Found '
              . $program_name
              . ' in miniconda directory: '
              . catdir( $parameter_href->{conda_prefix_path}, 'bin' ), "\n";

            if ( $parameter_href->{noupdate} ) {

                print STDERR 'Skipping writting installation process for '
                  . $program_name, "\n";
                return 1;
            }
            print STDERR 'Writting install instructions for ' . $program_name,
              "\n";
        }
        return 0;
    }
    else {

        print STDERR 'Writting install instructions for ' . $program_name, "\n";
        return 0;
    }
    return;
}

sub enable_executable {

##enable_executable

##Function : Make file executable
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE, $binary
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => FILEHANDLE to write to
##         : $binary         => The binary file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;
    my $binary;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
        binary =>
          { required => 1, defined => 1, strict_type => 1, store => \$binary },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Enable executable
    print $FILEHANDLE '## Enable executable', "\n";
    gnu_cd(
        {
            directory_path =>
              catdir( $parameter_href->{conda_prefix_path}, 'bin' ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n";

    print $FILEHANDLE 'chmod a+x ';
    print $FILEHANDLE $binary . q{ };
    print $FILEHANDLE "\n\n";

    ## Move to back
    print $FILEHANDLE '## Move to original working directory', "\n";
    gnu_cd(
        {
            directory_path => $pwd,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    return;
}

sub check_mt_codon_table {

##check_mt_codon_table

##Function : Check and if required add the vertebrate mitochondrial codon table to snpeff config
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE, $share_dir, $config_file, $genome_version_ref
##         : $parameter_href     => Holds all parameters
##         : $FILEHANDLE         => FILEHANDLE to write to
##         : $share_dir          => The conda env shared directory
##         : $config_file        => The config config_file
##         : $genome_version_ref => snpeff genome version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;
    my $share_dir;
    my $config_file;
    my $genome_version_ref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
        share_dir  => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$share_dir
        },
        config_file => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$config_file
        },
        genome_version_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$genome_version_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $pwd = cwd();
    my $detect_regexp =
        q?perl -nae 'if($_=~/?
      . $$genome_version_ref
      . q?.MT.codonTable/) {print 1}' ?;
    my $add_regexp =
        q?perl -nae 'if($_=~/?
      . $$genome_version_ref
      . q?.reference/) {print $_; print "?
      . $$genome_version_ref
      . q?.MT.codonTable : Vertebrate_Mitochondrial\n"} else {print $_;}' ?;
    my $ret;

    if ( -f catfile( $share_dir, $config_file ) ) {

        $ret = `$detect_regexp $share_dir/$config_file`;
    }
    if ( !$ret ) {    #No MT.codonTable in config

        print $FILEHANDLE '## Adding '
          . $$genome_version_ref
          . '.MT.codonTable : Vertebrate_Mitochondrial to '
          . $share_dir
          . $config_file, "\n";

        ## Add MT.codon Table to config
        print $FILEHANDLE $add_regexp . q{ }
          . catfile( $share_dir, $config_file ) . ' > '
          . catfile( $share_dir, $config_file . '.tmp' ), "\n";
        gnu_mv(
            {
                infile_path  => catfile( $share_dir, $config_file . '.tmp' ),
                outfile_path => catfile( $share_dir, $config_file ),
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        print $FILEHANDLE "\n\n";

    }
    else {

        print STDERR 'Found MT.codonTable in '
          . catfile( $share_dir, 'snpEff.config' )
          . '. Skipping addition to snpEff config', "\n";
    }
    return;
}

sub snpeff_download {

##snpeff_download

##Function : Write instructions to download snpeff database. This is done by install script to avoid race conditin when doing first analysis run in MIP
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE, $genome_version_ref
##         : $parameter_href     => Holds all parameters
##         : $FILEHANDLE         => FILEHANDLE to write to
##         : $genome_version_ref => snpeff genome version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;
    my $genome_version_ref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
        genome_version_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$genome_version_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    ## Only activate conda environment if supplied by user
    if ( exists $parameter_href->{conda_environment} ) {
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

    print $FILEHANDLE 'java -Xmx2g ';
    print $FILEHANDLE '-jar '
      . catfile( $parameter_href->{conda_prefix_path}, 'bin', 'snpEff.jar' )
      . q{ };
    print $FILEHANDLE 'download ';
    print $FILEHANDLE ' -v ';
    print $FILEHANDLE $$genome_version_ref . q{ };
    print $FILEHANDLE '-c '
      . catfile( $parameter_href->{conda_prefix_path}, 'bin', 'snpEff.config' )
      . q{ };
    print $FILEHANDLE "\n\n";

    ## Deactivate conda environment if conda_environment exists
    if ( exists $parameter_href->{conda_environment} ) {
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
    if ( exists $parameter_href->{conda_environment} ) {
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

    print STDERR 'Writting install instructions for references', "\n";

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
    if ( exists $parameter_href->{conda_environment} ) {
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
