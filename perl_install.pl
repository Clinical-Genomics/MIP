#!/usr/bin/env perl

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use 5.010000;    # Require perl v5.10.0 or higher
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };

use Getopt::Long;
use Cwd;
use IO::Handle;
use File::Basename qw{ basename fileparse };
use File::Spec::Functions qw{ catfile catdir rel2abs };
use ExtUtils::Installed;
use FindBin qw{ $Bin };

## Local modules
use lib catdir( $Bin, q{lib} );
use MIP::Script::Utils_v5_10 qw{ help set_default_array_parameters };
use MIP::Program::Download::Wget_v5_10 qw{ wget };
use MIP::Package_manager::Cpanm_v5_10 qw{ cpanm_install_module };
use MIP::Check::Modules_v5_10 qw{ check_perl_modules };

our $USAGE = build_usage( {} );

### Set parameter default
my %parameter;
my %array_parameter;

## Perl defaults
$parameter{perl_install_dir} = $ENV{HOME};
$parameter{perl_version}     = q{5.26.0};

## Module defaults
$array_parameter{perl_modules} = [
    q{Modern::Perl},              # MIP
    q{IPC::System::Simple},       # MIP
    q{Path::Iterator::Rule},      # MIP
    q{YAML},                      # MIP
    q{Log::Log4perl},             # MIP
    q{List::Util},                # MIP
    q{List::MoreUtils},           # MIP
    q{Readonly},                  # MIP
    q{Try::Tiny},                 # MIP
    q{Array::Utils},              # MIP
    q{IO::Uncompress::Gunzip},    # VEP
    q{HTML::Lint},                # VEP
    q{Archive::Zip},              # VEP
    q{Archive::Extract},          # VEP
    q{DBI},                       # VEP
    q{JSON},                      # VEP
    q{DBD::mysql},                # VEP
    q{CGI},                       # VEP
    q{Sereal::Encoder},           # VEP
    q{Sereal::Decoder},           # VEP
    q{Bio::Root::Version},        # VEP
    q{Module::Build},             # VEP
    q{File::Copy::Recursive},     # VEP
];

our $VERSION = q{1.0.2};

###User Options
GetOptions(
    q{pid|perl_install_dir=s}       => \$parameter{perl_install_dir},
    q{pev|perl_version=s}           => \$parameter{perl_version},
    q{psi|perl_skip_install}        => \$parameter{perl_skip_install},
    q{pfi|perl_force_install}       => \$parameter{perl_force_install},
    q{pst|perl_skip_test}           => \$parameter{perl_skip_test},
    q{pm|perl_modules:s}            => \@{ $parameter{perl_modules} },
    q{pmf|perl_modules_force}       => \$parameter{perl_modules_force},
    q{pma|perl_modules_append:s}    => \@{ $parameter{perl_modules_append} },
    q{ppd|print_parameters_default} => sub {

        # Display parameter defaults
        print_parameters(
            {
                parameter_href       => \%parameter,
                array_parameter_href => \%array_parameter
            }
        );
        exit;
    },
    q{q|quiet} => \$parameter{quiet},
    q{h|help}  => sub {

        # Display help text
        say STDOUT $USAGE;
        exit;
    },
    q{ver|version} => sub {

        # Display version number
        say STDOUT qq{\n} . basename($PROGRAM_NAME) . q{ } . $VERSION, qq{\n};
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

## Declaring hash references
my $parameter_href;
my $array_parameter_href;

## Set active parameters
set_default_array_parameters(
    {
        parameter_href       => \%parameter,
        array_parameter_href => \%array_parameter
    }
);

## Convert relative path to absolute
$parameter{perl_install_dir} = rel2abs( $parameter{perl_install_dir} );

##########
###MAIN###
##########

## Create install file
# Create anonymous filehandle
my $FILEHANDLE = IO::Handle->new();

# Installation instruction file
my $file_name_path = catfile( cwd(), q{perl_install.sh} );

open $FILEHANDLE, q{>}, $file_name_path
  or
  croak( q{Cannot write to '} . $file_name_path . q{' :} . $OS_ERROR . qq{\n} );

say STDERR q{## Writing perl install instructions to: } . $file_name_path;

# Add shebang
say {$FILEHANDLE} q{#!/usr/bin/env/bash};
say {$FILEHANDLE} qq{\n};

## Write perl install instructions
# Install Perl if upgrade needed
if ( $parameter{perl_skip_install} ) {
    say STDERR q{## Skipping perl installation};
    say {$FILEHANDLE} q{## Skipping installation of perl};
}
elsif ( $PERL_VERSION ge q{v} . $parameter{perl_version}
    and not $parameter{perl_force_install} )
{
    say STDERR q{## Perl version requirements already met, }
      . q{current perl version: }
      . $PERL_VERSION;
    say STDERR q{## Re-run the script with }
      . q{'--perl_force_install' flag to force perl installation};

    say {$FILEHANDLE} q{## Perl version requirement met};

    # Change the perl version string if a newer verseion is found
    $parameter{perl_version} = substr $PERL_VERSION, 1;
}
else {
    say STDERR q{## Writing perl installation recipe};
    say {$FILEHANDLE} q{## Perl installation recipe};

    install_perl(
        {
            parameter_href => \%parameter,
            FILEHANDLE     => $FILEHANDLE
        }
    );
}
say {$FILEHANDLE} qq{\n};

## Write cpanm install instructions
# Check if cpanm alread in path
my ($installed) = ExtUtils::Installed->new();
my (@modules)   = $installed->modules();

if (    grep { m/ cpanm /xms } @modules
    and $PERL_VERSION ge q{v} . $parameter{perl_version}
    and not $parameter{perl_force_install} )
{
    say STDERR q{## Cpanm already installed};
    say {$FILEHANDLE} q{## Cpanm already installed} . qq{\n};
}
else {
    say STDERR q{## Writing Cpanm install recipe};
    say {$FILEHANDLE} q{## Cpanm installation recipe};

    install_cpanm(
        {
            parameter_href => \%parameter,
            FILEHANDLE     => $FILEHANDLE
        }
    );
}
say {$FILEHANDLE} qq{\n};

## Install MIP requred modules via cpanm
install_cpanm_modules(
    {
        parameter_href => \%parameter,
        FILEHANDLE     => $FILEHANDLE
    }
);

close $FILEHANDLE
  or croak q{Couldn't close } . $file_name_path . q{: } . $OS_ERROR;

#################
###SubRoutines###
#################

sub install_perl {

## install_perl

## Function  : Write perl isntallation recipe to FILEHANDLE
## Returns   : ""
## Arguments : $parameter_href, $array_parameter_href
##           : $parameter_href       => Parameters hash {REF}
##           : $array_parameter_href => Hold the array parameter defaults as {REF}

    my ($arg_href) = @_;

    # Flatten arguments
    $parameter_href = $arg_href->{parameter_href};
    $FILEHANDLE     = $arg_href->{FILEHANDLE};

    my $perl_install_dir = $arg_href->{parameter_href}{perl_install_dir};
    my $perl_version     = $arg_href->{parameter_href}{perl_version};
    my $perl_skip_test   = $arg_href->{parameter_href}{perl_skip_test};

    my $pwd = cwd();

    ## Check installation path
    if ( not -d $perl_install_dir ) {

        # Create installation dir
        say {$FILEHANDLE} q{# Create perl installation dir};
        say {$FILEHANDLE} q{mkdir -p } . $perl_install_dir . qq{\n};
    }

    ## Writing wget command
    say {$FILEHANDLE} q{# Download perl from cpan};

    ## Download source
    my $url = q{http://www.cpan.org/src/5.0/perl-} . $perl_version . q{.tar.gz};

    ## Path to downloaded file
    my $outfile_path =
      catfile( $perl_install_dir, q{perl-} . $perl_version . q{.tar.gz} );

    my @commands = wget(
        {
            url          => $url,
            outfile_path => $outfile_path,
            quiet        => $parameter_href->{quiet},
            verbose      => $parameter_href->{verbose}
        }
    );

    ## Write commands to $FILEHANDLE
    print {$FILEHANDLE} join( q{ }, @commands ) . q{ };

    print {$FILEHANDLE} qq{\n};

    say {$FILEHANDLE} q{cd } . $perl_install_dir . qq{\n};

    ## Writing unpack command
    say {$FILEHANDLE} q{# Unpack and remove tar file};

    my $tar_file =
      catfile( $perl_install_dir, q{perl-} . $perl_version . q{.tar.gz} );

    say {$FILEHANDLE} q{tar -xzf }
      . $tar_file
      . q{ && rm }
      . $tar_file . qq{\n};

    ## Building perl
    say {$FILEHANDLE} q{# Build perl};

    my $perl_install_path =
      catdir( $perl_install_dir, q{perl-} . $perl_version );

    say {$FILEHANDLE} q{cd } . $perl_install_path;
    say {$FILEHANDLE} q{./Configure -des -Dprefix=} . $perl_install_path;
    say {$FILEHANDLE} q{make};

    if ( not $perl_skip_test ) {
        say {$FILEHANDLE} q{make test};
    }

    say {$FILEHANDLE} q{make install};
    say {$FILEHANDLE} q{cd } . $pwd . qq{\n};

    ## Editing PATH and .bashrc
    say {$FILEHANDLE}
      q{echo "# Added by MIP's perl installer $(date)" >> ~/.bashrc};
    say {$FILEHANDLE} q{# Edit PATH and .bashrc};
    say {$FILEHANDLE} q{echo 'export PATH=} . $perl_install_path
      . q{/:$PATH' >> ~/.bashrc};

    # Use newly installed perl
    say {$FILEHANDLE} q{export PATH=}
      . $perl_install_path
      . q{/:$PATH} . qq{\n};

    ## Editing bash_profile
    say {$FILEHANDLE} q{# Edit bash_profile};

    # Add at start-up
    say {$FILEHANDLE}
      q{echo "# Added by MIP's perl installer $(date)" >> ~/.bash_profile};
    say {$FILEHANDLE} q{echo 'eval `perl -I }
      . catdir( $perl_install_path, qw{ lib perl5 } )
      . q{ -Mlocal::lib=}
      . $perl_install_path
      . q{/`' >> ~/.bash_profile };

    # Handle Unicode
    say {$FILEHANDLE} q{echo 'export PERL_UNICODE=SAD' >> ~/.bash_profile }
      . qq{\n};

    return;
}

sub install_cpanm {

## install_cpanm

## Function  : Writes recipe for installing cpanm and setting up a local library for perl modules
## Returns   : ""
## Arguments : $parameter_href, $FILEHANDLE
##           : $parameter_href => Hash with paramters {REF}
##           : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    $parameter_href = $arg_href->{parameter_href};
    $FILEHANDLE     = $arg_href->{FILEHANDLE};

    ## Writing wget command and installing cpanm
    say {$FILEHANDLE} q{# Download Cpanm from cpan};

    my $url          = q{http://cpanmin.us};
    my $outfile_path = q{-};

    # Download wget and print to STDOUT
    my @commands = wget(
        {
            url          => $url,
            outfile_path => $outfile_path,
            quiet        => $parameter_href->{quiet},
            verbose      => $parameter_href->{verbose}
        }
    );

    ## Write commands to $FILEHANDLE
    print {$FILEHANDLE} join( q{ }, @commands ) . q{ };

    my $perl_install_path = catdir(
        $parameter_href->{perl_install_dir},
        q{perl-} . $parameter_href->{perl_version}
    );

    # Pipes the STDOUT output from wget to perl and run cpanm
    say {$FILEHANDLE} q{| perl - }

      # Tell cpanm where to install modules
      . q{-l } . $perl_install_path . q{/bin }

      # Install modules
      . q{App::cpanminus --local-lib=}
      . $perl_install_path
      . q{/ local::lib } . qq{\n};

    ## Setting up local cpan library
    say {$FILEHANDLE} q{# Set up local cpan library};

    ## Setting up environment variables to use local modules
    say {$FILEHANDLE} q{eval `perl -I }
      . catdir( $perl_install_path, qw{ lib perl5 } )
      . q{ -Mlocal::lib=}
      . $perl_install_path . q{/` } . qq{\n};

    say {$FILEHANDLE} q{echo 'eval `perl -I }
      . catdir( $perl_install_path, qw{ lib perl5 } )
      . q{ -Mlocal::lib=}
      . $perl_install_path
      . q{/`' >> ~/.bash_profile } . qq{\n};

    ## Add to .bashrc to enable man on local modules
    say {$FILEHANDLE} q{echo 'export MANPATH=}
      . catdir( $perl_install_path, q{man} )
      . q{:$MANPATH' >> ~/.bashrc} . qq{\n};

    ## Use newly installed perl
    say {$FILEHANDLE} q{PERL5LIB=}
      . catdir( $perl_install_path, qw{ lib perl5} ) . qq{\n};

    ## source .bashrc and .bash_profile
    say {$FILEHANDLE} q{echo 'source ~/.bash_profile'} . qq{\n};
    say {$FILEHANDLE} q{echo 'source ~/.bashrc'} . qq{\n};

    return;
}

sub install_cpanm_modules {

## install_cpanm_modules

## Function  : Perl wrapper for writing cpanm recipe to $FILEHANDLE.
## Returns   : ""
## Arguments : $parameter_href, $FILEHANDLE
##           : $parameter_href => Hash with paramters {REF}
##           : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    $parameter_href = $arg_href->{parameter_href};
    $FILEHANDLE     = $arg_href->{FILEHANDLE};

    my @perl_modules = @{ $arg_href->{parameter_href}{perl_modules} };
    my @perl_modules_append =
      @{ $arg_href->{parameter_href}{perl_modules_append} };

    ## Modules to install
    my %perl_modules;
    if (@perl_modules_append) {
        @perl_modules_append = split( /,/, join( q{,}, @perl_modules_append ) );

        @perl_modules = ( @perl_modules, @perl_modules_append );

        # Remove any duplicate modules
        %perl_modules = map { $_ => 1 } @perl_modules;
        @perl_modules = keys %perl_modules;
    }

    ## Check aginst what's already installed
    my @modules_to_install =
      check_perl_modules( { modules_ref => \@perl_modules } );

    if (@modules_to_install) {
        say STDERR q{## Writing recipe for installation of Cpanm modules};
        say {$FILEHANDLE} q{## Install modules required by MIP via Cpanm};
        my @commands = cpanm_install_module(
            {
                modules_ref => \@modules_to_install,
                force       => $parameter_href->{force},
                quiet       => $parameter_href->{quiet},
                verbose     => $parameter_href->{verbose}
            }
        );
        say {$FILEHANDLE} join q{ }, @commands;
    }
    else {
        say STDERR q{## Required Cpanm modules already installed};
        say {$FILEHANDLE} q{## Required Cpanm modules already installed};
    }

    return;
}

sub build_usage {

## build_usage

## Function : Build the USAGE instructions
## Returns  : ""
## Arguments: $script_name

    my ($arg_href) = @_;

    my $script_name = $arg_href->{script_name};

    $script_name = basename($PROGRAM_NAME);

    return <<"END_USAGE";
    $script_name [options]

    ## Perl installation
    -pid/--perl_install_dir     Set perl installation directory (default: home folder) 
    -pev/--perl_version         Set perl version (default: 5.18.2)
    -psi/--perl_skip_install    Skip perl installation (supply flag to enable)
    -pfi/--perl_force_install   Force perl installation (supply flag to enable)
    -pst/--perl_skip_test       Skip perl installation tests (supply flag to enable)


    ## Perl modules
    -pm/--perl_modules          Set the perl modules to be installed via cpanm (default: [ 
                                "Modern::Perl",         "List::Util",      "IPC::System::Simple", 
                                "Path::Iterator::Rule", "YAML",             "Log::Log4perl", 
                                "Set::IntervalTree",    "Net::SSLeay",      "LWP::Simple", 
                                "LWP::Protocol::https", "Archive::Zip",     "Archive::Extract", 
                                "DBI",                  "JSON",             "DBD::mysql", 
                                "CGI",                  "Sereal::Encoder",  "Sereal::Decoder", 
                                "Bio::Root::Version",   "Module::Build",    "Readonly"
                                ] )
    -pmf/--perl_modules_force   Force installation of perl modules (supply flag to enable)
    -pma/--perl_modules_append  Install additional perl modules aside from default (supply csv list)
   

    ## Utility
    -ppd/--print_parameters_default     Print the parameter defaults (supply flag to enabble)
    -q/--quiet                          No output from individual program that has a quiet flag)
    -h/--help                           Display this help message (supply flag to enable)
    -ver/--version                      Display version and exit (supply flag to enable)
    -v/--verbose                        Verbose output (supply flag to enable)
END_USAGE
}

sub print_parameters {

## print_parameters

## Function  : Print all parameters and the default values
## Returns   : ""
## Arguments : $parameter_href, $array_parameter_href
##           : $parameter_href => Holds all parameters {REF}
##           : $array_parameter_href => Hold the array parameter defaults as {REF}

    my ($arg_href) = @_;

    $parameter_href       = $arg_href->{parameter_href};
    $array_parameter_href = $arg_href->{array_parameter_href};

    ## Set default for array parameters
    set_default_array_parameters(
        {
            parameter_href       => $parameter_href,
            array_parameter_href => $array_parameter_href
        }
    );

  KEY:
    foreach my $key ( keys %{$parameter_href} ) {
        if ( ref( $parameter_href->{$key} ) !~ m/ ARRAY | HASH /xms ) {

            print STDOUT $key . q{ };
            if ( $parameter_href->{$key} ) {
                say STDOUT $parameter_href->{$key};
            }
            else {    ##Boolean value

                say STDOUT q{0};
            }
        }
        elsif ( ref( $parameter_href->{$key} ) =~ m/ HASH /xms ) {

          PROGRAM:
            foreach my $program ( keys %{ $parameter_href->{$key} } ) {

                say STDOUT $key . q{ } . $program . q{: }
                  . $parameter_href->{$key}{$program};
            }
        }
        elsif ( ref( $parameter_href->{$key} ) =~ m/ ARRAY /xms ) {

            say STDOUT $key . q{: } . join q{ }, @{ $parameter_href->{$key} };
        }
    }
    return;
}
