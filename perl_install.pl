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

our $USAGE = build_usage( {} );

### Set parameter default
my %parameter;
my %array_parameter;

## Perl defaults
$parameter{perl_install_dir} = $ENV{HOME};
$parameter{perl_version}     = q{5.18.2};

## Module defaults
$array_parameter{perl_modules} = [
    q{Modern::Perl},              # MIP
    q{IPC::System::Simple},       # MIP
    q{Path::Iterator::Rule},      # MIP
    q{YAML},                      # MIP
    q{Log::Log4perl},             # MIP
    q{List::Util},                # MIP
    q{List::MoreUtils},           # MIP
    q{Readonly},                  # VEP
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

our $VERSION = q{1.0.0};

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
        print_parameters(
            {
                parameter_href       => \%parameter,
                array_parameter_href => \%array_parameter
            }
        );
        exit;
    },    # Display parameter defaults
    q{q|quiet} => \$parameter{quiet},
    q{h|help}  => sub {
        say STDOUT $USAGE;
        exit;
    },    #Display help text
    q{ver|version} => sub {
        say STDOUT qq{\n} . basename($PROGRAM_NAME) . q{ } . $VERSION, qq{\n};
        exit;
    },    #Display version number
    q{v|verbose} => \$parameter{verbose},
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

## Write perl install instructions

## Install Perl if upgrade needed
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
    say STDERR
      q{## Re-run the script with '-pfi' flag to force perl installation};

    say {$FILEHANDLE} q{## Perl version requirement met};
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
my ($inst)    = ExtUtils::Installed->new();
my (@modules) = $inst->modules();

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
say {$FILEHANDLE} q{## Install modules reuired by MIP via cpanm};
install_modules(
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

sub set_default_array_parameters {

##set_default_array_parameters

## Function  : Set default for array parameters unless parameter already exists in parameter hash
## Returns   : ""
## Arguments : $parameter_href, $array_parameter_href
##           : $parameter_href       => Parameters hash {REF}
##           : $array_parameter_href => Hold the array parameter defaults as {REF}

    my ($arg_href) = @_;

    $parameter_href       = $arg_href->{parameter_href};
    $array_parameter_href = $arg_href->{array_parameter_href};
    $FILEHANDLE           = $arg_href->{FILEHANDLE};

  PARAMETER_NAME:
    foreach my $parameter_name ( keys %{$array_parameter_href} ) {
        if ( not @{ $parameter_href->{$parameter_name} } ) {
            $parameter_href->{$parameter_name} =
              $array_parameter_href->{$parameter_name};
        }
    }

    return;
}

sub install_perl {

## install_perl

## Function  : Write perl isntallation recipe to FILEHANDLE
## Returns   : ""
## Arguments : $parameter_href, $array_parameter_href
##           : $parameter_href       => Parameters hash {REF}
##           : $array_parameter_href => Hold the array parameter defaults as {REF}

    my ($arg_href) = @_;

    $parameter_href = $arg_href->{parameter_href};
    $FILEHANDLE     = $arg_href->{FILEHANDLE};

    my $pwd = cwd();

    ## Check installation path
    if ( not -d $parameter_href->{perl_install_dir} ) {

        # Create installation dir
        say {$FILEHANDLE} q{# Create perl installation dir};
        say {$FILEHANDLE} q{mkdir -p } . $parameter_href->{perl_install_dir}
          . qq{\n};
    }

    ## Writing wget command
    say {$FILEHANDLE} q{# Download perl from cpan};

    ## Download source
    my $url =
        q{http://www.cpan.org/src/5.0/perl-}
      . $parameter_href->{perl_version}
      . q{.tar.gz};

    ## Path to downloaded file
    my $outfile_path = catfile( $parameter_href->{perl_install_dir},
        q{perl-} . $parameter_href->{perl_version} . q{.tar.gz} );

    wget(
        {
            parameter_href => $parameter_href,
            url            => $url,
            outfile_path   => $outfile_path,
            FILEHANDLE     => $FILEHANDLE
        }
    );

    print {$FILEHANDLE} qq{\n};

    say {$FILEHANDLE} qq{cd } . $parameter_href->{perl_install_dir} . qq{\n};

    ## Writing unpack command
    say {$FILEHANDLE} q{# Unpack and remove tar file};

    my $tar_file = catfile( $parameter_href->{perl_install_dir},
        q{perl-} . $parameter_href->{perl_version} . q{.tar.gz} );

    say {$FILEHANDLE} q{tar -xzf }
      . $tar_file
      . q{ && rm }
      . $tar_file . qq{\n};

    ## Building perl
    say {$FILEHANDLE} q{# Build perl};

    my $perl_install_path = catdir(
        $parameter_href->{perl_install_dir},
        q{perl-} . $parameter_href->{perl_version}
    );

    say {$FILEHANDLE} q{cd } . $perl_install_path;
    say {$FILEHANDLE} q{./Configure -des -Dprefix=} . $perl_install_path;
    say {$FILEHANDLE} q{make};

    if ( not $parameter_href->{perl_skip_test} ) {
        say {$FILEHANDLE} q{make test};
    }
    say {$FILEHANDLE} q{make install};
    say {$FILEHANDLE} q{cd } . $pwd . qq{\n};

    ## Editing PATH and .bashrc
    say {$FILEHANDLE} q{# Edit PATH and .bashrc};
    say {$FILEHANDLE} q{echo 'export PATH=} . $perl_install_path
      . q{/:$PATH' >> ~/.bashrc};

    # Use newly installed perl
    say {$FILEHANDLE} q{export PATH=} . $perl_install_path . q{/:$PATH};

    # Add at start-up
    say {$FILEHANDLE} q{echo 'eval `perl -I }
      . catdir( $perl_install_path, qw{ lib perl5 } )
      . q{ -Mlocal::lib=}
      . $perl_install_path
      . q{/`' >> ~/.bash_profile } . qq{\n};

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
    wget(
        {
            parameter_href => $parameter_href,
            url            => $url,
            outfile_path   => $outfile_path,
            FILEHANDLE     => $FILEHANDLE
        }
    );

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

    return;
}

sub install_modules {

## insatall_modules

## Function  : Perl wrapper for writing cpanm recipe to $FILEHANDLE.
## Returns   : ""
## Arguments : $parameter_href, $FILEHANDLE
##           : $parameter_href => Hash with paramters {REF}
##           : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    $parameter_href = $arg_href->{parameter_href};
    $FILEHANDLE     = $arg_href->{FILEHANDLE};

    ## Base command
    my @commands = q{cpanm};

    ## Add optional force flag
    if ( $parameter_href->{force} ) {
        push @commands, q{--force};
    }

    ## Add optional quiet flag
    if ( $parameter_href->{quiet} ) {
        push @commands, q{--quiet};
    }

    ## Add optional verbose flag
    if ( $parameter_href->{verbose} ) {
        push @commands, q{--verbose};
    }

    ## Modules to install
    my @perl_modules;
    if ( @{ $parameter_href->{perl_modules_append} } ) {
        @perl_modules = @{ $parameter_href->{perl_modules_append} };
        @perl_modules = split( /,/, join( q{,}, @perl_modules ) );

        @perl_modules = ( @{ $parameter_href->{perl_modules} }, @perl_modules );

        # Remove any duplicate modules
        my %perl_modules = map { $_ => 1 } @perl_modules;
        @perl_modules = keys %perl_modules;
    }
    else {
        @perl_modules = @{ $parameter_href->{perl_modules} };
    }

    push @commands, join q{ }, @perl_modules;

    say {$FILEHANDLE} join q{ }, @commands;

    return;
}

sub wget {

## wget

## Function  : Perl wrapper for writing wget recipe to $FILEHANDLE or return commands array. Based on GNU Wget 1.12, a non-interactive network retriever.
## Returns   : "@commands"
## Arguments : $parameter_href, $FILEHANDLE
##           : $parameter_href => Hash that holds all parameters {REF}
##           : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten arguments
    my $url;
    my $outfile_path;

    ## Load arguments
    $parameter_href = $arg_href->{parameter_href};
    $url            = $arg_href->{url};
    $outfile_path   = $arg_href->{outfile_path};
    $FILEHANDLE     = $arg_href->{FILEHANDLE};

    # Stores commands depending on input parameters
    my @commands = qw(wget);

    if ( $parameter_href->{quiet} ) {
        push @commands, q{--quiet};
    }

    if ( $parameter_href->{verbose} ) {
        push @commands, q{--verbose};
    }
    else {
        push @commands, q{--no-verbose};
    }

    # Add url adress
    push @commands, $url;

    # Add outfile path
    push @commands, q{-O } . $outfile_path;

    print {$FILEHANDLE} join( q{ }, @commands ) . q{ };

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
