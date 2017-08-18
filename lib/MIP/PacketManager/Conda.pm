package MIP::PacketManager::Conda;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

use Getopt::Long;
use Cwd;
use Cwd qw{ abs_path };
use FindBin qw{ $Bin };               #Find directory of script
use IO::Handle;
use File::Basename qw{ dirname basename fileparse };
use File::Spec::Functions qw{ catfile catdir devnull };
use Readonly;

## MIPs lib/
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

## Constants
Readonly my $SPACE   => q{ };
Readonly my $NEWLINE => qq{\n};

BEGIN {

    use base qw{Exporter};
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{conda_create};
}

sub conda_create {

##conda_create

##Function : Create Conda environment
##Returns  : @commands
##Arguments: $packages_ref, $env_name, $python_version, $quiet, $no_confirmation,
##           $FILEHANDLE
##         : $packages_ref    => Packages to be installed
##         : $FILEHANDLE      => Filehandle to write to
##         : $env_name        => Name of environment to create
##         : $python_version  => Python version to install
##         : $quiet           => Do not display progress bar
##         : $no_confirmation => Do not ask for confirmation

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $env_name;
    my $python_version;
    my $quiet;
    my $no_confirmation;
    my $packages_ref;
    my $FILEHANDLE;

    my $tmpl = {
        packages_ref => {
            default     => [],
            strict_type => 1,
            store       => \$packages_ref
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE
        },
        env_name => {
            default     => q{},
            defined     => 1,
            strict_type => 1,
            store       => \$env_name
        },
        python_version => {
            default     => '2.7',
            allow       => [qr/\d+.\d+ | \d+.\d+.\d+/xms],
            strict_type => 1,
            store       => \$python_version
        },
        quiet => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$quiet
        },
        no_confirmation => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$no_confirmation
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ##Stores commands depending on input parameters
    # Basic command
    my @commands = q{conda create};

    # Add python version
    if ($python_version) {
        push @commands, q{python=} . $python_version;
    }

    if ($env_name) {
        push @commands, q{--name } . $env_name;
    }

    if ($quiet) {
        push @commands, q{--quiet};
    }

    if ($no_confirmation) {
        push @commands, q{--yes};
    }

    if ( @{$packages_ref} ) {
        push @commands, join $SPACE, @{$packages_ref};
    }

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

1;
