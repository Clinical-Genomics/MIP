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
use IO::Handle;
use Readonly;
use IPC::Cmd qw{ can_run run };

## MIPs lib/
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

## Constants
Readonly my $SPACE   => q{ };
Readonly my $NEWLINE => qq{\n};

BEGIN {

    use base qw{Exporter};
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ conda_create conda_source_activate conda_source_deactivate
      conda_update conda_check_env_status conda_install };
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
##         : $quiet           => Do not display progress bar
##         : $no_confirmation => Do not ask for confirmation

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $env_name;
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
            strict_type => 1,
            store       => \$env_name
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

    if ($env_name) {
        push @commands, q{--name} . $SPACE . $env_name;
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

sub conda_source_activate {

##conda_source_activate

##Function : Activate conda environment
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $env_name,
##         : $FILEHANDLE => Filehandle to write to
##         : $env_name   => Name of conda environment

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $env_name;

    my $tmpl = {
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE
        },
        env_name => {
            required    => 1,
            strict_type => 1,
            store       => \$env_name
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters

    # Basic command
    my @commands = q{source activate};

    # Activates env, default root
    push @commands, $env_name;

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub conda_source_deactivate {

## conda_source_deactivate

##Function : Deactivate conda environment
##Returns  : "@commands"
##Arguments: $FILEHANDLE
##         : $FILEHANDLE => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;

    my $tmpl = {
        FILEHANDLE => {
            required => 1,
            defined  => 1,
            store    => \$FILEHANDLE
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = q{source deactivate};

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
        },
    );

    return @commands;
}

sub conda_update {

## conda_update

## Function  : Update conda
## Returns   : @commands
## Arguments : $FILEHANDLE, $no_confirmation
##           : $FILEHANDLE      => Filehandle to write to
##           : $no_confirmation => Do not ask for confirmation

    my ($arg_href) = @_;

    ## Flatten arguments
    my $FILEHANDLE;
    my $no_confirmation;

    my $tmpl = {
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE
        },
        no_confirmation => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$no_confirmation
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = q{conda update};

    if ($no_confirmation) {
        push @commands, q{--yes};
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

sub conda_check_env_status {

## conda_check_env_status

## Function  : Check if a conda environment is active (exit if true).
## Returns   :
## Arguments :

    ## Deactivate any activate env prior to installation
    #   Perl options:
    #   n : loop over input
    #   a : automatically split input and store in array @F
    #   e : execute code
    #
    # Unless the active environment is root the expression will return true
    #   and print the environment name
    my $detect_active_conda_env =
      q?perl -nae 'if( ($_!~/^root/) && ($_=~/\*/) ) {print $F[0]}'?;

    # Pipes the output from the shell command "conda info --envs"
    #   to $detect_active_conda_env.
    #   Output is captured in $run_output.
    my $run_output;
    run(
        command => qq{conda info --envs | $detect_active_conda_env},
        buffer  => \$run_output
    );

    if ($run_output) {
        say STDOUT q{Found activated conda env:} . $SPACE . $run_output;
        say STDOUT q{Please exit conda env:} . $SPACE . $run_output . $SPACE
          . q{with 'source deactivate' before executing install script};
        exit 1;
    }

    return;
}

sub conda_install {

## conda_install

##Function : Install packages into conda environment
##Returns  : @commands
##Arguments: $packages_ref, $conda_channel, $FILEHANDLE, $env_name, $quiet, $no_cofirmation

##         : $packages_ref    => Packages to be installed
##         : $conda_channel   => Search for packages in specified conda channel
##         : $FILEHANDLE      => Filehandle to write to
##         : $env_name        => Name of environment to create
##         : $quiet           => Do not display progress bar
##         : $no_confirmation => Do not ask for confirmation

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $packages_ref;
    my $conda_channel;
    my $env_name;
    my $FILEHANDLE;
    my $quiet;
    my $no_confirmation;

    my $tmpl = {
        packages_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$packages_ref
        },
        conda_channel => {
            defined     => 1,
            strict_type => 1,
            store       => \$conda_channel
        },
        env_name => {
            default     => undef,
            strict_type => 1,
            store       => \$env_name
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE
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

    my @commands = q{conda install};

    if ($env_name) {
        push @commands, q{--name} . $SPACE . $env_name;
    }

    if ($quiet) {

        #Do not display progress bar
        push @commands, q{--quiet};
    }

    if ($no_confirmation) {
        push @commands, q{--yes};
    }

    if ($conda_channel) {
        push @commands, q{--channel} . $SPACE . $conda_channel;
    }

    push @commands, join $SPACE, @{$packages_ref};

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
