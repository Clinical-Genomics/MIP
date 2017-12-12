package MIP::Package_manager::Conda;

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
    our $VERSION = 1.08;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ conda_check_env_status conda_create conda_install conda_source_activate conda_source_deactivate conda_uninstall conda_update };
}

sub conda_create {

##Function : Create Conda environment
##Returns  : @commands
##Arguments: $conda_channel   => Search for packages in specified conda channel
##         : $env_name        => Name of environment to create
##         : $FILEHANDLE      => Filehandle to write to
##         : $no_confirmation => Do not ask for confirmation
##         : $packages_ref    => Packages to be installed
##         : $quiet           => Do not display progress bar

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_channel;
    my $env_name;
    my $FILEHANDLE;
    my $no_confirmation;
    my $packages_ref;
    my $quiet;

    my $tmpl = {
        conda_channel => {
            defined     => 1,
            strict_type => 1,
            store       => \$conda_channel
        },
        env_name => {
            default     => q{},
            strict_type => 1,
            store       => \$env_name
        },
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
        packages_ref => {
            default     => [],
            strict_type => 1,
            store       => \$packages_ref
        },
        quiet => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$quiet
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

    if ($conda_channel) {
        push @commands, q{--channel} . $SPACE . $conda_channel;
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

##Function : Activate conda environment
##Returns  : "@commands"
##Arguments: $FILEHANDLE => Filehandle to write to
##         : $env_name   => Name of conda environment

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $env_name;
    my $FILEHANDLE;

    my $tmpl = {
        env_name => {
            required    => 1,
            strict_type => 1,
            store       => \$env_name
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE
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

##Function : Deactivate conda environment
##Returns  : "@commands"
##Arguments: $FILEHANDLE => Filehandle to write to

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

## Function  : Update conda
## Returns   : @commands
## Arguments : $FILEHANDLE      => Filehandle to write to
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

## Function  : Check if a conda environment is active (returns name of env if true).
## Returns   :
## Arguments : $log => Log

    my ($arg_href) = @_;

    ## Flatten arguments
    my $log;

    my $tmpl = {
        log => {
            required => 1,
            defined  => 1,
            store    => \$log
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

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
    #   Output is captured in $env_status.
    my $env_status;
    run(
        command => qq{conda info --envs | $detect_active_conda_env},
        buffer  => \$env_status
    );

    # Kill script if a conda environment is active
    if ($env_status) {
        $log->fatal( q{Found activate conda env: } . $env_status );
        $log->fatal(
            q{Run 'source deactivate' prior to running installation script});
        exit 1;
    }

    return;
}

sub conda_install {

##Function : Install packages into conda environment
##Returns  : @commands
##Arguments: $conda_channel   => Search for packages in specified conda channel
##         : $env_name        => Name of environment to create
##         : $FILEHANDLE      => Filehandle to write to
##         : $no_confirmation => Do not ask for confirmation
##         : $packages_ref    => Packages to be installed
##         : $quiet           => Do not display progress bar

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_channel;
    my $env_name;
    my $FILEHANDLE;
    my $no_confirmation;
    my $packages_ref;
    my $quiet;

    my $tmpl = {
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
        no_confirmation => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$no_confirmation
        },
        packages_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$packages_ref
        },
        quiet => {
            default     => 1,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$quiet
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

sub conda_uninstall {

##Function : Uninstall packages from conda environment
##Returns  : @commands
##Arguments: $env_name        => Name of environment to create
##         : $FILEHANDLE      => Filehandle to write to
##         : $no_confirmation => Do not ask for confirmation
##         : $packages_ref    => Packages to be installed
##         : $quiet           => Do not display progress bar
##         : $verbose         => Verbose output

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $env_name;
    my $FILEHANDLE;
    my $no_confirmation;
    my $packages_ref;
    my $quiet;
    my $verbose;

    my $tmpl = {
        env_name => {
            default     => undef,
            strict_type => 1,
            store       => \$env_name
        },
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
        packages_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$packages_ref
        },
        quiet => {
            default     => 1,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$quiet
        },
        verbose => {
            default     => 1,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$verbose
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = q{conda uninstall};

    if ($env_name) {
        push @commands, q{--name} . $SPACE . $env_name;
    }

    # Do not display progress bar
    if ($quiet) {
        push @commands, q{--quiet};
    }

    if ($verbose) {
        push @commands, q{--verbose};
    }

    if ($no_confirmation) {
        push @commands, q{--yes};
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
