package MIP::Package_manager::Conda;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use Getopt::Long;
use IO::Handle;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## Cpanm
use IPC::Cmd qw{ run };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $NEWLINE $SPACE };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {

    use base qw{Exporter};
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.16;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ conda_activate conda_check_env_status conda_create conda_deactivate conda_install };
}

sub conda_activate {

##Function : Activate conda environment
##Returns  : @commands
##Arguments: $env_name   => Name of conda environment
##         : $FILEHANDLE => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $env_name;
    my $FILEHANDLE;

    my $tmpl = {
        env_name => {
            required    => 1,
            store       => \$env_name,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters

    # Basic command
    my @commands = qw{ conda activate };

    # Activates env, default base
    push @commands, $env_name;

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub conda_check_env_status {

## Function  : Check if a conda environment is active (returns name of env if true).
## Returns   :
## Arguments : $disable_env_check => Disable environment check

    my ($arg_href) = @_;

    ## Default(s)
    my $disable_env_check;

    my $tmpl = {
        disable_env_check => {
            default     => 0,
            allow       => [ 0, 1 ],
            store       => \$disable_env_check,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ### Require deactivate any activate env prior to installation

    ## Unless the active environment is root the expression will return true
    ##   and print the environment name
    my $detect_active_conda_env =
      q?perl -nae 'if( ($_!~/ ^root | ^base /xms) && ($_=~/\*/) ) {print $F[0]}'?;

    # Pipes the output from the shell command "conda info --envs"
    #   to $detect_active_conda_env.
    #   Output is captured in $env_status.
    my $env_status;
    run(
        command => qq{conda info --envs | $detect_active_conda_env},
        buffer  => \$env_status
    );

    # Exit if a conda environment is active
    if ($env_status) {

        $log->warn( q{Found activated conda env: } . $env_status );

        ## Mainly used for running test script in activated
        ## env not actual install
        if ( not $disable_env_check ) {

            $log->fatal(q{Run 'conda deactivate' prior to running installation script});
            exit 1;
        }
    }
    return;
}

sub conda_create {

##Function : Create Conda environment
##Returns  : @commands
##Arguments: $conda_channels_ref => Search for packages in specified conda channels {REF}
##         : $env_name           => Name of environment to create
##         : $FILEHANDLE         => Filehandle to write to
##         : $no_confirmation    => Do not ask for confirmation
##         : $packages_ref       => Packages to be installed
##         : $quiet              => Do not display progress bar

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_channels_ref;
    my $env_name;
    my $FILEHANDLE;
    my $no_confirmation;
    my $packages_ref;
    my $quiet;

    my $tmpl = {
        conda_channels_ref => {
            allow => sub {
                my $channels_ref = shift @_;
                ## Allow undef
                return 1 if scalar @{$channels_ref} == 0;
                ## Test values
                return _check_array_membership(
                    {
                        allowed_elements_ref => [qw{ bioconda conda-forge }],
                        test_elements_ref    => $channels_ref,
                    }
                );
            },
            default     => [],
            store       => \$conda_channels_ref,
            strict_type => 1,
        },
        env_name => {
            default     => q{},
            store       => \$env_name,
            strict_type => 1,
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE,
        },
        no_confirmation => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$no_confirmation,
            strict_type => 1,
        },
        packages_ref => {
            default     => [],
            store       => \$packages_ref,
            strict_type => 1,
        },
        quiet => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$quiet,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ##Stores commands depending on input parameters
    # Basic command
    my @commands = qw{ conda create };

    if ($env_name) {
        push @commands, q{--name} . $SPACE . $env_name;
    }

    if ($quiet) {
        push @commands, q{--quiet};
    }

    if ($no_confirmation) {
        push @commands, q{--yes};
    }

    if ( @{$conda_channels_ref} ) {
        push @commands,
          q{--channel} . $SPACE . join $SPACE . q{--channel} . $SPACE,
          @{$conda_channels_ref};
    }

    if ( @{$packages_ref} ) {
        push @commands, join $SPACE, @{$packages_ref};
    }

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );

    return @commands;
}

sub conda_deactivate {

##Function : Deactivate conda environment
##Returns  : @commands
##Arguments: $FILEHANDLE => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ conda deactivate };

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
        },
    );

    return @commands;
}

sub conda_install {

##Function : Install packages into conda environment
##Returns  : @commands
##Arguments: $conda_channels_ref => Search for packages in specified conda channels {REF}
##         : $env_name           => Name of environment to create
##         : $FILEHANDLE         => Filehandle to write to
##         : $no_confirmation    => Do not ask for confirmation
##         : $no_update_dep     => Only update dependencies that are required for the package to function
##         : $packages_ref       => Packages to be installed
##         : $quiet              => Do not display progress bar

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_channels_ref;
    my $env_name;
    my $FILEHANDLE;
    my $no_confirmation;
    my $no_update_dep;
    my $packages_ref;
    my $quiet;

    my $tmpl = {
        conda_channels_ref => {
            default => [],
            allow   => sub {
                my $channels_ref = shift @_;
                ## Allow undef
                return 1 if scalar @{$channels_ref} == 0;
                ## Test values
                return _check_array_membership(
                    {
                        allowed_elements_ref => [qw{ bioconda conda-forge }],
                        test_elements_ref    => $channels_ref,
                    }
                );
            },
            store       => \$conda_channels_ref,
            strict_type => 1,
        },
        env_name => {
            default     => undef,
            store       => \$env_name,
            strict_type => 1,
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE,
        },
        no_confirmation => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$no_confirmation,
            strict_type => 1,
        },
        no_update_dep => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$no_update_dep,
            strict_type => 1,
        },
        packages_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$packages_ref,
            strict_type => 1,
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$quiet,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ conda install };

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

    if ($no_update_dep) {

        push @commands, q{--no-update-deps};
    }

    if ( @{$conda_channels_ref} ) {

        push @commands,
          q{--channel} . $SPACE . join $SPACE . q{--channel} . $SPACE,
          @{$conda_channels_ref};
    }

    push @commands, join $SPACE, @{$packages_ref};

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub _check_array_membership {

##Function : Checks if all array elements are part of an array with allowed elements. Returns true/false
##Returns  : Boolean
##Arguments: $allowed_elements_ref => Allowed elements {REF}
##         : $test_elements_ref    => Array elements to test {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $allowed_elements_ref;
    my $test_elements_ref;

    my $tmpl = {
        allowed_elements_ref => {
            default     => [],
            required    => 1,
            store       => \$allowed_elements_ref,
            strict_type => 1,
        },
        test_elements_ref => {
            default     => [],
            required    => 1,
            store       => \$test_elements_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use List::Util qw{ any };

    ## Test values
  TEST_ELEMENT:
    foreach my $test_element ( @{$test_elements_ref} ) {

        return if not any { $_ eq $test_element } @{$allowed_elements_ref};
    }
    return 1;
}
1;
