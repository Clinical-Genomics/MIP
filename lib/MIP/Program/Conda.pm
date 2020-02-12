package MIP::Program::Conda;

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

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $NEWLINE $SPACE };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {

    use base qw{Exporter};
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.19;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ conda_activate conda_create conda_deactivate conda_install };
}

sub conda_activate {

##Function : Activate conda environment
##Returns  : @commands
##Arguments: $env_name   => Name of conda environment
##         : $filehandle => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $env_name;
    my $filehandle;

    my $tmpl = {
        env_name => {
            required    => 1,
            store       => \$env_name,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
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
            filehandle   => $filehandle,
        }
    );

    return @commands;
}

sub conda_create {

##Function : Create Conda environment
##Returns  : @commands
##Arguments: $conda_channels_ref => Search for packages in specified conda channels {REF}
##         : $env_name           => Name of environment to create
##         : $filehandle         => Filehandle to write to
##         : $no_confirmation    => Do not ask for confirmation
##         : $packages_ref       => Packages to be installed
##         : $quiet              => Do not display progress bar

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_channels_ref;
    my $env_name;
    my $filehandle;
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
        filehandle => {
            required => 1,
            store    => \$filehandle,
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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );

    return @commands;
}

sub conda_deactivate {

##Function : Deactivate conda environment
##Returns  : @commands
##Arguments: $filehandle => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ conda deactivate };

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
        },
    );

    return @commands;
}

sub conda_install {

##Function : Install packages into conda environment
##Returns  : @commands
##Arguments: $conda_channels_ref => Search for packages in specified conda channels {REF}
##         : $env_name           => Name of environment to create
##         : $filehandle         => Filehandle to write to
##         : $no_confirmation    => Do not ask for confirmation
##         : $no_update_dep     => Only update dependencies that are required for the package to function
##         : $packages_ref       => Packages to be installed
##         : $quiet              => Do not display progress bar

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_channels_ref;
    my $env_name;
    my $filehandle;
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
        filehandle => {
            required => 1,
            store    => \$filehandle,
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
            filehandle   => $filehandle,
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
