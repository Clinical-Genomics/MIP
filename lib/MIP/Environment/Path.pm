package MIP::Environment::Path;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib
use MIP::Constants qw{ $LOG_NAME };

## Constants
Readonly my $MINUS_TWO => -2;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_bin_file_path get_conda_path is_binary_in_path };
}

## Constants
Readonly my $MINUS_ONE => -1;

sub get_bin_file_path {

## Function : Get the absolute path to the binary file
##          : and the corresponding conda environment
## Returns  : $bin_file_path, $conda_environment
## Arguments: $bin_file         => Name of binary file
##          : $conda_path       => Path to conda directory
##          : $environment_href => Hash with programs and their environments {REF}
##          : $environment_key  => Key to the environment_href [program|recipe_name]

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bin_file;
    my $conda_path;
    my $environment_href;
    my $environment_key;

    my $tmpl = {
        bin_file => {
            defined     => 1,
            required    => 1,
            store       => \$bin_file,
            strict_type => 1,
        },
        conda_path => {
            defined     => 1,
            required    => 1,
            store       => \$conda_path,
            strict_type => 1,
        },
        environment_href => {
            default     => {},
            required    => 1,
            store       => \$environment_href,
            strict_type => 1,
        },
        environment_key => {
            defined     => 1,
            required    => 1,
            store       => \$environment_key,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Cwd qw{ abs_path };

    ## Get environment where binary file is located by using "conda activate <conda_env>"
    my $conda_environment = @{ $environment_href->{$environment_key} }[$MINUS_ONE];

    ## Place it in the conda environment binary path
    my $bin_file_path =
      catfile( $conda_path, q{envs}, $conda_environment, q{bin}, $bin_file );

    ## Return absolute path and conda environment for binary file
    return ( abs_path($bin_file_path), $conda_environment );
}

sub get_conda_path {

## Function : Get path to conda directory
## Returns  : $conda_path
## Arguments: $bin_file => Bin file to test

    my ($arg_href) = @_;

    ## Default(s)
    my $bin_file;

    my $tmpl = {
        bin_file => {
            default     => q{conda},
            store       => \$bin_file,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use IPC::Cmd qw{ can_run };

    ## Find path to conda bin
    my $conda_path = can_run($bin_file);

    return if ( not $conda_path );

    ## Split dirs to array
    my @conda_path_dirs = File::Spec->splitdir($conda_path);

    ## Traverse to conda dir from binary
    splice @conda_path_dirs, $MINUS_TWO;

    ## Return path to conda main directory
    return catdir(@conda_path_dirs);
}

sub is_binary_in_path {

## Function : Test if binary is in path
## Returns  : 1
## Arguments: $binary => Binary to test

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary;

    my $tmpl = {
        binary => {
            defined     => 1,
            required    => 1,
            store       => \$binary,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use IPC::Cmd qw{ can_run };

    if ( can_run($binary) ) {

        ## Broadcast successful scan through PATH for supplied binary
        _check_binary_broadcast_pass(
            {
                binary => $binary,
            }
        );
        return 1;
    }

    ## Broadcast scan through PATH for supplied binary when not found
    _check_binary_broadcast_fail(
        {
            binary => $binary,
        }
    );
    exit 1;

}

sub _check_binary_broadcast_fail {

## Function : Broadcast scan through PATH for supplied binary when not found
## Returns  :
## Arguments: $binary => Binary to search for

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary;

    my $tmpl = {
        binary => {
            defined     => 1,
            required    => 1,
            store       => \$binary,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Broadcast binary not found
    $log->fatal( q{Could not detect } . $binary . q{ in PATH} );

    return;
}

sub _check_binary_broadcast_pass {

## Function  : Broadcast successful scan through PATH for supplied binary
## Returns   :
## Arguments : $binary => Binary to search for

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary;

    my $tmpl = {
        binary => {
            defined     => 1,
            required    => 1,
            store       => \$binary,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Broadcast if found
    $log->info( q{Program check: } . $binary . q{ in PATH} );
    return 1;

}

1;
