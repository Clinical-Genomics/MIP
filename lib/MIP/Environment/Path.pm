package MIP::Environment::Path;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile splitdir };
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
use MIP::Constants qw{ $LOG_NAME $SPACE };

## Constants
Readonly my $MINUS_TWO => -2;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_binary_in_path
      get_bin_file_path
      get_conda_bin_dir_path
      get_conda_path
      is_binary_in_path
      reduce_dir_paths
    };
}

## Constants
Readonly my $MINUS_ONE => -1;

sub check_binary_in_path {

## Function : Scans through PATH for supplied binary
## Returns  : $binary_path
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $binary                => Binary to search for
##          : $program_name          => MIP program name (Analysis recipe switch)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $binary;
    my $program_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        binary => {
            defined     => 1,
            required    => 1,
            store       => \$binary,
            strict_type => 1,
        },
        program_name => {
            store       => \$program_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Search for binary in PATH in any MIP conda env defined by config
    ## or conda base
    my $env_binary_path = get_conda_bin_dir_path(
        {
            active_parameter_href => $active_parameter_href,
            bin_file              => $binary,
            environment_key       => $program_name,
        }
    );

    ## Potential full path to binary
    my $binary_path = catfile( $env_binary_path, $binary );

    ## Already tested
    return if ( exists $active_parameter_href->{binary_path}{$binary} );

    ## Test binary
    is_binary_in_path( { binary => $binary_path, } );

    return $binary_path;
}

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

sub get_conda_bin_dir_path {

## Function : Attempts to find path to directory with binary in conda env
## Returns  : $conda_bin_dir_path
## Arguments: $active_parameter_href => Active parameter hash {REF}
##          : $bin_file              => Bin file to test
##          : $environment_key       => Key to conda environment

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $bin_file;
    my $environment_key;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        bin_file => {
            defined     => 1,
            required    => 1,
            store       => \$bin_file,
            strict_type => 1,
        },
        environment_key => {
            store       => \$environment_key,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ get_package_env_attributes };
    use MIP::Environment::Manager qw{ get_env_method_cmds };
    use MIP::Environment::Path qw{ get_bin_file_path };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack
    my $conda_path = $active_parameter_href->{conda_path};
    my ( $env_name, $env_method );

    ## Get environment name and manager in use for $environment_key
  ENV_KEY:
    foreach my $env_key ( $environment_key, qw{ mip } ) {

        ( $env_name, $env_method ) = get_package_env_attributes(
            {
                load_env_href => $active_parameter_href->{load_env},
                package_name  => $env_key,
            }
        );
        ## Found program|recipe within env
        last if ($env_name);
    }

    ## Get env load command
    my @env_method_cmds = get_env_method_cmds(
        {
            action     => q{load},
            env_method => $env_method,
            env_name   => $env_name,
        }
    );

    ## Add to environment hash with "recipe_name" as keys and "source env command" as value
    my %environment = ( $environment_key => [@env_method_cmds] );

    ## Get the bin file path
    my ( $bin_file_path, $conda_env ) = get_bin_file_path(
        {
            bin_file         => $bin_file,
            conda_path       => $conda_path,
            environment_href => \%environment,
            environment_key  => $environment_key,
        }
    );

    ## Test if path exists
    if (   not $bin_file_path
        or not -f $bin_file_path )
    {
        $log->logcroak( q{Failed to find default path for}
              . $SPACE
              . $bin_file
              . $SPACE
              . q{in conda environment}
              . $SPACE
              . $conda_env );
    }

    ## Remove bin file from path
    my $conda_bin_dir_path = dirname($bin_file_path);

    return $conda_bin_dir_path;
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

sub reduce_dir_paths {

## Function : Parses directory paths and reduces them to a non-overlapping array. No check for existing files or directories
## Returns  : @reduced_dir_paths
## Arguments: $dir_paths_ref => Directory paths to parse {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $dir_paths_ref;

    my $tmpl = {
        dir_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$dir_paths_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @dir_paths;

    ## Split to dir path to array
  DIR_PATH:
    foreach my $dir_path ( @{$dir_paths_ref} ) {

        next DIR_PATH if ( not defined $dir_path );

        push @dir_paths, [ splitdir($dir_path) ];
    }

    ## Sort according to size
    @dir_paths = sort { @{$a} <=> @{$b} } @dir_paths;

    ## Reformat to strings
    @dir_paths = map { catdir( @{$_} ) } @dir_paths;

    my @reduced_dir_paths;

  BIND_PATH:
    while (@dir_paths) {

        ## Shift array
        my $dir_path = shift @dir_paths;

        ## Save path
        push @reduced_dir_paths, $dir_path;

        ## Get indexes of all the paths in the array that have an identical beginning to one we are testing
        ## The \Q and \E in the regex turns of interpolation
        my @match_idxs =
          grep { $dir_paths[$_] =~ / ^\Q$dir_path\E.* /xms } 0 .. $#dir_paths;

      MATCH_IDX:
        foreach my $match_idx ( reverse @match_idxs ) {

            ## Remove those paths with matching starts
            splice @dir_paths, $match_idx, 1;
        }
    }

    return @reduced_dir_paths;
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
