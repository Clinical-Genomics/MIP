package MIP::Recipes::Install::Singularity;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile devnull };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Check::Path qw{ check_future_filesystem_for_directory };
use MIP::Constants
  qw{ $AT $BACKWARD_SLASH $DOLLAR_SIGN $DOUBLE_QUOTE $COLON $EMPTY_STR $LOG_NAME $NEWLINE $SINGLE_QUOTE $SPACE };
use MIP::Gnu::Bash qw{ gnu_unset };
use MIP::Gnu::Coreutils qw{ gnu_chmod gnu_echo gnu_mkdir };
use MIP::Language::Shell qw{ build_shebang };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Program::Singularity qw{ singularity_exec singularity_pull };
use MIP::Recipes::Install::Htslib qw{ install_htslib };
use MIP::Recipes::Install::Cadd qw{ install_cadd };
use MIP::Recipes::Install::Vep qw{ install_vep };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_singularity_containers };
}

sub install_singularity_containers {

## Function : Pull container from singularity hub or docker hub
## Returns  :
## Arguments: $active_parameter_href => Active parameter hash {REF}
##          : $conda_env_path        => Path to conda environment
##          : $container_href        => Hash with container {REF}
##          : $filehandle            => Filehandle
##          : $quiet                 => Optionally turn on quiet output
##          : $verbose               => Log debug messages

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $conda_env_path;
    my $container_href;
    my $filehandle;
    my $quiet;
    my $verbose;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        conda_env_path => {
            required    => 1,
            store       => \$conda_env_path,
            strict_type => 1,
        },
        container_href => {
            default     => {},
            required    => 1,
            store       => \$container_href,
            strict_type => 1,
        },
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            store       => \$quiet,
            strict_type => 1,
        },
        verbose => {
            allow       => [ undef, 0, 1 ],
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => $LOG_NAME,
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    ## Return if no containers
    return if not keys %{$container_href};

    say {$filehandle} q{## Pull containers with Singularity};

    ## Set dir path for containers
    my $container_dir_path = catdir( $conda_env_path, qw{ share containers } );

    ## Containers requiring something extra
    my %finish_container_installation = (
        htslib => \&install_htslib,
        cadd   => \&install_cadd,
        vep    => \&install_vep,
    );

    ## Write check command to filehandle
    say {$filehandle} q{## Check for container path};
    my $dir_check = check_future_filesystem_for_directory(
        {
            directory_path => $container_dir_path,
        }
    );
    say {$filehandle} $dir_check . $NEWLINE;

  CONTAINER:
    foreach my $container ( keys %{$container_href} ) {

        $log->info( q{Writing instructions for pulling container}
              . $COLON
              . $SPACE
              . $container );

        say {$filehandle} q{## Setting up } . $container . q{ container};

        ## Create placeholder
        $container_href->{program_bind_paths} = [];

        my $container_path = catfile( $container_dir_path, $container . q{.sif} );

        ## Place relative to conda proxy bin
        my $relative_container_path =
          catdir( qw{ \"\$CONDA_ENV_DIR\" share containers }, $container . q{.sif} );

        singularity_pull(
            {
                container_uri => $container_href->{$container}{uri},
                filehandle    => $filehandle,
                force         => 1,
                outfile_path  => $container_path,
            }
        );
        print {$filehandle} $NEWLINE;

        ## Finishing touches for certain containers
        if ( $finish_container_installation{$container} ) {

            $finish_container_installation{$container}->(
                {
                    active_parameter_href => $active_parameter_href,
                    container_href        => $container_href,
                    container_path        => $container_path,
                    filehandle            => $filehandle,
                }
            );
        }

        ## Make available as exeuctable in bin
        setup_singularity_executable(
            {
                conda_env_path         => $conda_env_path,
                container_path         => $relative_container_path,
                executable_href        => $container_href->{$container}{executable},
                filehandle             => $filehandle,
                program_bind_paths_ref => $container_href->{program_bind_paths},
            }
        );
        print {$filehandle} $NEWLINE;

    }
    return 1;
}

sub setup_singularity_executable {

## Function : Make singularity program executable available in conda bin
## Returns  :
## Arguments: $conda_env_path         => Path to conda environment
##          : $container_path         => Path to container
##          : $executable_href        => Hash with executables and their path in the container (if not in PATH) {REF}
##          : $filehandle             => Filehandle
##          : $program_bind_paths_ref => Extra static bind paths

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_env_path;
    my $container_path;
    my $executable_href;
    my $filehandle;
    my $program_bind_paths_ref;

    my $tmpl = {
        conda_env_path => {
            required    => 1,
            store       => \$conda_env_path,
            strict_type => 1,
        },
        container_path => {
            required    => 1,
            store       => \$container_path,
            strict_type => 1,
        },
        executable_href => {
            default     => {},
            required    => 1,
            store       => \$executable_href,
            strict_type => 1,
        },
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        program_bind_paths_ref => {
            default     => [],
            store       => \$program_bind_paths_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @shebang = build_shebang(
        {
            bash_bin_path => catfile( dirname( dirname( devnull() ) ), qw{ bin bash } ),
        }
    );

    my $bash_command =
        $BACKWARD_SLASH
      . $DOUBLE_QUOTE
      . $BACKWARD_SLASH
      . $DOLLAR_SIGN
      . $AT
      . $BACKWARD_SLASH
      . $DOUBLE_QUOTE;

    ## Find proxy bin directory
    my $executable_dir_cmd = q{DIR=$(dirname "$(readlink -f "$0")")};
    ## Assign conda env dir to enable relative path to proxy bin within conda env
    my $conda_env_dir_cmd = q{CONDA_ENV_DIR="$(dirname "$DIR")"};

  EXECUTABLE:
    foreach my $executable ( keys %{$executable_href} ) {

        my $container_executable = $executable;
        if ( $executable_href->{$executable} ) {
            $container_executable = $executable_href->{$executable};
        }

        if (    $executable_href->{$executable}
            and $executable_href->{$executable} eq q{no_executable_in_image} )
        {
            $container_executable = $EMPTY_STR;
        }
        my @singularity_cmds = singularity_exec(
            {
                bind_paths_ref                 => $program_bind_paths_ref,
                singularity_container          => $container_path,
                singularity_container_cmds_ref => [$container_executable],
            }
        );
        push @singularity_cmds, $bash_command;

        my $proxy_executable_path = catfile( $conda_env_path, q{bin}, $executable );

        gnu_echo(
            {
                filehandle     => $filehandle,
                outfile_path   => $proxy_executable_path,
                string_wrapper => $SINGLE_QUOTE,
                strings_ref    => \@shebang,
            }
        );
        print {$filehandle} $NEWLINE;

      CMD:
        foreach my $cmd ( $executable_dir_cmd, $conda_env_dir_cmd ) {

            gnu_echo(
                {
                    filehandle             => $filehandle,
                    stdoutfile_path_append => $proxy_executable_path,
                    string_wrapper         => $SINGLE_QUOTE,
                    strings_ref            => [$cmd],
                }
            );
            print {$filehandle} $NEWLINE;
        }
        gnu_echo(
            {
                filehandle             => $filehandle,
                no_trailing_newline    => 1,
                stdoutfile_path_append => $proxy_executable_path,
                strings_ref            => [ join $SPACE, @singularity_cmds ],
            }
        );
        print {$filehandle} $NEWLINE;
        gnu_chmod(
            {
                file_path  => $proxy_executable_path,
                filehandle => $filehandle,
                permission => q{a+x},
            }
        );
        print {$filehandle} $NEWLINE;

        ## Unset VEP_VERSION
        if ( $executable eq q{vep} ) {

            gnu_unset(
                {
                    bash_variable => q{VEP_VERSION},
                    filehandle    => $filehandle,
                }
            );
            print {$filehandle} $NEWLINE;
        }
    }
    return 1;
}

1;
