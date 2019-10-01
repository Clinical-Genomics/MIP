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
  qw{ $AT $DOLLAR_SIGN $DOUBLE_QUOTE $COLON $LOG_NAME $NEWLINE $SINGLE_QUOTE $SPACE };
use MIP::Gnu::Coreutils qw{ gnu_chmod gnu_mkdir gnu_echo };
use MIP::Language::Shell qw{ build_shebang };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Program::Singularity qw{ singularity_exec singularity_pull };
use MIP::Recipes::Install::Vep qw{ install_vep };
use MIP::Recipes::Install::Cadd qw{ install_cadd };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_singularity_containers };
}

sub install_singularity_containers {

## Function : Pull container from singularity hub or docker hub
## Returns  :
## Arguments: $active_parameter_href => Active parameter hash {REF}
##          : $conda_env             => Conda environemnt name
##          : $conda_env_path        => Path to conda environment
##          : $container_dir_path    => Pull containers to this path
##          : $container_href        => Hash with container {REF}
##          : $FILEHANDLE            => Filehandle
##          : $quiet                 => Optionally turn on quiet output
##          : $verbose               => Log debug messages

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $conda_env;
    my $conda_env_path;
    my $container_dir_path;
    my $container_href;
    my $FILEHANDLE;
    my $quiet;
    my $verbose;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        conda_env => {
            required    => 1,
            store       => \$conda_env,
            strict_type => 1,
        },
        conda_env_path => {
            required    => 1,
            store       => \$conda_env_path,
            strict_type => 1,
        },
        container_dir_path => {
            required    => 1,
            store       => \$container_dir_path,
            strict_type => 1,
        },
        container_href => {
            default     => {},
            required    => 1,
            store       => \$container_href,
            strict_type => 1,
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE,
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

    say {$FILEHANDLE} q{## Pull containers with Singularity};

    ## Set default path for containers
    if ( not $container_dir_path ) {
        $container_dir_path = catdir( $conda_env_path, qw{ share containers } );
    }

    ## Containers requiring something extra
    my %finish_container_installation = (
        cadd => \&install_cadd,
        vep  => \&install_vep,
    );

    ## Write check command to FILEHANDLE
    say {$FILEHANDLE} q{## Check for container path};
    my $dir_check = check_future_filesystem_for_directory(
        {
            directory_path => $container_dir_path,
        }
    );
    say {$FILEHANDLE} $dir_check . $NEWLINE;

  CONTAINER:
    foreach my $container ( keys %{$container_href} ) {

        $log->info( q{Writing instructions for pulling container}
              . $COLON
              . $SPACE
              . $container );

        say {$FILEHANDLE} q{## Setting up } . $container . q{ container};

        ## Create placeholder
        $container_href->{program_bind_paths} = [];

        my $container_path = catfile( $container_dir_path, $container . q{.sif} );
        singularity_pull(
            {
                container_uri => $container_href->{$container}{uri},
                FILEHANDLE    => $FILEHANDLE,
                force         => 1,
                outfile_path  => $container_path,
            }
        );
        print {$FILEHANDLE} $NEWLINE;

        ## Finishing touches for certain containers
        if ( $finish_container_installation{$container} ) {

            $finish_container_installation{$container}->(
                {
                    active_parameter_href => $active_parameter_href,
                    conda_env             => $conda_env,
                    conda_env_path        => $conda_env_path,
                    container_href        => $container_href,
                    container_path        => $container_path,
                    FILEHANDLE            => $FILEHANDLE,
                }
            );
        }

        ## Make available as exeutable in bin
        setup_singularity_executable(
            {
                conda_env_path         => $conda_env_path,
                container_path         => $container_path,
                executable_href        => $container_href->{$container}{executable},
                FILEHANDLE             => $FILEHANDLE,
                program_bind_paths_ref => $container_href->{program_bind_paths},
            }
        );
        print {$FILEHANDLE} $NEWLINE;

    }
    return 1;
}

sub setup_singularity_executable {

## Function : Make singularity program executable available in conda bin
## Returns  :
## Arguments: $conda_env_path         => Path to conda environment
##          : $container_path         => Path to container
##          : $executable_href        => Hash with executables and their path in the container (if not in PATH) {REF}
##          : $FILEHANDLE             => Filehandle
##          : $program_bind_paths_ref => Extra static bind paths

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_env_path;
    my $container_path;
    my $executable_href;
    my $FILEHANDLE;
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
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE,
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

    my $bash_command = $DOUBLE_QUOTE . $DOLLAR_SIGN . $AT . $DOUBLE_QUOTE;

  EXECUTABLE:
    foreach my $executable ( keys %{$executable_href} ) {

        my $container_executable = $executable;
        if ( $executable_href->{$executable} ) {
            $container_executable = $executable_href->{$executable};
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
                FILEHANDLE     => $FILEHANDLE,
                outfile_path   => $proxy_executable_path,
                string_wrapper => $SINGLE_QUOTE,
                strings_ref    => \@shebang,
            }
        );
        print {$FILEHANDLE} $NEWLINE;
        gnu_echo(
            {
                FILEHANDLE             => $FILEHANDLE,
                no_trailing_newline    => 1,
                stdoutfile_path_append => $proxy_executable_path,
                string_wrapper         => $SINGLE_QUOTE,
                strings_ref            => [ join $SPACE, @singularity_cmds ],
            }
        );
        print {$FILEHANDLE} $NEWLINE;
        gnu_chmod(
            {
                file_path  => $proxy_executable_path,
                FILEHANDLE => $FILEHANDLE,
                permission => q{a+x},
            }
        );
        print {$FILEHANDLE} $NEWLINE;
    }
    return 1;
}

1;
