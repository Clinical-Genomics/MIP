package MIP::Recipes::Install::Rhocall;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use Cwd;
use File::Spec::Functions qw{ catdir catfile };

## Cpanm
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_rhocall };
}

## Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub install_rhocall {

## Function : Install rhocall
## Returns  : ""
## Arguments: $program_parameters_href => Hash with rhocall specific parameters {REF}
##          : $conda_prefix_path       => Conda prefix path
##          : $conda_environment       => Conda environment
##          : $noupdate                => Do not update
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity
##          : $FILEHANDLE              => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $rhocall_parameters_href;
    my $conda_prefix_path;
    my $conda_environment;
    my $noupdate;
    my $quiet;
    my $verbose;
    my $FILEHANDLE;

    my $tmpl = {
        program_parameters_href => {
            required    => 1,
            default     => {},
            strict_type => 1,
            store       => \$rhocall_parameters_href
        },
        conda_prefix_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$conda_prefix_path
        },
        conda_environment => {
            strict_type => 1,
            store       => \$conda_environment
        },
        noupdate => {
            strict_type => 1,
            store       => \$noupdate
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$quiet
        },
        verbose => {
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$verbose
        },
        FILEHANDLE => {
            required => 1,
            defined  => 1,
            store    => \$FILEHANDLE
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Modules
    use MIP::Gnu::Coreutils qw{ gnu_rm gnu_mkdir };
    use MIP::Gnu::Bash qw{ gnu_cd };
    use MIP::Program::Download::Wget qw{ wget };
    use MIP::Program::Compression::Zip qw{ unzip };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };
    use MIP::Package_manager::Conda
      qw{ conda_source_activate conda_source_deactivate };
    use MIP::Package_manager::Pip qw{ pip_install };
    use MIP::Check::Installation qw{ check_existing_installation };

    ## Unpack parameters
    my $rhocall_version = $rhocall_parameters_href->{version};
    my $rhocall_path    = $rhocall_parameters_href->{path};

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => q{mip_install::install_rhocall},
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    ## Store original working directory
    my $pwd = cwd();

    say {$FILEHANDLE} q{### Install rhocall};

    ## Check if installation exists and remove directory unless a noupdate flag is provided
    my $install_check = check_existing_installation(
        {
            program_directory_path => $rhocall_path,
            program                => q{Rhocall},
            conda_environment      => $conda_environment,
            conda_prefix_path      => $conda_prefix_path,
            noupdate               => $noupdate,
            log                    => $log,
            FILEHANDLE             => $FILEHANDLE,
        }
    );

    # Return if the directory is found and a noupdate flag has been provided
    if ($install_check) {
        say {$FILEHANDLE} $NEWLINE;
        return;
    }

    if ($conda_environment) {
        ## Activate conda environment
        say $FILEHANDLE q{## Activate conda environment};
        conda_source_activate(
            {
                FILEHANDLE => $FILEHANDLE,
                env_name   => $conda_environment,
            }
        );
        say $FILEHANDLE $NEWLINE;
    }

    ## Creating install directory if not present already
    if ( not -d $rhocall_path ) {
        say {$FILEHANDLE} q{## Create rhocall install directory};
        gnu_mkdir(
            {
                indirectory_path => $rhocall_path,
                FILEHANDLE       => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Download
    say {$FILEHANDLE} q{## Download rhocall};
    my $url =
        q{https://github.com/dnil/rhocall/archive/}
      . $rhocall_version
      . $DOT . q{zip};
    my $rhocall_zip_path =
      catfile( $rhocall_path, q{rhocall-} . $rhocall_version . $DOT . q{zip} );
    wget(
        {
            url          => $url,
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $quiet,
            verbose      => $verbose,
            outfile_path => $rhocall_zip_path
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    unzip(
        {
            infile_path => $rhocall_zip_path,
            outdir_path => $rhocall_path,
            force       => 1,
            quiet       => $quiet,
            verbose     => $verbose,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Move to rhocall directory
    say {$FILEHANDLE} q{## Move to rhocall directory};
    my $rhocall_install_path =
      catdir( $rhocall_path, q{rhocall-} . $rhocall_version );
    gnu_cd(
        {
            directory_path => $rhocall_install_path,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Configure and install via pip
    say {$FILEHANDLE} q{## Configure and install};
    pip_install(
        {
            packages_ref => [qw{ numpy Cython }],
            quiet        => $quiet,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    print {$FILEHANDLE} $NEWLINE;
    pip_install(
        {
            requirement => q{requirements.txt},
            quiet       => $quiet,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    print {$FILEHANDLE} $NEWLINE;
    pip_install(
        {
            editable   => $DOT,
            quiet      => $quiet,
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Go back to subroutine origin
    say {$FILEHANDLE} q{## Moving back to original working directory};
    gnu_cd(
        {
            directory_path => $pwd,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Deactivate conda environment if conda_environment exists
    if ($conda_environment) {
        say $FILEHANDLE q{## Deactivate conda environment};
        conda_source_deactivate(
            {
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say $FILEHANDLE $NEWLINE;
    }

    print {$FILEHANDLE} $NEWLINE;

    return;
}

1;
