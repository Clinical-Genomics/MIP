package MIP::Recipes::Install::Tiddit;

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
    our @EXPORT_OK = qw{ install_tiddit };
}

## Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub install_tiddit {

## Function : Install TIDDIT
## Returns  : ""
## Arguments: $program_parameters_href => Hash with TIDDIT specific parameters {REF}
##          : $conda_prefix_path       => Conda prefix path
##          : $conda_environment       => Conda environment
##          : $noupdate                => Do not update
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity
##          : $FILEHANDLE              => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $tiddit_parameters_href;
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
            store       => \$tiddit_parameters_href
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
    use MIP::Gnu::Coreutils qw{ gnu_rm gnu_mkdir gnu_ln gnu_chmod };
    use MIP::Gnu::Bash qw{ gnu_cd };
    use MIP::Program::Download::Wget qw{ wget };
    use MIP::Program::Compression::Zip qw{ unzip };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };
    use MIP::Package_manager::Conda
      qw{ conda_source_activate conda_source_deactivate };
    use MIP::Gnu::Software::Gnu_make qw{ gnu_make };
    use MIP::Check::Installation qw{ check_existing_installation };
    use MIP::Compile::Cmake qw{ cmake };

    ## Unpack parameters
    my $tiddit_version = $tiddit_parameters_href->{version};

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => q{mip_install::install_tiddit},
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    ## Store original working directory
    my $pwd = cwd();

    say {$FILEHANDLE} q{### Install TIDDIT};

    ## Check if installation exists and remove directory unless a noupdate flag is provided
    my $tiddit_dir = catdir( $conda_prefix_path, q{TIDDIT-} . $tiddit_version );
    my $install_check = check_existing_installation(
        {
            program_directory_path => $tiddit_dir,
            program_name           => q{TIDDIT},
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

    ## Only activate conda environment if supplied by user
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

    ## Move to miniconda environment
    gnu_cd(
        {
            directory_path => $conda_prefix_path,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Download
    say {$FILEHANDLE} q{## Download TIDDIT};
    my $url =
        q{https://github.com/J35P312/TIDDIT/archive/}
      . $tiddit_version
      . $DOT . q{zip};
    my $tiddit_zip_path =
      catfile( q{TIDDIT-} . $tiddit_version . $DOT . q{zip} );
    wget(
        {
            url          => $url,
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $quiet,
            verbose      => $verbose,
            outfile_path => $tiddit_zip_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    unzip(
        {
            infile_path => $tiddit_zip_path,
            quiet       => $quiet,
            verbose     => $verbose,
            FILEHANDLE  => $FILEHANDLE,
            force       => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Move to TIDDIT directory
    say {$FILEHANDLE} q{## Move to TIDDIT directory};
    gnu_cd(
        {
            directory_path => q{TIDDIT-} . $tiddit_version,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    gnu_mkdir(
        {
            indirectory_path => q{build},
            parents          => 1,
            FILEHANDLE       => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    gnu_cd(
        {
            directory_path => q{build},
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Configure
    say {$FILEHANDLE} q{## Configure TIDDIT};
    cmake(
        {
            makefile_dir => q{..},
            FILEHANDLE   => $FILEHANDLE
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Compile
    say {$FILEHANDLE} q{## Compile};
    gnu_make(
        {
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Move to TIDDIT bin directory};
    gnu_cd(
        {
            directory_path => catdir(qw{ .. bin }),
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    gnu_chmod(
        {
            file_path  => q{TIDDIT},
            permission => q{a+x},
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Make available from conda environment
    say {$FILEHANDLE} q{## Make available from conda environment};
    my $target_path = catfile(
        $conda_prefix_path,
        q{TIDDIT-} . $tiddit_version,
        qw{ bin TIDDIT }
    );
    my $link_path = catfile( $conda_prefix_path, qw{ bin TIDDIT } );
    gnu_ln(
        {
            FILEHANDLE  => $FILEHANDLE,
            target_path => $target_path,
            link_path   => $link_path,
            symbolic    => 1,
            force       => 1,
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

    ## Clean-up
    gnu_rm(
        {
            infile_path => catfile( $conda_prefix_path, $tiddit_zip_path ),
            force       => 1,
            FILEHANDLE  => $FILEHANDLE,
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
