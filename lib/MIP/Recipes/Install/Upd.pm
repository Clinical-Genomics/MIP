package MIP::Recipes::Install::Upd;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## CPAN
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $DASH $DOT $LOG $NEWLINE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_upd };
}

sub install_upd {

## Function : Install UPD
## Returns  :
## Arguments: $conda_environment       => Conda environment
##          : $conda_prefix_path       => Conda prefix path
##          : $FILEHANDLE              => Filehandle to write to
##          : $noupdate                => Do not update
##          : $program_parameters_href => Hash with Program specific parameters {REF}
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment;
    my $conda_prefix_path;
    my $FILEHANDLE;
    my $noupdate;
    my $upd_parameters_href;
    my $quiet;
    my $verbose;

    my $tmpl = {
        conda_environment => {
            store       => \$conda_environment,
            strict_type => 1,
        },
        conda_prefix_path => {
            defined     => 1,
            required    => 1,
            store       => \$conda_prefix_path,
            strict_type => 1,
        },
        FILEHANDLE => {
            defined  => 1,
            required => 1,
            store    => \$FILEHANDLE,
        },
        noupdate => {
            store       => \$noupdate,
            strict_type => 1,
        },
        program_parameters_href => {
            default     => {},
            required    => 1,
            store       => \$upd_parameters_href,
            strict_type => 1,
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

    ## Modules
    use MIP::Check::Installation qw{ check_existing_installation };
    use MIP::Gnu::Bash qw{ gnu_cd };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };
    use MIP::Package_manager::Conda qw{ conda_activate conda_deactivate };
    use MIP::Package_manager::Pip qw{ pip_install };
    use MIP::Program::Download::Wget qw{ wget };
    use MIP::Program::Compression::Zip qw{ unzip };

    ## Unpack parameters
    my $program_version = $upd_parameters_href->{version};

    ## Set program specific parameters
    my $program_name = q{UPD};
    my $executable   = q{upd};
    my $program_directory_path =
      catdir( $conda_prefix_path, q{share}, $executable . $DASH . $program_version );
    my $program_url =
      q{https://github.com/bjhall/upd/archive/} . $program_version . $DOT . q{zip};

    ## Store original working directory
    my $pwd = cwd();

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => $LOG,
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    say {$FILEHANDLE} q{### Install} . $SPACE . $program_name;

    ## Check if installation exists and remove directory unless a noupdate flag is provided
    my $install_check = check_existing_installation(
        {
            conda_environment      => $conda_environment,
            conda_prefix_path      => $conda_prefix_path,
            FILEHANDLE             => $FILEHANDLE,
            log                    => $log,
            noupdate               => $noupdate,
            program_directory_path => $program_directory_path,
            program_name           => $program_name,
        }
    );

    ## Return if the directory is found and a noupdate flag has been provided
    if ($install_check) {
        say {$FILEHANDLE} $NEWLINE;
        return;
    }

    ## Only activate conda environment if supplied by user
    if ($conda_environment) {
        ## Activate conda environment
        say {$FILEHANDLE} q{## Activate conda environment};
        conda_activate(
            {
                env_name   => $conda_environment,
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Download
    say {$FILEHANDLE} q{## Download} . $SPACE . $program_name;
    my $program_zip_path =
      catfile( $conda_prefix_path,
        $program_name . $DASH . $program_version . $DOT . q{zip} );
    wget(
        {
            FILEHANDLE   => $FILEHANDLE,
            outfile_path => $program_zip_path,
            quiet        => $quiet,
            url          => $program_url,
            verbose      => $verbose,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    unzip(
        {
            FILEHANDLE  => $FILEHANDLE,
            force       => 1,
            infile_path => $program_zip_path,
            outdir_path => catdir( $conda_prefix_path, q{share} ),
            quiet       => $quiet,
            verbose     => $verbose,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Move to UPD directory
    say {$FILEHANDLE} q{## Move to} . $SPACE . $program_name . $SPACE . q{directory};
    gnu_cd(
        {
            directory_path => $program_directory_path,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Install
    pip_install(
        {
            editable      => $DOT,
            FILEHANDLE    => $FILEHANDLE,
            python_module => 1,
            quiet         => $quiet,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Go back to starting directoriy
    say {$FILEHANDLE} q{## Go back to starting directory};
    gnu_cd(
        {
            directory_path => $pwd,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Deactivate conda environment if conda_environment exists
    if ($conda_environment) {
        say {$FILEHANDLE} q{## Deactivate conda environment};
        conda_deactivate(
            {
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    print {$FILEHANDLE} $NEWLINE;
    return;
}

1;
