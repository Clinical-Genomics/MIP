package MIP::Recipes::Install::Vt;

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
use MIP::Check::Installation qw{ check_existing_installation };
use MIP::Constants qw{ $DASH $DOT $LOG_NAME $NEWLINE $SPACE $UNDERSCORE };
use MIP::Gnu::Coreutils qw{ gnu_ln gnu_mkdir gnu_mv gnu_rm };
use MIP::Gnu::Software::Gnu_make qw{ gnu_make };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Script::Utils qw{ create_temp_dir };
use MIP::Versionmanager::Git qw{ git_clone };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_vt };
}

sub install_vt {

## Function : Install vt
## Returns  :
## Arguments: $conda_environment       => Conda environment
##          : $conda_prefix_path       => Conda prefix path
##          : $FILEHANDLE              => Filehandle to write to
##          : $program_parameters_href => Hash with vt specific parameters {REF}
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment;
    my $conda_prefix_path;
    my $FILEHANDLE;
    my $quiet;
    my $verbose;
    my $vt_parameters_href;

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
        program_parameters_href => {
            default     => {},
            required    => 1,
            store       => \$vt_parameters_href,
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

    ## Unpack parameters
    my $vt_version = $vt_parameters_href->{version};

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => $LOG_NAME,
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    ## Store original working directory
    my $pwd = cwd();

    say {$FILEHANDLE} q{### Install VT};
    $log->info(qq{Writing instructions for VT installation via SHELL});

    ## Check if installation exists and remove directory
    my $vt_dir = catdir( $conda_prefix_path, q{Vt} );
    check_existing_installation(
        {
            conda_environment      => $conda_environment,
            conda_prefix_path      => $conda_prefix_path,
            FILEHANDLE             => $FILEHANDLE,
            log                    => $log,
            program_directory_path => $vt_dir,
            program_name           => q{VT},
        }
    );

    ## Creating temporary install directory
    say {$FILEHANDLE} q{## Create temporary VT install directory};
    my $temp_dir = create_temp_dir( { FILEHANDLE => $FILEHANDLE } );
    say {$FILEHANDLE} $NEWLINE;

    ## Download
    say {$FILEHANDLE} q{## Download VT};
    git_clone(
        {
            FILEHANDLE  => $FILEHANDLE,
            outdir_path => $temp_dir,
            quiet       => $quiet,
            url         => q{https://github.com/atks/vt.git},
            verbose     => $verbose,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Make
    say {$FILEHANDLE} q{## Recompile};
    gnu_make(
        {
            FILEHANDLE   => $FILEHANDLE,
            makefile_dir => $temp_dir,
        }
    );
    print {$FILEHANDLE} $NEWLINE;
    gnu_make(
        {
            FILEHANDLE   => $FILEHANDLE,
            makefile_dir => $temp_dir,
            test         => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Create conda directory
    say {$FILEHANDLE} q{## Make VT directory in conda env};
    gnu_mkdir(
        {
            FILEHANDLE       => $FILEHANDLE,
            indirectory_path => $vt_dir,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Move to Conda env
    say {$FILEHANDLE} q{## Make available from conda environment};
    gnu_mv(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => catfile( $temp_dir, q{vt} ),
            outfile_path => $vt_dir,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Create link in conda's bin directory to binary
    say {$FILEHANDLE} q{## Create softlink to binary};
    gnu_ln(
        {
            FILEHANDLE  => $FILEHANDLE,
            force       => 1,
            link_path   => catfile( $conda_prefix_path, q{bin} ),
            symbolic    => 1,
            target_path => catfile( $vt_dir, q{vt} ),
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Remove temporary install directory};
    ## Remove the temporary install directory
    gnu_rm(
        {
            FILEHANDLE  => $FILEHANDLE,
            force       => 1,
            infile_path => $temp_dir,
            recursive   => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    return;
}
1;
