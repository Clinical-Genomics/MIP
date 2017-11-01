package MIP::Recipes::Install::Picard;

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
    our @EXPORT_OK = qw{ install_picard };
}

## Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub install_picard {

## Function : Install Picard
## Returns  : ""
## Arguments: $program_parameters_href => Hash with Picard specific parameters {REF}
##          : $conda_prefix_path       => Conda prefix path
##          : $conda_environment       => Conda environment
##          : $noupdate                => Do not update
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity
##          : $FILEHANDLE              => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $picard_parameters_href;
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
            store       => \$picard_parameters_href
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
    use MIP::Gnu::Coreutils qw{ gnu_rm gnu_mv gnu_mkdir gnu_ln };
    use MIP::Program::Download::Wget qw{ wget };
    use MIP::Program::Compression::Zip qw{ unzip };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };
    use MIP::Check::Installation qw{ check_existing_installation };

    ## Unpack parameters
    my $picard_version = $picard_parameters_href->{version};

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => q{mip_install::install_picard},
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    ## Store original working directory
    my $pwd = cwd();

    say {$FILEHANDLE} q{### Install Picard};

    ## Check if installation exists and remove directory unless a noupdate flag is provided
    my $picard_dir_path = catdir( $conda_prefix_path, q{share},
        q{picard-tools-} . $picard_version );
    my $install_check = check_existing_installation(
        {
            program_directory_path => $picard_dir_path,
            program                => q{Picard},
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

    ## Creating temporary install directory
    say {$FILEHANDLE} q{## Create temporary picard install directory};
    my $temp_dir = catdir( $pwd, q{picard_temp} . $UNDERSCORE . int rand 1000 );
    gnu_mkdir(
        {
            indirectory_path => $temp_dir,
            FILEHANDLE       => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Download
    say {$FILEHANDLE} q{## Download Picard};
    my $url =
        q{https://github.com/broadinstitute/picard/releases/download/}
      . $picard_version
      . q{/picard-tools-}
      . $picard_version
      . $DOT . q{zip};
    my $picard_zip_path =
      catfile( $temp_dir, q{picard-tools-} . $picard_version . $DOT . q{zip} );

    wget(
        {
            url          => $url,
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $quiet,
            verbose      => $verbose,
            outfile_path => $picard_zip_path
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    unzip(
        {
            infile_path => $picard_zip_path,
            outdir_path => $temp_dir,
            quiet       => $quiet,
            verbose     => $verbose,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Make available from conda environment
    say {$FILEHANDLE} q{## Make available from conda environment};
    gnu_mv(
        {
            infile_path =>
              catdir( $temp_dir, q{picard-tools-} . $picard_version ),
            outfile_path => catdir( $conda_prefix_path, q{share} ),
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Specifying target and link paths
    say {$FILEHANDLE} q{## Linking picard.jar};
    my $target_path =
      catfile( $conda_prefix_path, q{share}, q{picard-tools-} . $picard_version,
        q{picard.jar} );
    my $link_path = catfile( $conda_prefix_path, qw{ bin picard.jar } );
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

    ## Remove the temporary install directory
    gnu_rm(
        {
            infile_path => $temp_dir,
            recursive   => 1,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE . $NEWLINE;

    return;
}

1;
