package MIP::Recipes::Install::Sambamba;

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
    our @EXPORT_OK = qw{ install_sambamba };
}

## Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub install_sambamba {

## Function : Install sambamba
## Returns  : ""
## Arguments: $program_parameters_href => Hash with sambamba specific parameters {REF}
##          : $conda_prefix_path       => Conda prefix path
##          : $conda_environment       => Conda environment
##          : $noupdate                => Do not update
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity
##          : $FILEHANDLE              => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sambamba_parameters_href;
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
            store       => \$sambamba_parameters_href
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
    use MIP::Gnu::Coreutils qw{ gnu_rm gnu_mv  gnu_ln gnu_chmod};
    use MIP::Program::Download::Wget qw{ wget };
    use MIP::Program::Compression::Bzip2 qw{ bzip2 };
    use MIP::Program::Compression::Tar qw{ tar };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };
    use MIP::Check::Installation qw{ check_existing_installation };
    use MIP::Script::Utils qw{ create_temp_dir };

    ## Unpack parameters
    my $sambamba_version = $sambamba_parameters_href->{version};

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => q{mip_install::install_sambamba},
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    ## Store original working directory
    my $pwd = cwd();

    say {$FILEHANDLE} q{### Install Sambamba};

    ## Check if installation exists and remove directory unless a noupdate flag is provided
    my $sambamba_dir = catdir( $conda_prefix_path, q{Sambamba} );
    my $install_check = check_existing_installation(
        {
            program_directory_path => $sambamba_dir,
            program_name           => q{Sambamba},
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
    say {$FILEHANDLE} q{## Create temporary sambamba install directory};
    my $temp_dir = create_temp_dir( { FILEHANDLE => $FILEHANDLE } );
    say {$FILEHANDLE} $NEWLINE;

    ## Download
    say {$FILEHANDLE} q{## Download Sambamba};
    my $url =
        q{https://github.com/lomereiter/sambamba/releases/download/v}
      . $sambamba_version
      . q{/sambamba_v}
      . $sambamba_version
      . $UNDERSCORE
      . q{linux.tar.bz2};
    my $sambamba_download_path = catfile( $temp_dir,
        q{sambamba_v} . $sambamba_version . $UNDERSCORE . q{linux.tar.bzp2} );
    wget(
        {
            url          => $url,
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $quiet,
            verbose      => $verbose,
            outfile_path => $sambamba_download_path
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Decompress
    say {$FILEHANDLE} q{## Decompress sambamba file};
    my $outfile_path = catfile( $temp_dir,
        q{sambamba_v} . $sambamba_version . $UNDERSCORE . q{linux.tar} );
    bzip2(
        {
            infile_path  => $sambamba_download_path,
            stdout       => 1,
            outfile_path => $outfile_path,
            decompress   => 1,
            force        => 1,
            quiet        => $quiet,
            verbose      => $verbose,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    $sambamba_download_path = catfile( $temp_dir,
        q{sambamba_v} . $sambamba_version . $UNDERSCORE . q{linux.tar} );
    tar(
        {
            file_path         => $sambamba_download_path,
            outdirectory_path => $temp_dir,
            extract           => 1,
            FILEHANDLE        => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Make executable
    say {$FILEHANDLE} q{## Make executable};
    gnu_chmod(
        {
            file_path =>
              catfile( $temp_dir, q{sambamba_v} . $sambamba_version ),
            permission => q{+x},
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Move directory to conda environment
    say {$FILEHANDLE} q{## Make available from conda environment};
    gnu_mv(
        {
            infile_path  => $temp_dir,
            outfile_path => $sambamba_dir,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Specifying target and link paths
    say {$FILEHANDLE} q{## Linking sambamba version};
    my $target_path =
      catfile( $sambamba_dir, q{sambamba_v} . $sambamba_version ),
      my $link_path = catfile( $conda_prefix_path, qw{ bin sambamba } );
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

    ## Remove compressed files
    say {$FILEHANDLE} q{## Remove compressed files};
    gnu_rm(
        {
            infile_path => catfile( $sambamba_dir, q{*tar*} ),
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE . $NEWLINE;

    return;
}

1;
