package MIP::Recipes::Install::Sambamba;

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
use MIP::Constants qw{ $LOG_NAME $NEWLINE $SPACE $UNDERSCORE };
use MIP::Gnu::Coreutils qw{ gnu_rm gnu_mv  gnu_ln gnu_chmod};
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Program::Compression::Bzip2 qw{ bzip2 };
use MIP::Program::Compression::Tar qw{ tar };
use MIP::Program::Download::Wget qw{ wget };
use MIP::Script::Utils qw{ create_temp_dir };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_sambamba };
}

sub install_sambamba {

## Function : Install sambamba
## Returns  :
## Arguments: $conda_environment       => Conda environment
##          : $conda_prefix_path       => Conda prefix path
##          : $FILEHANDLE              => Filehandle to write to
##          : $program_parameters_href => Hash with sambamba specific parameters {REF}
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment;
    my $conda_prefix_path;
    my $FILEHANDLE;
    my $quiet;
    my $sambamba_parameters_href;
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
        program_parameters_href => {
            default     => {},
            required    => 1,
            store       => \$sambamba_parameters_href,
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
    my $sambamba_version = $sambamba_parameters_href->{version};

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

    say {$FILEHANDLE} q{### Install Sambamba};
    $log->info(qq{Writing instructions for Sambamba installation via SHELL});

    ## Check if installation exists and remove directory
    my $sambamba_dir = catdir( $conda_prefix_path, q{Sambamba} );
    check_existing_installation(
        {
            conda_environment      => $conda_environment,
            conda_prefix_path      => $conda_prefix_path,
            FILEHANDLE             => $FILEHANDLE,
            log                    => $log,
            program_directory_path => $sambamba_dir,
            program_name           => q{Sambamba},
        }
    );

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
            FILEHANDLE   => $FILEHANDLE,
            outfile_path => $sambamba_download_path,
            quiet        => $quiet,
            url          => $url,
            verbose      => $verbose,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Decompress
    say {$FILEHANDLE} q{## Decompress sambamba file};
    my $outfile_path = catfile( $temp_dir,
        q{sambamba_v} . $sambamba_version . $UNDERSCORE . q{linux.tar} );
    bzip2(
        {
            decompress   => 1,
            FILEHANDLE   => $FILEHANDLE,
            force        => 1,
            infile_path  => $sambamba_download_path,
            outfile_path => $outfile_path,
            quiet        => $quiet,
            stdout       => 1,
            verbose      => $verbose,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    $sambamba_download_path = catfile( $temp_dir,
        q{sambamba_v} . $sambamba_version . $UNDERSCORE . q{linux.tar} );
    tar(
        {
            extract           => 1,
            file_path         => $sambamba_download_path,
            FILEHANDLE        => $FILEHANDLE,
            outdirectory_path => $temp_dir,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Make executable
    say {$FILEHANDLE} q{## Make executable};
    gnu_chmod(
        {
            file_path  => catfile( $temp_dir, q{sambamba_v} . $sambamba_version ),
            FILEHANDLE => $FILEHANDLE,
            permission => q{+x},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Move directory to conda environment
    say {$FILEHANDLE} q{## Make available from conda environment};
    gnu_mv(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $temp_dir,
            outfile_path => $sambamba_dir,
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
            force       => 1,
            link_path   => $link_path,
            symbolic    => 1,
            target_path => $target_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Remove compressed files
    say {$FILEHANDLE} q{## Remove compressed files};
    gnu_rm(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => catfile( $sambamba_dir, q{*tar*} ),
        }
    );
    say {$FILEHANDLE} $NEWLINE . $NEWLINE;

    return;
}

1;
