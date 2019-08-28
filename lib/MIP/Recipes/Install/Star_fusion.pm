package MIP::Recipes::Install::Star_fusion;

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
use MIP::Constants qw{ $ASTERISK $DOT $LOG $NEWLINE $SPACE $UNDERSCORE };
use MIP::Gnu::Coreutils qw{ gnu_ln gnu_rm };
use MIP::Gnu::Software::Gnu_make qw{ gnu_make };
use MIP::Gnu::Software::Gnu_sed qw{ gnu_sed };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Program::Compression::Tar qw{ tar };
use MIP::Program::Download::Wget qw{ wget };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_star_fusion };
}

sub install_star_fusion {

## Function : Install STAR-Fusion
## Returns  :
## Arguments: $conda_environment       => Conda environment
##          : $conda_prefix_path       => Conda prefix path
##          : $FILEHANDLE              => Filehandle to write to
##          : $program_parameters_href => Hash with star_fusion specific parameters {REF}
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment;
    my $conda_prefix_path;
    my $FILEHANDLE;
    my $quiet;
    my $star_fusion_parameters_href;
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
            store       => \$star_fusion_parameters_href,
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
    my $star_fusion_version = $star_fusion_parameters_href->{version};

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => $LOG,
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    ## Store original working directory
    my $pwd = cwd();

    say {$FILEHANDLE} q{### Install Star-Fusion};
    $log->info(qq{Writing instructions for Star-Fusion installation via SHELL});

    ## Check if installation exists and remove directory
    my $star_fusion_dir =
      catdir( $conda_prefix_path, q{share}, q{STAR-Fusion-v} . $star_fusion_version );
    check_existing_installation(
        {
            conda_environment      => $conda_environment,
            conda_prefix_path      => $conda_prefix_path,
            FILEHANDLE             => $FILEHANDLE,
            log                    => $log,
            program_directory_path => $star_fusion_dir,
            program_name           => q{Star-Fusion},
        }
    );

    ## Download
    say {$FILEHANDLE} q{## Download Star-Fusion};
    my $url =
        q{https://github.com/STAR-Fusion/STAR-Fusion/releases/download/STAR-Fusion-v}
      . $star_fusion_version
      . q{/STAR-Fusion-v}
      . $star_fusion_version
      . $DOT
      . q{FULL.tar.gz};
    my $star_fusion_download_path =
      catfile( $conda_prefix_path, q{share},
        q{STAR-Fusion-v} . $star_fusion_version . $DOT . q{FULL.tar.gz} );
    wget(
        {
            FILEHANDLE   => $FILEHANDLE,
            outfile_path => $star_fusion_download_path,
            quiet        => $quiet,
            url          => $url,
            verbose      => $verbose,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    tar(
        {
            extract           => 1,
            file_path         => $star_fusion_download_path,
            FILEHANDLE        => $FILEHANDLE,
            outdirectory_path => catdir( $conda_prefix_path, q{share} ),
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Make
    say {$FILEHANDLE} q{## Recompile};
    gnu_make(
        {
            makefile_dir => $star_fusion_dir,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Replace occurrences of $FindBin::Bin with $FindBin::RealBin in STAR-Fusion script
    ## Done so that Star-Fusion can be executed from the link in the conda's bin folder.
    my $sed_script = q{'s/\$FindBin::Bin/\$FindBin::RealBin/g'};
    say {$FILEHANDLE} q{## Find and replace Bin with RealBin};
    gnu_sed(
        {
            FILEHANDLE   => $FILEHANDLE,
            script       => $sed_script,
            infile_path  => catfile( $star_fusion_dir, q{STAR-Fusion} ),
            inplace_edit => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Create link in conda environment bin to binary in Star-Fusion folder
    say {$FILEHANDLE} q{## Create softlinks to binary};
    gnu_ln(
        {
            link_path   => catfile( $conda_prefix_path, q{bin} ),
            target_path => catfile( $star_fusion_dir,   q{STAR-Fusion} ),
            symbolic    => 1,
            force       => 1,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Remove the downloaded tar file
    say {$FILEHANDLE} q{## Remove tar file};
    gnu_rm(
        {
            infile_path => $star_fusion_download_path,
            recursive   => 1,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE x 2;

    return;
}

1;
