package MIP::Recipes::Install::Bedtools;

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
use MIP::Constants qw{ $ASTERISK $DOT $LOG $NEWLINE $SPACE $UNDERSCORE };
use MIP::Gnu::Coreutils qw{ gnu_ln gnu_mkdir gnu_mv gnu_rm };
use MIP::Gnu::Software::Gnu_make qw{ gnu_make };
use MIP::Check::Installation qw{ check_existing_installation };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Program::Download::Wget qw{ wget };
use MIP::Program::Compression::Tar qw{ tar };
use MIP::Script::Utils qw{ create_temp_dir };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_bedtools };
}

sub install_bedtools {

## Function : Install bedtools
## Returns  :

## Arguments: $conda_environment       => Conda environment
##          : $conda_prefix_path       => Conda prefix path
##          : $FILEHANDLE              => Filehandle to write to
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity
##          : $program_parameters_href => Hash with bedtools specific parameters {REF}
    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bedtools_parameters_href;
    my $conda_environment;
    my $conda_prefix_path;
    my $FILEHANDLE;
    my $quiet;
    my $verbose;

    my $tmpl = {
        conda_prefix_path => {
            defined     => 1,
            required    => 1,
            store       => \$conda_prefix_path,
            strict_type => 1,
        },
        conda_environment => {
            store       => \$conda_environment,
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
            store       => \$bedtools_parameters_href,
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
    my $bedtools_version      = $bedtools_parameters_href->{version};
    my $bedtools_main_version = substr $bedtools_version, 0, 1;

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

    say {$FILEHANDLE} q{### Install Bedtools};
    $log->info(qq{Writing instructions for Bedtools installation via SHELL});

    ## Check if installation exists and remove directory
    my $bedtools_dir = catdir( $conda_prefix_path, q{bedtools} . $bedtools_main_version );
    check_existing_installation(
        {
            conda_environment      => $conda_environment,
            conda_prefix_path      => $conda_prefix_path,
            FILEHANDLE             => $FILEHANDLE,
            log                    => $log,
            program_directory_path => $bedtools_dir,
            program_name           => q{Bedtools},
        }
    );

    ## Creating temporary install directory
    say {$FILEHANDLE} q{## Create temporary bedtools install directory};
    my $temp_dir = create_temp_dir( { FILEHANDLE => $FILEHANDLE } );
    say {$FILEHANDLE} $NEWLINE;

    ## Download
    say {$FILEHANDLE} q{## Download Bedtools};
    my $url =
        q{https://github.com/arq5x/bedtools}
      . $bedtools_main_version
      . q{/releases/download/v}
      . $bedtools_version
      . q{/bedtools-}
      . $bedtools_version
      . $DOT
      . q{tar.gz};
    my $bedtools_download_path =
      catfile( $temp_dir, q{bedtools-} . $bedtools_version . $DOT . q{tar.gz} );
    wget(
        {
            FILEHANDLE   => $FILEHANDLE,
            outfile_path => $bedtools_download_path,
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
            file_path         => $bedtools_download_path,
            FILEHANDLE        => $FILEHANDLE,
            outdirectory_path => $temp_dir,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Make
    say {$FILEHANDLE} q{## Recompile};
    my $makefile_dir = catdir( $temp_dir, q{bedtools} . $bedtools_main_version );
    gnu_make(
        {
            FILEHANDLE   => $FILEHANDLE,
            makefile_dir => $makefile_dir,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Move folder to conda_environment
    say {$FILEHANDLE} q{## Move to conda environment};
    gnu_mv(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $makefile_dir,
            outfile_path => $bedtools_dir,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Create links in conda environment bin to binaries in  bedtools folder
    say {$FILEHANDLE} q{## Create softlinks to binaries};
    gnu_ln(
        {
            FILEHANDLE  => $FILEHANDLE,
            force       => 1,
            link_path   => catfile( $conda_prefix_path, q{bin} ),
            symbolic    => 1,
            target_path => catfile( $bedtools_dir, q{bin}, $ASTERISK ),
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Remove the temporary install directory
    say {$FILEHANDLE} q{## Remove temporary install directory};
    gnu_rm(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $temp_dir,
            recursive   => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE x 2;

    return;
}

1;
