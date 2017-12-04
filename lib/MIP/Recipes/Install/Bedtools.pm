package MIP::Recipes::Install::Bedtools;

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
    our @EXPORT_OK = qw{ install_bedtools };
}

## Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};
Readonly my $ASTERIX    => q{*};

sub install_bedtools {

## Function : Install bedtools
## Returns  :
## Arguments: $program_parameters_href => Hash with bedtools specific parameters {REF}
##          : $conda_prefix_path       => Conda prefix path
##          : $conda_environment       => Conda environment
##          : $noupdate                => Do not update
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity
##          : $FILEHANDLE              => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bedtools_parameters_href;
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
            store       => \$bedtools_parameters_href
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
    use MIP::Gnu::Coreutils qw{ gnu_ln gnu_mkdir gnu_mv gnu_rm };
    use MIP::Gnu::Software::Gnu_make qw{ gnu_make };
    use MIP::Check::Installation qw{ check_existing_installation };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };
    use MIP::Program::Download::Wget qw{ wget };
    use MIP::Program::Compression::Tar qw{ tar };
    use MIP::Script::Utils qw{ create_temp_dir };

    ## Unpack parameters
    my $bedtools_version = $bedtools_parameters_href->{version};
    my $bedtools_main_version = substr $bedtools_version, 0, 1;

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => q{mip_install::install_bedtools},
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    ## Store original working directory
    my $pwd = cwd();

    say {$FILEHANDLE} q{### Install Bedtools};

    ## Check if installation exists and remove directory unless a noupdate flag is provided
    my $bedtools_dir =
      catdir( $conda_prefix_path, q{bedtools} . $bedtools_main_version );
    my $install_check = check_existing_installation(
        {
            program_directory_path => $bedtools_dir,
            program_name           => q{Bedtools},
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
            url          => $url,
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $quiet,
            verbose      => $verbose,
            outfile_path => $bedtools_download_path
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    tar(
        {
            file_path         => $bedtools_download_path,
            outdirectory_path => $temp_dir,
            extract           => 1,
            FILEHANDLE        => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Make
    say {$FILEHANDLE} q{## Recompile};
    my $makefile_dir =
      catdir( $temp_dir, q{bedtools} . $bedtools_main_version );
    gnu_make(
        {
            makefile_dir => $makefile_dir,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Move folder to conda_environment
    say {$FILEHANDLE} q{## Move to conda environment};
    gnu_mv(
        {
            infile_path  => $makefile_dir,
            outfile_path => $bedtools_dir,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Create links in conda environment bin to binaries in  bedtools folder
    say {$FILEHANDLE} q{## Create softlinks to binaries};
    gnu_ln(
        {
            link_path   => catfile( $conda_prefix_path, q{bin} ),
            target_path => catfile( $bedtools_dir,      q{bin}, $ASTERIX ),
            symbolic   => 1,
            force      => 1,
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Remove the temporary install directory
    say {$FILEHANDLE} q{## Remove temporary install directory};
    gnu_rm(
        {
            infile_path => $temp_dir,
            recursive   => 1,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE x 2;

    return;
}

1;
