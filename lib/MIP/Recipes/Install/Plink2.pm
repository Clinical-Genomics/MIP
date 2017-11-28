package MIP::Recipes::Install::Plink2;

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
    our @EXPORT_OK = qw{ install_plink2 };
}

## Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub install_plink2 {

## Function : Install Plink2
## Returns  : ""
## Arguments: $program_parameters_href => Hash with Plink2 specific parameters {REF}
##          : $conda_prefix_path       => Conda prefix path
##          : $conda_environment       => Conda environment
##          : $noupdate                => Do not update
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity
##          : $FILEHANDLE              => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $plink2_parameters_href;
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
            store       => \$plink2_parameters_href
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
    use MIP::Check::Installation qw{ check_existing_installation };
    use MIP::Gnu::Coreutils qw{ gnu_ln gnu_mv gnu_rm };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };
    use MIP::Program::Compression::Zip qw{ unzip };
    use MIP::Program::Download::Wget qw{ wget };
    use MIP::Script::Utils qw{ create_temp_dir };

    ## Unpack parameters
    my $plink2_version = $plink2_parameters_href->{version};

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => q{mip_install::install_plink2},
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    ## Store original working directory
    my $pwd = cwd();

    say {$FILEHANDLE} q{### Install Plink2};

    ## Check if installation exists and remove directory unless a noupdate flag is provided
    my $plink2_dir = catdir( $conda_prefix_path, q{Plink2} );
    my $install_check = check_existing_installation(
        {
            program_directory_path => $plink2_dir,
            program_name           => q{Plink2},
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
    say {$FILEHANDLE} q{## Create temporary Plink2 install directory};
    my $temp_dir = create_temp_dir( { FILEHANDLE => $FILEHANDLE } );
    say {$FILEHANDLE} $NEWLINE;

    ## Download
    say {$FILEHANDLE} q{## Download Plink2};
    my $url =
        q{https://www.cog-genomics.org/static/bin/plink}
      . $plink2_version
      . q{/plink_linux_x86_64.zip};
    my $plink2_zip_path = catfile( $temp_dir, q{plink_linux_x86_64.zip} );
    wget(
        {
            url          => $url,
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $quiet,
            verbose      => $verbose,
            outfile_path => $plink2_zip_path
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    unzip(
        {
            infile_path => $plink2_zip_path,
            outdir_path => $temp_dir,
            quiet       => $quiet,
            verbose     => $verbose,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Remove the zip file
    say {$FILEHANDLE} q{## Remove zip file};
    gnu_rm(
        {
            infile_path => $plink2_zip_path,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Make available from conda environment
    say {$FILEHANDLE} q{## Make available from conda environment};
    gnu_mv(
        {
            infile_path  => $temp_dir,
            outfile_path => $plink2_dir,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Create link in conda's bin directory to binary
    say {$FILEHANDLE} q{## Create softlink to binary};
    gnu_ln(
        {
            link_path   => catfile( $conda_prefix_path, q{bin} ),
            target_path => catfile( $plink2_dir,        q{plink2} ),
            symbolic    => 1,
            force       => 1,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE . $NEWLINE;

    return;
}

1;
