package MIP::Recipes::Install::Cnvnator;

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
    our @EXPORT_OK = qw{ install_cnvnator };
}

## Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub install_cnvnator {

## Function : Install CNVnator
## Returns  : ""
## Arguments: $program_parameters_href => Hash with CNVnator specific parameters {REF}
##          : $conda_prefix_path       => Conda prefix path
##          : $conda_environment       => Conda environment
##          : $noupdate                => Do not update
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity
##          : $FILEHANDLE              => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cnvnator_parameters_href;
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
            store       => \$cnvnator_parameters_href
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
    use MIP::Gnu::Coreutils qw{ gnu_rm gnu_mkdir gnu_mv gnu_ln};
    use MIP::Gnu::Bash qw{ gnu_cd };
    use MIP::Program::Download::Wget qw{ wget };
    use MIP::Program::Compression::Tar qw{ tar };
    use MIP::Program::Compression::Zip qw{ unzip };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };
    use MIP::Package_manager::Conda
      qw{ conda_source_activate conda_source_deactivate };
    use MIP::Gnu::Software::Gnu_make qw{ gnu_make };
    use MIP::Check::Installation qw{ check_existing_installation };

    ## Unpack parameters
    my $cnvnator_version = $cnvnator_parameters_href->{version};
    my $cnvnator_root_binary =
      $cnvnator_parameters_href->{cnvnator_root_binary};

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => q{mip_install::install_cnvnator},
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    ## Store original working directory
    my $pwd = cwd();

    say {$FILEHANDLE} q{### Install CNVnator/ROOT};

    ## ROOT installation (prerequisite for CNVnator)

    ## Check if installation exists and remove directory unless a noupdate flag is provided
    my $root_bin_dir = catdir( $conda_prefix_path, q{root} );
    my $root_install_check = check_existing_installation(
        {
            program_directory_path => $root_bin_dir,
            program                => q{ROOT (CNVnator prequisite)},
            conda_environment      => $conda_environment,
            conda_prefix_path      => $conda_prefix_path,
            noupdate               => $noupdate,
            log                    => $log,
            FILEHANDLE             => $FILEHANDLE,
        }
    );

    # Return if the directory is found and a noupdate flag has been provided
    if ($root_install_check) {
        goto CNVNATOR;
    }

    ## Move to miniconda environment
    gnu_cd(
        {
            directory_path => $conda_prefix_path,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say $FILEHANDLE $NEWLINE;

    ## Download ROOT
    say {$FILEHANDLE} q{## Download ROOT};
    my $root_url = q{https://root.cern.ch/download/} . $cnvnator_root_binary;
    wget(
        {
            url          => $root_url,
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $quiet,
            verbose      => $verbose,
            outfile_path => $cnvnator_root_binary,
        }
    );

    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    tar(
        {
            file_path  => $cnvnator_root_binary,
            extract    => 1,
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    my $binary_regexp = catdir( $conda_prefix_path, qw{ root bin } );
    if ( not $ENV{PATH} =~ m{$binary_regexp}xms ) {

        ## Export path
        say {$FILEHANDLE} q{## Export path};
        say {$FILEHANDLE}
          q{echo '# Added by mip_installer ' "$(date)" >> ~/.bashrc};
        say {$FILEHANDLE} q{echo 'source}
          . $SPACE
          . catfile( $conda_prefix_path, qw{ root bin thisroot.sh } )
          . q{' >> ~/.bashrc}
          . $NEWLINE;

        ## Use newly installed root
        say {$FILEHANDLE} q{source}
          . $SPACE
          . catfile( $conda_prefix_path, qw{ root bin thisroot.sh } )
          . $NEWLINE;
    }

    ## Remove ROOT archive
    say {$FILEHANDLE} q{## Removing ROOT archive};
    gnu_rm(
        {
            infile_path => $cnvnator_root_binary,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Moving back to original working directory};
    gnu_cd(
        {
            directory_path => $pwd,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

  CNVNATOR:

    ## Check if installation exists and remove directory unless a noupdate flag is provided
    my $cnvnator_bin_dir = catdir( $conda_prefix_path, q{CVNnator} );
    my $cnvnator_install_check = check_existing_installation(
        {
            program_directory_path => $cnvnator_bin_dir,
            program                => q{CNVnator},
            conda_environment      => $conda_environment,
            conda_prefix_path      => $conda_prefix_path,
            noupdate               => $noupdate,
            log                    => $log,
            FILEHANDLE             => $FILEHANDLE,
        }
    );

    # Return if the directory is found and a noupdate flag has been provided
    if ($cnvnator_install_check) {
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

    ## Creating temporary install directory
    say {$FILEHANDLE} q{## Create temporary CNVnator install directory};
    my $temp_dir =
      catdir( $pwd, q{cnvnator_temp} . $UNDERSCORE . int rand 1000 );
    gnu_mkdir(
        {
            indirectory_path => $temp_dir,
            FILEHANDLE       => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Download
    say {$FILEHANDLE} q{## Download CNVNator};
    my $cnvnator_url =
        q{https://github.com/abyzovlab/CNVnator/releases/download/v}
      . $cnvnator_version
      . q{/CNVnator_v}
      . $cnvnator_version
      . $DOT . q{zip};
    my $cnvnator_zip_path =
      catfile( $temp_dir, q{CNVnator_v} . $cnvnator_version . $DOT . q{zip} );
    wget(
        {
            url          => $cnvnator_url,
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $quiet,
            verbose      => $verbose,
            outfile_path => $cnvnator_zip_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    unzip(
        {
            infile_path => $cnvnator_zip_path,
            outdir_path => $temp_dir,
            quiet       => $quiet,
            verbose     => $verbose,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Move to CNVnator directory
    say {$FILEHANDLE} q{## Move to CNVnator directory};
    my $cnvnator_configure_path =
      catdir( $temp_dir, q{CNVnator_v} . $cnvnator_version,
        qw{ src samtools } );
    gnu_cd(
        {
            directory_path => $cnvnator_configure_path,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Configure
    say {$FILEHANDLE} q{## Configure CNVnator samTools specific version};
    say {$FILEHANDLE} q{./configure --without-curses} . $NEWLINE;

    ## Compile
    say {$FILEHANDLE} q{## Compile};
    gnu_make(
        {
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Move to CNVnator directory};
    gnu_cd(
        {
            directory_path => q{..},
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    gnu_make(
        {
            FILEHANDLE => $FILEHANDLE,
        }
    );

    # Add no parallel support argument to make command
    say {$FILEHANDLE} q{OMP=no} . $NEWLINE;

    ## Make available from conda environment
    say {$FILEHANDLE} q{## Move to conda environment};
    my $cnvnator_path = catdir( $temp_dir, q{CNVnator_v} . $cnvnator_version );
    gnu_mv(
        {
            infile_path  => $cnvnator_path,
            outfile_path => $cnvnator_bin_dir,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Link binaries to conda bin directory
    say {$FILEHANDLE} q{## Link binaries to conda bin directory};
    gnu_ln(
        {
            link_path   => catfile( $conda_prefix_path, q{bin} ),
            target_path => catfile( $cnvnator_bin_dir,  qw{ src cnvnator } ),
            symbolic    => 1,
            force       => 1,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    print {$FILEHANDLE} $NEWLINE;
    gnu_ln(
        {
            link_path   => catfile( $conda_prefix_path, q{bin} ),
            target_path => catfile( $cnvnator_bin_dir,  q{cnvnator2VCF.pl} ),
            symbolic    => 1,
            force       => 1,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Go back to starting directory
    say {$FILEHANDLE} q{## Remove temporary install directory};
    gnu_cd(
        {
            directory_path => $pwd,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Remove the temporary install directory
    gnu_rm(
        {
            infile_path => $temp_dir,
            FILEHANDLE  => $FILEHANDLE,
            recursive   => 1,
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
