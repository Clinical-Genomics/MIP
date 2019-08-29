package MIP::Recipes::Install::Cnvnator;

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
use MIP::Constants qw{ $DOT $LOG $NEWLINE $SPACE $UNDERSCORE };
use MIP::Gnu::Bash qw{ gnu_cd };
use MIP::Gnu::Coreutils qw{ gnu_ln gnu_mv gnu_rm};
use MIP::Gnu::Software::Gnu_make qw{ gnu_make };
use MIP::Gnu::Software::Gnu_sed qw{ gnu_sed };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Package_manager::Conda qw{ conda_activate conda_deactivate };
use MIP::Program::Compression::Tar qw{ tar };
use MIP::Program::Compression::Zip qw{ unzip };
use MIP::Program::Download::Wget qw{ wget };
use MIP::Recipes::Install::Root qw{ install_root };
use MIP::Script::Utils qw{ create_temp_dir };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_cnvnator };
}

sub install_cnvnator {

## Function : Install CNVnator
## Returns  :
## Arguments: $conda_environment       => Conda environment
##          : $conda_prefix_path       => Conda prefix path
##          : $FILEHANDLE              => Filehandle to write to
##          : $program_parameters_href => Hash with CNVnator specific parameters {REF}
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cnvnator_parameters_href;
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
            store       => \$cnvnator_parameters_href,
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
    my $cnvnator_version     = $cnvnator_parameters_href->{version};
    my $cnvnator_root_binary = $cnvnator_parameters_href->{cnvnator_root_binary};

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

    say {$FILEHANDLE} q{### Install CNVnator/ROOT};
    $log->info(qq{Writing instructions for CNVnator/ROOT installation via SHELL});

    ## Check if ROOT (CNVnator requirement) installation exists and remove directory
    my $root_bin_dir = catdir( $conda_prefix_path, q{root} );
    check_existing_installation(
        {
            conda_environment      => $conda_environment,
            conda_prefix_path      => $conda_prefix_path,
            FILEHANDLE             => $FILEHANDLE,
            log                    => $log,
            program_directory_path => $root_bin_dir,
            program_name           => q{ROOT (CNVnator prequisite)},
        }
    );

    ## Install Root
    install_root(
        {
            conda_prefix_path => $conda_prefix_path,
            FILEHANDLE        => $FILEHANDLE,
            quiet             => $quiet,
            root_binary       => $cnvnator_root_binary,
            verbose           => $verbose,
        }
    );

    say {$FILEHANDLE} q{## Install CNVnator};

    ## Check if CNVnator installation exists and remove directory
    my $cnvnator_bin_dir = catdir( $conda_prefix_path, q{CNVnator} );
    check_existing_installation(
        {
            conda_environment      => $conda_environment,
            conda_prefix_path      => $conda_prefix_path,
            FILEHANDLE             => $FILEHANDLE,
            log                    => $log,
            program_directory_path => $cnvnator_bin_dir,
            program_name           => q{CNVnator},
        }
    );

    ## Activate conda environment
    say {$FILEHANDLE} q{## Activate conda environment};
    conda_activate(
        {
            env_name   => $conda_environment,
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Creating temporary install directory
    say {$FILEHANDLE} q{## Create temporary CNVnator install directory};
    my $temp_dir_path = create_temp_dir( { FILEHANDLE => $FILEHANDLE } );
    say {$FILEHANDLE} $NEWLINE;

    ## Download
    say {$FILEHANDLE} q{## Download CNVnator};
    my $cnvnator_url =
        q{https://github.com/abyzovlab/CNVnator/releases/download/v}
      . $cnvnator_version
      . q{/CNVnator_v}
      . $cnvnator_version
      . $DOT . q{zip};
    my $cnvnator_zip_path =
      catfile( $temp_dir_path, q{CNVnator_v} . $cnvnator_version . $DOT . q{zip} );
    wget(
        {
            FILEHANDLE   => $FILEHANDLE,
            outfile_path => $cnvnator_zip_path,
            quiet        => $quiet,
            url          => $cnvnator_url,
            verbose      => $verbose,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    unzip(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $cnvnator_zip_path,
            outdir_path => $temp_dir_path,
            quiet       => $quiet,
            verbose     => $verbose,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Move to CNVnator directory
    say {$FILEHANDLE} q{## Move to CNVnator directory};
    my $cnvnator_configure_path =
      catdir( $temp_dir_path, q{CNVnator_v} . $cnvnator_version, qw{ src samtools } );
    gnu_cd(
        {
            directory_path => $cnvnator_configure_path,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Configure
    say {$FILEHANDLE} q{## Configure CNVnator samtools specific version};
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

    ## Modify Makefile according to https://github.com/abyzovlab/CNVnator/issues/15#issuecomment-370376682
    say {$FILEHANDLE} q{## Modify CNVnator makefile};
    my $sed_script = q{'s/-std=c++11/-std=c++11 -lpthread/g'};
    gnu_sed(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => q{Makefile},
            inplace_edit => 1,
            script       => $sed_script,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Get make command
    my @make_commands = gnu_make( {} );
    ## Add no parallel support argument to make command
    push @make_commands, q{OMP=no};
    unix_write_to_file(
        {
            commands_ref => \@make_commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Make available from conda environment
    say {$FILEHANDLE} q{## Move to conda environment};
    my $cnvnator_path = catdir( $temp_dir_path, q{CNVnator_v} . $cnvnator_version );
    gnu_mv(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $cnvnator_path,
            outfile_path => $cnvnator_bin_dir,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Link binaries to conda bin directory
    say {$FILEHANDLE} q{## Link binaries to conda bin directory};
    my @executable_paths = (
        catfile( $cnvnator_bin_dir, qw{ src cnvnator } ),
        catfile( $cnvnator_bin_dir, q{cnvnator2VCF.pl} )
    );
  EXECUTABLE_PATH:
    foreach my $executable_path (@executable_paths) {
        gnu_ln(
            {
                FILEHANDLE  => $FILEHANDLE,
                force       => 1,
                link_path   => catfile( $conda_prefix_path, q{bin} ),
                symbolic    => 1,
                target_path => $executable_path,
            }
        );
        print {$FILEHANDLE} $NEWLINE;
    }
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
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $temp_dir_path,
            recursive   => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Deactivate conda environment if conda_environment exists
    say {$FILEHANDLE} q{## Deactivate conda environment};
    conda_deactivate(
        {
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    return;
}

1;
