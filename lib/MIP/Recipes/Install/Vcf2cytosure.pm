package MIP::Recipes::Install::Vcf2cytosure;

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
use MIP::Constants qw{ $DASH $DOT $LOG $NEWLINE $SPACE };
use MIP::Gnu::Bash qw{ gnu_cd };
use MIP::Gnu::Coreutils qw{ gnu_rm };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Package_manager::Conda qw{ conda_activate conda_deactivate };
use MIP::Package_manager::Pip qw{ check_pip_package pip_install };
use MIP::Program::Compression::Zip qw{ unzip };
use MIP::Program::Download::Wget qw{ wget };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_vcf2cytosure };
}

sub install_vcf2cytosure {

## Function : Install Vcf2cytosure
## Returns  :
## Arguments: $conda_environment       => Conda environment
##          : $conda_prefix_path       => Conda prefix path
##          : $FILEHANDLE              => Filehandle to write to
##          : $noupdate                => Do not update
##          : $program_parameters_href => Hash with Program specific parameters {REF}
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment;
    my $conda_prefix_path;
    my $FILEHANDLE;
    my $noupdate;
    my $quiet;
    my $vcf2cytosure_parameters_href;
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
        noupdate => {
            store       => \$noupdate,
            strict_type => 1,
        },
        program_parameters_href => {
            default     => {},
            required    => 1,
            store       => \$vcf2cytosure_parameters_href,
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
    my $program_version = $vcf2cytosure_parameters_href->{version};
    my $program_url =
      q{https://github.com/NBISweden/vcf2cytosure/archive/} . $program_version . q{.zip};

    ## Set program specific parameters
    my $program_name = q{vcf2cytosure};

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

    say {$FILEHANDLE} q{### Install Vcf2cytosure};

    ## Check if vcf2cytosure has been installed via pip
    my $vcf2cytosure_status = check_pip_package(
        {
            conda_environment => $conda_environment,
            conda_prefix_path => $conda_prefix_path,
            package           => $program_name,
            version           => $program_version,
        }
    );

    ## Check if installation exists and is executable
    if ( -x catfile( $conda_prefix_path, qw{ bin }, $program_name )
        || $vcf2cytosure_status )
    {
        $log->info( $program_name
              . q{ is already installed in the specified conda environment.} );

        if ($noupdate) {

            say {$FILEHANDLE} q{## }
              . $program_name
              . q{ is already installed, skippping reinstallation};
            say {$FILEHANDLE} $NEWLINE;
            return;
        }
        $log->warn(
            q{This will overwrite the current } . $program_name . q{ installation.} );

    }

    ## Only activate conda environment if supplied by user
    if ($conda_environment) {

        ## Activate conda environment
        say {$FILEHANDLE} q{## Activate conda environment};
        conda_activate(
            {
                env_name   => $conda_environment,
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Move to miniconda environment
    say {$FILEHANDLE} q{## Move to conda environment};
    gnu_cd(
        {
            directory_path => $conda_prefix_path,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Download
    say {$FILEHANDLE} q{## Download} . $SPACE . $program_name;
    my $program_zip_path =
      catfile( $conda_prefix_path,
        $program_name . $DASH . $program_version . $DOT . q{zip} );
    wget(
        {
            FILEHANDLE   => $FILEHANDLE,
            outfile_path => $program_zip_path,
            quiet        => $quiet,
            url          => $program_url,
            verbose      => $verbose,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    unzip(
        {
            FILEHANDLE  => $FILEHANDLE,
            force       => 1,
            infile_path => $program_zip_path,
            quiet       => $quiet,
            verbose     => $verbose,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Move to vcf2cytosure directory
    say {$FILEHANDLE} q{## Move to vcf2cytosure directory};
    gnu_cd(
        {
            directory_path => $program_name . $DASH . $program_version,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Pip install the downloaded vcf2cytosure package
    say {$FILEHANDLE} q{## Install};
    pip_install(
        {
            editable      => $DOT,
            FILEHANDLE    => $FILEHANDLE,
            python_module => 1,
            quiet         => $quiet,
            verbose       => $verbose,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Go back to subroutine origin
    say {$FILEHANDLE} q{## Moving back to original working directory};
    gnu_cd(
        {
            directory_path => $pwd,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Deactivate conda environment if conda_environment exists
    if ($conda_environment) {

        say {$FILEHANDLE} q{## Deactivate conda environment};
        conda_deactivate(
            {
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }
    return;
}

1;
