package MIP::Recipes::Install::BootstrapAnn;

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
use MIP::Constants qw{ $LOG_NAME $NEWLINE $SPACE };
use MIP::Gnu::Bash qw{ gnu_cd };
use MIP::Gnu::Coreutils qw{ gnu_chmod gnu_ln };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Package_manager::Conda qw{ conda_activate conda_deactivate };
use MIP::Versionmanager::Git qw{ git_clone };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_bootstrapann };
}

sub install_bootstrapann {

## Function : Install BootstrapAnn
## Returns  :
## Arguments: $conda_environment       => Conda environment
##          : $conda_prefix_path       => Conda prefix path
##          : $FILEHANDLE              => Filehandle to write to
##          : $program_parameters_href => Hash with Program specific parameters {REF}
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment;
    my $conda_prefix_path;
    my $FILEHANDLE;
    my $quiet;
    my $bootstrapann_parameters_href;
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
            store       => \$bootstrapann_parameters_href,
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
    my $program_url     = $bootstrapann_parameters_href->{url};
    my $program_version = $bootstrapann_parameters_href->{version};

    ## Set parameters
    my $program_name           = q{BootstrapAnn};
    my $program_directory_path = catdir( $conda_prefix_path, qw{ share BootstrapAnn } );
    my $executable             = q{BootstrapAnn.py};

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => $LOG_NAME,
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    say {$FILEHANDLE} q{### Install} . $SPACE . $program_name;
    $log->info(qq{Writing instructions for $program_name installation via SHELL});

    ## Check if installation exists and remove directory
    check_existing_installation(
        {
            conda_environment      => $conda_environment,
            conda_prefix_path      => $conda_prefix_path,
            FILEHANDLE             => $FILEHANDLE,
            log                    => $log,
            program_directory_path => $program_directory_path,
            program_name           => $program_name,
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

    ## Download
    say {$FILEHANDLE} q{## Download} . $SPACE . $program_name;
    git_clone(
        {
            FILEHANDLE  => $FILEHANDLE,
            outdir_path => $program_directory_path,
            quiet       => $quiet,
            url         => $program_url,
            verbose     => $verbose,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Change mode
    say {$FILEHANDLE} q{## Make file executable};
    gnu_chmod(
        {
            FILEHANDLE => $FILEHANDLE,
            file_path  => catfile( $program_directory_path, $executable ),
            permission => q{a+x},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Place symlink in bin
    say {$FILEHANDLE} q{## Linking executable};
    gnu_ln(
        {
            FILEHANDLE  => $FILEHANDLE,
            link_path   => catfile( $conda_prefix_path, q{bin}, $executable ),
            target_path => catfile( $program_directory_path, $executable ),
            symbolic    => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

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
