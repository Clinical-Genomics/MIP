package MIP::Recipes::Install::Picard;

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
use MIP::Gnu::Coreutils qw{ gnu_ln gnu_mv gnu_rm };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Program::Compression::Zip qw{ unzip };
use MIP::Program::Download::Wget qw{ wget };
use MIP::Script::Utils qw{ create_temp_dir };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_picard };
}

sub install_picard {

## Function : Install Picard
## Returns  : ""
## Arguments: $conda_environment       => Conda environment
##          : $conda_prefix_path       => Conda prefix path
##          : $FILEHANDLE              => Filehandle to write to
##          : $program_parameters_href => Hash with Picard specific parameters {REF}
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment;
    my $conda_prefix_path;
    my $FILEHANDLE;
    my $picard_parameters_href;
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
            required => 1,
            store    => \$FILEHANDLE,
            defined  => 1,
        },
        program_parameters_href => {
            default     => {},
            required    => 1,
            store       => \$picard_parameters_href,
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
    my $picard_version = $picard_parameters_href->{version};

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

    say {$FILEHANDLE} q{### Install Picard};
    $log->info(qq{Writing instructions for Picard installation via SHELL});

    ## Check if installation exists and remove directory
    my $picard_dir_path =
      catdir( $conda_prefix_path, q{share}, q{picard-tools-} . $picard_version );
    check_existing_installation(
        {
            conda_environment      => $conda_environment,
            conda_prefix_path      => $conda_prefix_path,
            FILEHANDLE             => $FILEHANDLE,
            log                    => $log,
            program_directory_path => $picard_dir_path,
            program_name           => q{Picard},
        }
    );

    ## Creating temporary install directory
    say {$FILEHANDLE} q{## Create temporary picard install directory};
    my $temp_dir = create_temp_dir( { FILEHANDLE => $FILEHANDLE } );
    say {$FILEHANDLE} $NEWLINE;

    ## Download
    say {$FILEHANDLE} q{## Download Picard};
    my $url =
        q{https://github.com/broadinstitute/picard/releases/download/}
      . $picard_version
      . q{/picard-tools-}
      . $picard_version
      . $DOT . q{zip};
    my $picard_zip_path =
      catfile( $temp_dir, q{picard-tools-} . $picard_version . $DOT . q{zip} );
    wget(
        {
            FILEHANDLE   => $FILEHANDLE,
            outfile_path => $picard_zip_path,
            quiet        => $quiet,
            url          => $url,
            verbose      => $verbose,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    unzip(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $picard_zip_path,
            outdir_path => $temp_dir,
            quiet       => $quiet,
            verbose     => $verbose,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Make available from conda environment
    say {$FILEHANDLE} q{## Make available from conda environment};
    gnu_mv(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => catdir( $temp_dir, q{picard-tools-} . $picard_version ),
            outfile_path => catdir( $conda_prefix_path, q{share} ),
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Specifying target and link paths
    say {$FILEHANDLE} q{## Linking picard.jar};
    my $target_path =
      catfile( $conda_prefix_path, q{share}, q{picard-tools-} . $picard_version,
        q{picard.jar} );
    my $link_path = catfile( $conda_prefix_path, qw{ bin picard.jar } );
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

    ## Remove the temporary install directory
    gnu_rm(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $temp_dir,
            recursive   => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    return;
}

1;
