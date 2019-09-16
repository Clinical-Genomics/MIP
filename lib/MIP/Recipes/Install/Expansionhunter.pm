package MIP::Recipes::Install::Expansionhunter;

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

## MIPs lib/
use MIP::Check::Installation qw{ check_existing_installation };
use MIP::Constants qw{ $LOG_NAME $NEWLINE };
use MIP::Gnu::Coreutils qw{ gnu_ln gnu_rm };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Program::Compression::Tar qw{ tar };
use MIP::Program::Download::Wget qw{ wget };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_expansionhunter };
}

sub install_expansionhunter {

## Function : Install Expansion Hunter
## Returns  :
## Arguments: $program_parameters_href => Hash with Expansion Hunter specific parameters {REF}
##          : $conda_prefix_path       => Conda prefix path
##          : $conda_environment       => Conda environment
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity
##          : $FILEHANDLE              => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment;
    my $conda_prefix_path;
    my $FILEHANDLE;
    my $quiet;
    my $expansionhunter_parameters_href;
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
            store       => \$expansionhunter_parameters_href,
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
    my $expansionhunter_version = $expansionhunter_parameters_href->{version};

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => $LOG_NAME,
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    say {$FILEHANDLE} q{### Install Expansion Hunter};
    $log->info(qq{Writing instructions for Expansion Hunter installation via SHELL});

    ## Check if installation exists and remove directory
    my $expansionhunter_dir =
      catdir( $conda_prefix_path, q{share},
        q{ExpansionHunter-v} . $expansionhunter_version . q{-linux_x86_64} );
    check_existing_installation(
        {
            conda_environment      => $conda_environment,
            conda_prefix_path      => $conda_prefix_path,
            FILEHANDLE             => $FILEHANDLE,
            log                    => $log,
            program_directory_path => $expansionhunter_dir,
            program_name           => q{Expansion Hunter},
        }
    );

    ## Download
    say {$FILEHANDLE} q{## Download Expansion Hunter};
    my $url =
        q{https://github.com/Illumina/ExpansionHunter/releases/download/v}
      . $expansionhunter_version
      . q{/ExpansionHunter-v}
      . $expansionhunter_version
      . q{-linux_x86_64.tar.gz};
    my $expansionhunter_download_path = catfile( $conda_prefix_path, q{share},
        q{ExpansionHunter-v} . $expansionhunter_version . q{-linux_x86_64.tar.gz} );
    wget(
        {
            FILEHANDLE   => $FILEHANDLE,
            outfile_path => $expansionhunter_download_path,
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
            file_path         => $expansionhunter_download_path,
            FILEHANDLE        => $FILEHANDLE,
            outdirectory_path => catdir( $conda_prefix_path, q{share} ),
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Create link in conda environment bin to binary in Expansion Hunter folder
    say {$FILEHANDLE} q{## Create softlinks to binary};
    gnu_ln(
        {
            link_path   => catfile( $conda_prefix_path,   q{bin} ),
            target_path => catfile( $expansionhunter_dir, qw{ bin ExpansionHunter } ),
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
            infile_path => $expansionhunter_download_path,
            recursive   => 1,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE x 2;

    return;
}

1;
