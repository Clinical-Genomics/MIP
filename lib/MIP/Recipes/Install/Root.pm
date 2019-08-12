package MIP::Recipes::Install::Root;

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
use MIP::Constants qw{ $NEWLINE $SPACE };
use MIP::Gnu::Bash qw{ gnu_cd };
use MIP::Gnu::Coreutils qw{ gnu_rm };
use MIP::Program::Compression::Tar qw{ tar };
use MIP::Program::Download::Wget qw{ wget };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_root };
}

sub install_root {

## Function : Install ROOT, a CNVnator requirement
## Returns  :
## Arguments: $conda_environment => Conda environment
##          : $conda_prefix_path => Conda prefix path
##          : $FILEHANDLE        => Filehandle to write to
##          : $noupdate          => Do not update
##          : $quiet             => Be quiet
##          : $root_binary       => Name of ROOT binary
##          : $verbose           => Set verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment;
    my $conda_prefix_path;
    my $FILEHANDLE;
    my $noupdate;
    my $quiet;
    my $root_binary;
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
        quiet => {
            allow       => [ undef, 0, 1 ],
            store       => \$quiet,
            strict_type => 1,
        },
        root_binary => {
            required    => 1,
            store       => \$root_binary,
            strict_type => 1,
        },
        verbose => {
            allow       => [ undef, 0, 1 ],
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};
    ## Store original working directory
    my $pwd = cwd();

    ### ROOT installation (prerequisite for CNVnator)
    say {$FILEHANDLE} q{### Install ROOT};

    ## Move to miniconda environment
    gnu_cd(
        {
            directory_path => $conda_prefix_path,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Download ROOT
    say {$FILEHANDLE} q{## Download ROOT};
    my $root_url = q{https://root.cern.ch/download/} . $root_binary;
    wget(
        {
            FILEHANDLE   => $FILEHANDLE,
            outfile_path => $root_binary,
            quiet        => $quiet,
            url          => $root_url,
            verbose      => $verbose,
        }
    );

    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    tar(
        {
            extract    => 1,
            file_path  => $root_binary,
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Use newly installed root during installation
    say {$FILEHANDLE} q{source}
      . $SPACE
      . catfile( $conda_prefix_path, qw{ root bin thisroot.sh } )
      . $NEWLINE;

    ## Remove ROOT archive
    say {$FILEHANDLE} q{## Removing ROOT archive};
    gnu_rm(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $root_binary,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Moving back to original working directory};
    gnu_cd(
        {
            FILEHANDLE     => $FILEHANDLE,
            directory_path => $pwd,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    return;

}

1;
