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
use MIP::Program::Tar qw{ tar };
use MIP::Program::Download::Wget qw{ wget };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_root };
}

sub install_root {

## Function : Install ROOT, a CNVnator requirement
## Returns  :
## Arguments: $conda_environment => Conda environment
##          : $conda_prefix_path => Conda prefix path
##          : $filehandle        => Filehandle to write to
##          : $quiet             => Be quiet
##          : $root_binary       => Name of ROOT binary
##          : $verbose           => Set verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment;
    my $conda_prefix_path;
    my $filehandle;
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
        filehandle => {
            defined  => 1,
            required => 1,
            store    => \$filehandle,
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
    say {$filehandle} q{### Install ROOT};

    ## Move to miniconda environment
    gnu_cd(
        {
            directory_path => $conda_prefix_path,
            filehandle     => $filehandle,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Download ROOT
    say {$filehandle} q{## Download ROOT};
    my $root_url = q{https://root.cern.ch/download/} . $root_binary;
    wget(
        {
            filehandle   => $filehandle,
            outfile_path => $root_binary,
            quiet        => $quiet,
            url          => $root_url,
            verbose      => $verbose,
        }
    );

    say {$filehandle} $NEWLINE;

    ## Extract
    say {$filehandle} q{## Extract};
    tar(
        {
            extract    => 1,
            file_path  => $root_binary,
            filehandle => $filehandle,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Use newly installed root during installation
    say {$filehandle} q{source}
      . $SPACE
      . catfile( $conda_prefix_path, qw{ root bin thisroot.sh } )
      . $NEWLINE;

    ## Remove ROOT archive
    say {$filehandle} q{## Removing ROOT archive};
    gnu_rm(
        {
            filehandle  => $filehandle,
            infile_path => $root_binary,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Moving back to original working directory};
    gnu_cd(
        {
            filehandle     => $filehandle,
            directory_path => $pwd,
        }
    );
    say {$filehandle} $NEWLINE;

    return;

}

1;
