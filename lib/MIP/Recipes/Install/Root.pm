package MIP::Recipes::Install::Root;

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
    our @EXPORT_OK = qw{ install_root };
}

## Constants
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

sub install_root {

## Function : Install ROOT, a CNVnator requirement
## Returns  : ""
## Arguments: $root_binary       => Name of ROOT binary
##          : $conda_prefix_path => Conda prefix path
##          : $conda_environment => Conda environment
##          : $noupdate          => Do not update
##          : $quiet             => Be quiet
##          : $verbose           => Set verbosity
##          : $FILEHANDLE        => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $root_binary;
    my $conda_prefix_path;
    my $conda_environment;
    my $noupdate;
    my $quiet;
    my $verbose;
    my $FILEHANDLE;

    my $tmpl = {
        root_binary => {
            required    => 1,
            strict_type => 1,
            store       => \$root_binary
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
    use MIP::Gnu::Coreutils qw{ gnu_rm };
    use MIP::Gnu::Bash qw{ gnu_cd };
    use MIP::Program::Download::Wget qw{ wget };
    use MIP::Program::Compression::Tar qw{ tar };

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
    say $FILEHANDLE $NEWLINE;

    ## Download ROOT
    say {$FILEHANDLE} q{## Download ROOT};
    my $root_url = q{https://root.cern.ch/download/} . $root_binary;
    wget(
        {
            url          => $root_url,
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $quiet,
            verbose      => $verbose,
            outfile_path => $root_binary,
        }
    );

    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    tar(
        {
            file_path  => $root_binary,
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
            infile_path => $root_binary,
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

    return;

}

1;
