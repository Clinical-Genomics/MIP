package MIP::Check::Installation;

use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## Cpanm
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_existing_installation };
}

## Constants
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

sub check_existing_installation {

## Function : Checks if the program has already been installed and optionally removes the current installation.
##          : Returns "1" if the program is found and a noupdate flag has been provided
## Returns  : $install_check
## Arguments: $conda_environment      => Conda environment
##          : $conda_prefix_path      => Path to conda environment
##          : $FILEHANDLE             => Filehandle to write to
##          : $log                    => Log to write messages to
##          : $noupdate               => Do not update
##          : $program_directory_path => Path to program directory
##          : $program_name                => Program name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment;
    my $conda_prefix_path;
    my $FILEHANDLE;
    my $log;
    my $noupdate;
    my $program_directory_path;
    my $program_name;

    my $tmpl = {
        program_directory_path => {
            required    => 1,
            store       => \$program_directory_path,
            strict_type => 1,
        },
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
            strict_type => 1,
        },
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
        noupdate => {
            allow       => [ undef, 0, 1 ],
            store       => \$noupdate,
            strict_type => 1,
        },
        FILEHANDLE => {
            defined  => 1,
            required => 1,
            store    => \$FILEHANDLE,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Modules
    use MIP::Gnu::Coreutils qw{ gnu_rm };
    use MIP::Gnu::Findutils qw{ gnu_find };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };

    ## Set default for conda environment if undef
    if ( not $conda_environment ) {
        $conda_environment = q{root/base};
    }

    ## Check if installation directory exists
    if ( -d $program_directory_path ) {
        $log->info( $program_name
              . $SPACE
              . q{is already installed in conda environment: }
              . $conda_environment );

        if ($noupdate) {
            $log->info( q{Skipping writting installation instructions for }
                  . $program_name );
            say {$FILEHANDLE}
              q{## Skipping writting installation instructions for }
              . $program_name;
            say {$FILEHANDLE} $NEWLINE;

            return 1;
        }

        $log->warn(
            qq{This will overwrite the current $program_name installation});

        say {$FILEHANDLE} qq{## Removing old $program_name directory};
        gnu_rm(
            {
                FILEHANDLE  => $FILEHANDLE,
                force       => 1,
                infile_path => $program_directory_path,
                recursive   => 1,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        say {$FILEHANDLE} qq{## Removing old $program_name links};
        gnu_find(
            {
                action        => q{-delete},
                FILEHANDLE    => $FILEHANDLE,
                search_path   => catdir( $conda_prefix_path, q{bin} ),
                test_criteria => q{-xtype l},
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    $log->info(
        qq{Writing instructions for $program_name installation via SHELL});

    return 0;
}

1;
