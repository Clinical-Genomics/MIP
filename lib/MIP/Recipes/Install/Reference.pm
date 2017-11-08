package MIP::Recipes::Install::Reference;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Spec::Functions qw{ catdir catfile };

## Cpanm
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ download_genome_references };
}

## Constants
Readonly my $NEWLINE => qq{\n};

sub download_genome_references {

## Function : Recipe for writing instructions to download genome references
## Returns  : 
## Arguments: $reference_genome_versions_ref => Array with genome versions to download {REF}
##          : $reference_dir_path            => Path to reference directory
##          : $conda_prefix_path             => Conda prefix path
##          : $conda_environment             => Conda environment
##          : $quiet                         => Be quiet
##          : $verbose                       => Set verbosity
##          : $FILEHANDLE                    => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $reference_genome_versions_ref;
    my $reference_dir_path;
    my $conda_prefix_path;
    my $conda_environment;
    my $quiet;
    my $verbose;
    my $FILEHANDLE;

    my $tmpl = {
        reference_genome_versions_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$reference_genome_versions_ref
        },
        reference_dir_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$reference_dir_path
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
    use MIP::Log::MIP_log4perl qw{ retrieve_log };
    use MIP::Package_manager::Conda
      qw{ conda_source_activate conda_source_deactivate };
    use MIP::Program::Download::Download_reference qw{ download_reference };

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => q{mip_install::download_genome_references},
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    say {$FILEHANDLE} q{### Download genome references};
    $log->info(q{Writting download instructions for references});

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

    say {$FILEHANDLE} q{## Generate shell script for reference download};
    download_reference(
        {
            reference_genome_versions_ref => $reference_genome_versions_ref,
            reference_dir_path            => $reference_dir_path,
            FILEHANDLE                    => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Launch bash
    say {$FILEHANDLE} q{## Download references};
    say {$FILEHANDLE} q{bash download_reference.sh} . $NEWLINE;

    ## Cleanup
    say {$FILEHANDLE} q{## Remove generated shell script};
    gnu_rm(
        {
            infile_path => q{download_reference.sh},
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Deactivate conda environment if environment exists
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
