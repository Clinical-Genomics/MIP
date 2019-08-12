package MIP::Recipes::Install::Reference;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
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
use MIP::Constants qw{ $DOT $LOG $NEWLINE $UNDERSCORE };
use MIP::Gnu::Coreutils qw{ gnu_rm };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Package_manager::Conda qw{ conda_activate conda_deactivate };
use MIP::Program::Download::Download_reference qw{ download_reference };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ download_genome_references };
}

sub download_genome_references {

## Function : Recipe for writing instructions to download genome references
## Returns  :
## Arguments: $conda_environment             => Conda environment
##          : $conda_prefix_path             => Conda prefix path
##          : $FILEHANDLE                    => Filehandle to write to
##          : $pipeline                      => Pipeline for which to download references
##          : $quiet                         => Be quiet
##          : $reference_dir_path            => Path to reference directory
##          : $reference_genome_versions_ref => Array with genome versions to download {REF}
##          : $verbose                       => Set verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment;
    my $conda_prefix_path;
    my $FILEHANDLE;
    my $pipeline;
    my $quiet;
    my $reference_dir_path;
    my $reference_genome_versions_ref;
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
        pipeline => {
            allow       => [qw{ rare_disease rna }],
            defined     => 1,
            required    => 1,
            store       => \$pipeline,
            strict_type => 1,
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            store       => \$quiet,
            strict_type => 1,
        },
        reference_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$reference_dir_path,
            strict_type => 1,
        },
        reference_genome_versions_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$reference_genome_versions_ref,
            strict_type => 1,
        },
        verbose => {
            allow       => [ undef, 0, 1 ],
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => $LOG,
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    say {$FILEHANDLE} q{### Download genome references};
    $log->info(q{Writting download instructions for references});

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

    ## Set path to config file
    my $config_file_path = catfile( $Bin, q{definitions},
            q{download}
          . $UNDERSCORE
          . $pipeline
          . $UNDERSCORE
          . q{parameters}
          . $DOT
          . q{yaml} );

    say {$FILEHANDLE} q{## Generate shell script for reference download};
    download_reference(
        {
            FILEHANDLE                    => $FILEHANDLE,
            config_file_path              => $config_file_path,
            pipeline                      => $pipeline,
            reference_dir_path            => $reference_dir_path,
            reference_genome_versions_ref => $reference_genome_versions_ref,
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
            FILEHANDLE  => $FILEHANDLE,
            infile_path => q{download_reference.sh},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Deactivate conda environment if environment exists
    if ($conda_environment) {

        say {$FILEHANDLE} q{## Deactivate conda environment};
        conda_deactivate(
            {
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    print {$FILEHANDLE} $NEWLINE;

    return;
}

1;
