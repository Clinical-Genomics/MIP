package MIP::Recipes::Install::Vep;

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
use List::MoreUtils qw{ any };

## Cpanm
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_vep };
}

## Constants
Readonly my $DOT     => q{.};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

sub install_vep {

## Function : Install varianteffectpredictor
## Returns  : ""
## Arguments: $program_parameters_href => Hash with vep specific parameters {REF}
##          : $conda_prefix_path       => Conda prefix path
##          : $conda_environment       => Conda environment
##          : $noupdate                => Do not update
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity
##          : $FILEHANDLE              => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vep_parameters_href;
    my $conda_prefix_path;
    my $conda_environment;
    my $noupdate;
    my $quiet;
    my $verbose;
    my $FILEHANDLE;

    my $tmpl = {
        program_parameters_href => {
            required    => 1,
            default     => {},
            strict_type => 1,
            store       => \$vep_parameters_href
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
    use MIP::Gnu::Bash qw(gnu_cd);
    use MIP::Gnu::Coreutils qw{ gnu_mkdir gnu_mv gnu_rm };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };
    use MIP::Package_manager::Conda
      qw{ conda_source_activate conda_source_deactivate };
    use MIP::Program::Compression::Tar qw{ tar };
    use MIP::Program::Download::Wget qw{ wget };

    use MIP::Program::Variantcalling::Vep
      qw{ variant_effect_predictor_install };
    use MIP::Versionmanager::Git qw{ git_checkout git_clone };

    ## Unpack parameters
    # Plugins
    my @plugins = @{ $vep_parameters_href->{vep_plugins} };

    # Assembly names to use during --AUTO
    my @assemblies = @{ $vep_parameters_href->{vep_assemblies} };

    # Vep version
    my $vep_version = $vep_parameters_href->{version};

# Run installer without user prompts. Use "a" (API + Faidx/htslib),"l" (Faidx/htslib only), "c" (cache), "f" (FASTA), "p" (plugins) to specify parts to install.
    my $auto = $vep_parameters_href->{vep_auto_flag};

    # Set destination directory for cache files
    my $cache_directory = $vep_parameters_href->{vep_cache_dir};

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => q{mip_install::install_varianteffectpredictor},
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    ## Store original working directory
    my $pwd = cwd();

    my $miniconda_bin_dir = catdir( $conda_prefix_path, q{ensembl-vep} );

    ## Check for existing vep installation
    if ( -d $miniconda_bin_dir ) {

        $log->info(
q{Varianteffectpredictor is already installed in conda environment: }
              . $miniconda_bin_dir );

        if ($noupdate) {

            $log->info(
q{Skipping writting installation instructions for varianteffectpredictor}
            );
            return;
        }
        else {

            $log->warn(
q{This will overwrite the current varianteffectpredictor installation}
            );

            ## Removing varianteffectpredictor
            say {$FILEHANDLE}
              q{## Removing old varianteffectpredictor directory};
            gnu_rm(
                {
                    infile_path => $miniconda_bin_dir,
                    force       => 1,
                    recursive   => 1,
                    FILEHANDLE  => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }
    }
    else {

        $log->info(q{Writting install instructions for varianteffectpredictor});
    }

    ## Install VEP
    say {$FILEHANDLE} q{### Install varianteffectpredictor};

    ## Only activate conda environment if supplied by user
    if ($conda_environment) {

        ## Activate conda environment
        say {$FILEHANDLE} q{## Activate conda environment};
        conda_source_activate(
            {
                FILEHANDLE => $FILEHANDLE,
                env_name   => $conda_environment,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Make sure that the cache directory exists
    gnu_mkdir(
        {
            indirectory_path => $cache_directory,
            parents          => 1,
            FILEHANDLE       => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Move to miniconda environment
    gnu_cd(
        {
            directory_path => catdir($conda_prefix_path),
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Git clone
    say {$FILEHANDLE} q{## Git clone VEP};
    git_clone(
        {
            url        => q{https://github.com/Ensembl/ensembl-vep.git},
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Move to vep directory
    gnu_cd(
        {
            directory_path => catdir(q{ensembl-vep}),
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Checkout release branch
    say {$FILEHANDLE} q{## Checkout release branch};
    git_checkout(
        {
            branch     => catdir( q{release}, $vep_version ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Install VEP
    say {$FILEHANDLE} q{## Install VEP};

    variant_effect_predictor_install(
        {
            plugins_ref     => \@plugins,
            species_ref     => [qw{ homo_sapiens }],
            auto            => $auto,
            cache_directory => $cache_directory,
            assembly        => $assemblies[0],
            FILEHANDLE      => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## If more than one assembly requested
    if (   ( scalar @assemblies > 1 )
        && ( $auto =~ / [cf] /xsm ) )
    {

        # Find last index of array and initate
        Readonly my $NUMBER_OF_ASSEMBLIES => $#assemblies;

      ASSEMBLY:
        for my $assembly_version ( 1 .. $NUMBER_OF_ASSEMBLIES ) {
            ## Skip first assembly since it is already installed above

            say {$FILEHANDLE}
              q{## Install additional VEP cache assembly version};

            variant_effect_predictor_install(
                {
                    species_ref     => [qw{ homo_sapiens }],
                    auto            => q{cf},
                    cache_directory => $cache_directory,
                    assembly        => $assemblies[$assembly_version],
                    FILEHANDLE      => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }
    }

    # Initate
    my $vep_plugin_dir = catdir( q{$HOME}, $DOT, qw{vep Plugins} );

    if (@plugins) {

        if ( any { $_ eq q{MaxEntScan} } @plugins ) {

            ## Add MaxEntScan required text file
            say {$FILEHANDLE} q{## Add MaxEntScan required text file};
            wget(
                {
                    url =>
q{http://genes.mit.edu/burgelab/maxent/download/fordownload.tar.gz},
                    FILEHANDLE => $FILEHANDLE,
                    quiet      => $quiet,
                    verbose    => $verbose,
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            # Unpack
            tar(
                {
                    extract     => 1,
                    filter_gzip => 1,
                    file_path   => catfile(q{fordownload.tar.gz}),
                    FILEHANDLE  => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            gnu_mv(
                {
                    infile_path  => q{fordownload},
                    outfile_path => catfile( $cache_directory, ),
                    force        => 1,
                    FILEHANDLE   => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} $NEWLINE;

        }
        if ( any { $_ eq q{LoFtool} } @plugins ) {

            ## Add LofTool required text file
            say {$FILEHANDLE} q{## Add LofTool required text file};
            wget(
                {
                    url =>
q{https://raw.githubusercontent.com/Ensembl/VEP_plugins/master/LoFtool_scores.txt},
                    FILEHANDLE => $FILEHANDLE,
                    quiet      => $quiet,
                    verbose    => $verbose,
                    outfile_path =>
                      catfile( $vep_plugin_dir, q{LoFtool_scores.txt} ),
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }
    }

    ## Clean up
    say {$FILEHANDLE} q{## Clean up};
    gnu_rm(
        {
            infile_path => catdir(
                $conda_prefix_path,
                q{VariantEffectPredictor-} . $vep_version . $DOT . q{zip}
            ),
            force      => 1,
            FILEHANDLE => $FILEHANDLE,
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
        conda_source_deactivate(
            {
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE x 2;
    }
    return;
}

1;
