package MIP::Recipes::Install::Vep;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use List::MoreUtils qw{ any };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## CPAN
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Check::Installation qw{ check_existing_installation };
use MIP::Constants qw{ $DASH $DOT $LOG $NEWLINE $SPACE };
use MIP::Gnu::Bash qw{ gnu_cd gnu_unset };
use MIP::Gnu::Coreutils qw{ gnu_ln gnu_mkdir gnu_rm };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Package_manager::Conda qw{ conda_activate conda_deactivate };
use MIP::Program::Compression::Tar qw{ tar };
use MIP::Program::Download::Wget qw{ wget };
use MIP::Program::Variantcalling::Vep qw{ variant_effect_predictor_install };
use MIP::Versionmanager::Git qw{ git_checkout git_clone };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.14;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_vep };
}

sub install_vep {

## Function : Install varianteffectpredictor
## Returns  :
## Arguments: $conda_environment       => Conda environment
##          : $conda_prefix_path       => Conda prefix path
##          : $FILEHANDLE              => Filehandle to write to
##          : $noupdate                => Do not update
##          : $program_parameters_href => Hash with vep specific parameters {REF}
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment;
    my $conda_prefix_path;
    my $FILEHANDLE;
    my $noupdate;
    my $quiet;
    my $vep_parameters_href;
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
        program_parameters_href => {
            default     => {},
            required    => 1,
            store       => \$vep_parameters_href,
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
    # Assembly names to use during --AUTO
    my @assemblies = @{ $vep_parameters_href->{vep_assemblies} };

    # Plugins
    my @plugins = @{ $vep_parameters_href->{vep_plugins} };

    # Species names to use during --AUTO
    my @species = @{ $vep_parameters_href->{vep_species} };

    # Vep version
    my $vep_version = $vep_parameters_href->{version};

# Run installer without user prompts. Use "a" (API + Faidx/htslib),"l" (Faidx/htslib only), "c" (cache), "f" (FASTA), "p" (plugins) to specify parts to install.
    my $auto = $vep_parameters_href->{vep_auto_flag};

    # Set destination directory for cache files
    my $cache_directory;
    if ( not $vep_parameters_href->{vep_cache_dir} ) {
        $cache_directory =
          catdir( $conda_prefix_path, q{ensembl-tools-release-} . $vep_version,
            q{cache} );
    }
    else {
        $cache_directory = catdir( $vep_parameters_href->{vep_cache_dir} );
    }

    # Set vep api installation directory
    my $vep_dir_path = catdir( $conda_prefix_path, q{ensembl-vep} );

    # Set vep plugin directory
    my $vep_plugin_dir = catdir( $cache_directory, q{Plugins} );

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

    ## Install VEP
    say {$FILEHANDLE} q{### Install varianteffectpredictor};

    ## Vipe the api if a reinstallation has been requested
    my $install_check;
    if ( $auto =~ m/[al]/xms ) {
        ## Check if installation exists and remove directory unless a noupdate flag is provided
        $install_check = check_existing_installation(
            {
                conda_environment      => $conda_environment,
                conda_prefix_path      => $conda_prefix_path,
                FILEHANDLE             => $FILEHANDLE,
                log                    => $log,
                noupdate               => $noupdate,
                program_directory_path => $vep_dir_path,
                program_name           => q{VEP-api},
            }
        );
        print {$FILEHANDLE} $NEWLINE;
    }

    # Return if all the directories is found and a noupdate flag has been provided
    if ($install_check) {
        ## Skip api installation but continue with rest
        $auto =~ s/[al]//gxms;
        if ( not $auto =~ m/[cfp]/xms ) {
            say {$FILEHANDLE} $NEWLINE;
            return;
        }
    }
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

    ## Make sure that the VEP:s INSTALL.pl exist if the user has selected to skip api installation
    if ( ( $auto =~ m/[cfp]/xms ) && ( not $auto =~ m/[al]/xms ) ) {

        if ( not -s catfile( $vep_dir_path, q{INSTALL.pl} ) ) {
            $log->fatal(
q{Can't install cache, fasta files or plugins without a VEP installation file!}
            );
            $log->fatal(
q{Please add the [a] and/or [l] flag to --vep_auto_flag when running mip_install.pl}
            );
            exit 1;
        }
    }

    ## Make sure that the cache directory exists
    if ( not -d $cache_directory ) {
        say {$FILEHANDLE} q{## Create cache directory};
        gnu_mkdir(
            {
                FILEHANDLE       => $FILEHANDLE,
                indirectory_path => $cache_directory,
                parents          => 1,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Set LD_LIBRARY_PATH for VEP installation
    say {$FILEHANDLE} q{LD_LIBRARY_PATH=}
      . $conda_prefix_path
      . q{/lib/:$LD_LIBRARY_PATH};
    say {$FILEHANDLE} q{export LD_LIBRARY_PATH} . $NEWLINE;

    ## Download VEP
    if ( $auto =~ m/[al]/xms ) {
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
                FILEHANDLE => $FILEHANDLE,
                url        => q{https://github.com/Ensembl/ensembl-vep.git},
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Move to vep directory
    gnu_cd(
        {
            directory_path => $vep_dir_path,
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

    ## Don't install plugins unless specified in the auto flag
    if ( not $auto =~ m/p/xms ) {
        @plugins = qw{};
    }
    variant_effect_predictor_install(
        {
            assembly        => $assemblies[0],
            auto            => $auto,
            cache_directory => $cache_directory,
            FILEHANDLE      => $FILEHANDLE,
            no_htslib       => 0,
            plugins_ref     => \@plugins,
            species_ref     => \@species,
            version         => $vep_version,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## If more than one assembly requested
    if (   ( scalar @assemblies > 1 )
        && ( $auto =~ / [cf] /xsm ) )
    {

        ## Remove api and plugins from the auto flag
        my $cf_auto = $auto;
        $cf_auto =~ s/[alp]//gxms;

        # Find last index of array and initate
        Readonly my $NUMBER_OF_ASSEMBLIES => $#assemblies;

      ASSEMBLY:
        for my $assembly_version ( 1 .. $NUMBER_OF_ASSEMBLIES ) {
            ## Skip first assembly since it is already installed above

            say {$FILEHANDLE} q{## Install additional VEP cache assembly version};

            variant_effect_predictor_install(
                {
                    assembly        => $assemblies[$assembly_version],
                    auto            => $cf_auto,
                    cache_directory => $cache_directory,
                    cache_version   => $vep_version,
                    FILEHANDLE      => $FILEHANDLE,
                    species_ref     => \@species,
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }
    }

    ## Install and download extra plugin files
    if ( @plugins && $auto =~ m/p/xms ) {

        if ( any { $_ eq q{MaxEntScan} } @plugins ) {

            ## Add MaxEntScan required text file
            say {$FILEHANDLE} q{## Add MaxEntScan required text file};
            my $maxent_file_path = catfile( $vep_plugin_dir, q{fordownload.tar.gz} );
            wget(
                {
                    FILEHANDLE   => $FILEHANDLE,
                    outfile_path => $maxent_file_path,
                    quiet        => $quiet,
                    url =>
                      q{http://genes.mit.edu/burgelab/maxent/download/fordownload.tar.gz},
                    verbose => $verbose,
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            # Unpack
            tar(
                {
                    extract           => 1,
                    outdirectory_path => dirname($maxent_file_path),
                    FILEHANDLE        => $FILEHANDLE,
                    filter_gzip       => 1,
                    file_path => catfile( $vep_plugin_dir, q{fordownload.tar.gz} ),
                }
            );
            say {$FILEHANDLE} $NEWLINE;

            ## Clean up
            say {$FILEHANDLE} q{## Clean up};
            gnu_rm(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    force       => 1,
                    infile_path => catfile( $vep_plugin_dir, q{fordownload.tar.gz} ),
                }
            );
            say {$FILEHANDLE} $NEWLINE;

        }
        if ( any { $_ eq q{LoFtool} } @plugins ) {

            ## Add LofTool required text file
            say {$FILEHANDLE} q{## Add LofTool required text file};
            wget(
                {
                    FILEHANDLE   => $FILEHANDLE,
                    outfile_path => catfile( $vep_plugin_dir, q{LoFtool_scores.txt} ),
                    quiet        => $quiet,
                    url =>
q{https://raw.githubusercontent.com/Ensembl/VEP_plugins/master/LoFtool_scores.txt},
                    verbose => $verbose,
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }
        if ( any { $_ eq q{ExACpLI} } @plugins ) {

            ## Add ExACpLI required text file
            say {$FILEHANDLE} q{## Add pLI required value file};
            wget(
                {
                    FILEHANDLE   => $FILEHANDLE,
                    outfile_path => catfile( $vep_plugin_dir, q{ExACpLI_values.txt} ),
                    quiet        => $quiet,
                    url =>
q{https://raw.githubusercontent.com/Ensembl/VEP_plugins/master/ExACpLI_values.txt},
                    verbose => $verbose,
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }
    }

    if ( $auto =~ m/[al]/xms ) {

        ## Make VEP-api available from conda environment
        say {$FILEHANDLE} q{## Make VEP-api available from conda environment};
        my $target_path = catfile( $vep_dir_path,      q{vep} );
        my $link_path   = catfile( $conda_prefix_path, qw{ bin vep } );
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
    }

    ## Unset LD_LIBRARY_PATH as to not pollute the rest of the installation
    gnu_unset(
        {
            bash_variable => q{LD_LIBRARY_PATH},
            FILEHANDLE    => $FILEHANDLE,
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
        conda_deactivate(
            {
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE x 2;
    }
    return;
}

1;
