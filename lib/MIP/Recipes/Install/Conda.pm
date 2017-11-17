package MIP::Recipes::Install::Conda;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Basename qw{ dirname };

## CPANM
use Readonly;

## Constants
Readonly my $SPACE      => q{ };
Readonly my $NEWLINE    => qq{\n};
Readonly my $DOT        => q{.};
Readonly my $UNDERSCORE => q{_};
Readonly my $BACKTICK   => q{`};

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.0.6;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ setup_conda_env install_bioconda_packages };

}

sub setup_conda_env {

## Function  : Creates necessary conda environment and install package(s) from the default channel.
## Returns   :
## Arguments : $conda_packages_href => Hash with conda packages and their version numbers {REF}
##           : $conda_env           => Name of conda environment
##           : $conda_env_path      => Path to conda environment (default: conda root)
##           : $FILEHANDLE          => Filehandle to write to
##           : $conda_update        => Update Conda if defined
##           : $quiet               => Log only warnings and above
##           : $verbose             => Log debug messages

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_packages_href;
    my $conda_env;
    my $conda_env_path;
    my $FILEHANDLE;
    my $conda_update;
    my $quiet;
    my $verbose;

    my $tmpl = {
        conda_packages_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$conda_packages_href,
        },
        conda_env => {
            strict_type => 1,
            store       => \$conda_env,
        },
        conda_env_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$conda_env_path,
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE,
        },
        conda_update => {
            store => \$conda_update,
        },
        quiet => {
            allow => [ undef, 0, 1 ],
            store => \$quiet,
        },
        verbose => {
            allow => [ undef, 0, 1 ],
            store => \$verbose,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Unix qw{ check_binary_in_path };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };
    use MIP::Package_manager::Conda
      qw{ conda_check_env_status conda_create conda_install conda_update };

    ## Get logger
    my $log = retrieve_log(
        {
            log_name => q{mip_install::setup_conda_env},
            verbose  => $verbose,
            quiet    => $quiet,
        }
    );

    ## Scan the PATH for conda
    check_binary_in_path(
        {
            binary => q{conda},
            log    => $log
        }
    );

    ## Check for active conda environment (exit if true)
    my $env_status = conda_check_env_status();
    if ($env_status) {
        $log->fatal( q{Found activate conda env: } . $env_status );
        $log->fatal(
            q{Run 'source deactivate' prior to running installation script});
        exit 1;
    }

    ## Optionally update conda
    if ($conda_update) {
        say {$FILEHANDLE} q{## Updating Conda};
        conda_update(
            {
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Setup the conda environment
    # Create arrray of packages that are to be installed
    my @packages = _create_package_array(
        {
            package_href              => $conda_packages_href,
            package_version_separator => q{=},
        }
    );

    if ( defined $conda_env ) {
        ## Check for existing conda environment
        if ( not -d $conda_env_path ) {
            ## Create conda environment and install packages
            $log->info(
                q{Writing installtion instructions for environment: }
                  . $conda_env );
            say {$FILEHANDLE} q{## Creating conda environment: }
              . $conda_env
              . q{ and install packages};
            conda_create(
                {
                    env_name     => $conda_env,
                    packages_ref => \@packages,
                    FILEHANDLE   => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }
        else {
            $log->warn( q{Conda environment: }
                  . $conda_env
                  . $SPACE
                  . q{already exists} );
            $log->warn(
                q{Will try to install packages into existing environment});
            $log->info(
q{Writing installtion instructions for conda packages to environment: }
                  . $conda_env );
            say {$FILEHANDLE}
              q{## Installing conda packages into existing environment};
            conda_install(
                {
                    packages_ref => \@packages,
                    FILEHANDLE   => $FILEHANDLE,
                    env_name     => $conda_env,
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }
    }
    else {
        $log->info(
q{Writing instructions for installing and/or updating packages in conda root}
        );
        say {$FILEHANDLE}
          q{## Installing and/or updating python and packages in conda root};
        conda_install(
            {
                packages_ref => \@packages,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }
    return;
}

sub install_bioconda_packages {

## Function  : Install conda packages from the bioconda channel into a conda environment.
## Returns   :
## Arguments : $bioconda_packages_href     => Hash holding bioconda packages and their version numbers {REF}
##           : $snpeff_genome_versions_ref => Array with the genome versions for the snpeff databases {REF}
##           : $conda_env                  => Name of conda environment
##           : $conda_env_path             => Path to conda environment (default: conda root)
##           : $FILEHANDLE                 => Filehandle to write to
##           : $quiet                      => Log only warnings and above
##           : $verbose                    => Log debug messages

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bioconda_packages_href;
    my $snpeff_genome_versions_ref;
    my $conda_env;
    my $conda_env_path;
    my $FILEHANDLE;
    my $quiet;
    my $verbose;

    my $tmpl = {
        bioconda_packages_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$bioconda_packages_href,
        },
        snpeff_genome_versions_ref => {
            required    => 1,
            default     => [],
            defined     => 1,
            strict_type => 1,
            store       => \$snpeff_genome_versions_ref
        },
        conda_env => {
            strict_type => 1,
            store       => \$conda_env,
        },
        conda_env_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$conda_env_path,
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE,
        },
        quiet => {
            allow => [ undef, 0, 1 ],
            store => \$quiet,
        },
        verbose => {
            allow => [ undef, 0, 1 ],
            store => \$verbose,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments};

    use MIP::Gnu::Coreutils qw{ gnu_ln };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };
    use MIP::Package_manager::Conda qw{ conda_install };

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => q{mip_install::install_bioconda_packages},
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    ## Create an array for bioconda packages that are to be installed from provided hash
    my @packages = _create_package_array(
        {
            package_href              => $bioconda_packages_href,
            package_version_separator => q{=},
        }
    );

    ## Install bioconda packages
    $log->info(q{Writing installation instructions for Bioconda packages});
    say {$FILEHANDLE} q{## Installing bioconda modules in conda environment};
    conda_install(
        {
            FILEHANDLE    => $FILEHANDLE,
            packages_ref  => \@packages,
            conda_channel => q{bioconda},
            env_name      => $conda_env,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Link bioconda packages
    # Creating target-link paths
    my %target_link_paths = _create_target_link_paths(
        {
            bioconda_packages_href => $bioconda_packages_href,
            conda_env_path         => $conda_env_path,
            FILEHANDLE             => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} q{## Creating symbolic links for bioconda packages};
  TARGET_AND_LINK_PATHS:
    while ( my ( $target_path, $link_path ) = each %target_link_paths ) {
        gnu_ln(
            {
                FILEHANDLE  => $FILEHANDLE,
                target_path => $target_path,
                link_path   => $link_path,
                symbolic    => 1,
                force       => 1,
            }
        );
        print {$FILEHANDLE} $NEWLINE;
    }
    print {$FILEHANDLE} $NEWLINE;

    ## Custom solutions for BWA, SnpEff and Manta
    ## Copying files, downloading necessary databases and make files executable
    finish_bioconda_package_install(
        {
            bioconda_packages_href     => $bioconda_packages_href,
            conda_env                  => $conda_env,
            conda_env_path             => $conda_env_path,
            FILEHANDLE                 => $FILEHANDLE,
            snpeff_genome_versions_ref => $snpeff_genome_versions_ref,
            verbose                    => $verbose,
            quiet                      => $quiet,
        }
    );

    ## Unset variables
    say {$FILEHANDLE} q{## Unset variables};
    my %program_path_aliases = (
        bwakit  => q{BWAKIT_PATH},
        snpeff  => q{SNPEFF_PATH},
        snpsift => q{SNPSIFT_PATH},
        manta   => q{MANTA_PATH},
        picard  => q{PICARD_PATH},
    );

  PROGRAM:
    foreach my $program ( keys %program_path_aliases ) {
        # Check if the program has been set to be installed via shell and
        # thus has been removed from the bioconda_packages hash
        next PROGRAM if ( not $bioconda_packages_href->{$program} );
        say {$FILEHANDLE} q{unset} . $SPACE . $program_path_aliases{$program};
    }
    say {$FILEHANDLE} $NEWLINE;

    return;
}

sub finish_bioconda_package_install {

## Function  : Custom solutions to finish the install of BWA, SnpEff and Manta
## Returns   :
## Arguments : $bioconda_packages_href     => Hash with bioconda packages {REF}
##           : $snpeff_genome_versions_ref => Array with the genome versions for the snpeff databases {REF}
##           : $conda_env_path             => Path to conda environment
##           : $FILEHANDLE                 => Filehandle to write to
##           : $conda_env                  => Name of conda env

    my ($arg_href) = @_;

    ## Flatten arguments
    my $bioconda_packages_href;
    my $snpeff_genome_versions_ref;
    my $conda_env_path;
    my $FILEHANDLE;
    my $conda_env;
    my $quiet;
    my $verbose;

    my $tmpl = {
        bioconda_packages_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$bioconda_packages_href
        },
        snpeff_genome_versions_ref => {
            required    => 1,
            default     => [],
            defined     => 1,
            strict_type => 1,
            store       => \$snpeff_genome_versions_ref
        },
        conda_env_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$conda_env_path,
        },
        FILEHANDLE => {
            required => 1,
            defined  => 1,
            store    => \$FILEHANDLE
        },
        conda_env => {
            required => 1,
            store    => \$conda_env,
        },
        verbose => {
            allow => [ undef, 0, 1 ],
            store => \$verbose,
        },
        quiet => {
            allow => [ undef, 0, 1 ],
            store => \$quiet,
        }
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use File::Spec::Functions qw{ catdir catfile };
    use IPC::Cmd qw{ run };
    use MIP::Gnu::Coreutils qw{ gnu_cp gnu_chmod gnu_rm };
    use MIP::Package_manager::Conda
      qw{ conda_source_activate conda_source_deactivate };
    use MIP::Program::Variantcalling::SnpEff qw{ snpeff_download };
    use MIP::Recipes::Install::Gatk qw{ gatk_download };
    use MIP::Recipes::Install::SnpEff qw{ check_mt_codon_table };

    ## Only activate conda environment if supplied by user
    if ($conda_env) {
        ## Activate conda environment
        say {$FILEHANDLE} q{## Activate conda environment};
        conda_source_activate(
            {
                FILEHANDLE => $FILEHANDLE,
                env_name   => $conda_env,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Custom BWA
    say {$FILEHANDLE} q{## Custom BWA solutions};
    my $infile_path  = catdir( q/${BWAKIT_PATH}/, q{resource-human-HLA} );
    my $outfile_path = catdir( $conda_env_path,     q{bin} );
    gnu_cp(
        {
            FILEHANDLE   => $FILEHANDLE,
            recursive    => 1,
            force        => 1,
            infile_path  => $infile_path,
            outfile_path => $outfile_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Check if snpeff has been set to be installed via shell
    if ( $bioconda_packages_href->{snpeff} ) {
        ## Custom snpeff - Download necessary databases
        ## Check and if required add the vertebrate mitochondrial codon table to snpeff config
        say {$FILEHANDLE} q{## Custom SnpEff solutions};

        ## Get the full version, including patch for snpeff
        my $command =
            qq{conda install --dry-run snpeff=$bioconda_packages_href->{snpeff}}
          . $SPACE
          . q{| grep 'snpeff' | awk '{ print $2 }'};
        my $version;
        run(
            command => $command,
            buffer  => \$version
        );
        chomp $version;

      SNPEFF_GENOME_VERSION:
        foreach my $genome_version ( @{$snpeff_genome_versions_ref} ) {

            my $share_dir =
              catdir( $conda_env_path, q{share}, q{snpeff-} . $version );
            check_mt_codon_table(
                {
                    FILEHANDLE     => $FILEHANDLE,
                    share_dir      => $share_dir,
                    config_file    => q{snpEff.config},
                    genome_version => $genome_version,
                    verbose        => $verbose,
                    quiet          => $quiet,
                }
            );

            my $snpeff_genome_dir =
              catdir( $share_dir, q{data}, $genome_version );
            next SNPEFF_GENOME_VERSION if ( -d $snpeff_genome_dir );

            ## Write instructions to download snpeff database.
            ## This is done by install script to avoid race conditin when doing first analysis run in MIP
            say {$FILEHANDLE} q{## Downloading snpeff database};
            my $jar_path = catfile( $conda_env_path, qw{ bin snpEff.jar} );
            my $config_file_path =
              catfile( $conda_env_path, qw{bin snpEff.config} );
            snpeff_download(
                {
                    FILEHANDLE              => $FILEHANDLE,
                    genome_version_database => $genome_version,
                    jar_path                => $jar_path,
                    config_file_path        => $config_file_path,
                    temp_directory          => 1,
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }
    }

    ## Custom manta
    # Make file executable
    say {$FILEHANDLE} q{## Changing mode of configManta.py to executable};
    my $file_path = catfile( $conda_env_path, qw{bin configManta.py} );
    gnu_chmod(
        {
            file_path  => $file_path,
            permission => q{a+x},
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Custom GATK
    say {$FILEHANDLE} q{## Custom GATK solutions};

    ## Download gatk .tar.bz2
    my $gatk_tar_path = gatk_download(
        {
            gatk_version => $bioconda_packages_href->{gatk},
            verbose      => $verbose,
            quiet        => $quiet,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    ## Hard coding here since GATK 4.0 will be open source.
    ## Then this step will be unnecessary
    say {$FILEHANDLE} q{gatk-register} . $SPACE . $gatk_tar_path . $NEWLINE;

    gnu_rm(
        {
            infile_path => dirname($gatk_tar_path),
            force       => 1,
            recursive   => 1,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE x 2;

    ## Deactivate conda environment if conda_environment exists
    if ($conda_env) {
        say {$FILEHANDLE} q{## Deactivate conda environment};
        conda_source_deactivate(
            {
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }
    return;
}

sub _create_package_array {

##Function  : Takes a reference to hash of packages and creates an array with
##          : package and version joined with a supplied separator if value is defined.
##          : Also checks that the version number makes sense
##Returns   : "@packages"
##Arguments : $package_href              => Hash with packages {Hash}
##          : $package_version_separator => Scalar separating the package and the version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $package_href;
    my $package_version_separator;

    my $tmpl = {
        package_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$package_href
        },
        package_version_separator => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$package_version_separator
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @packages;

  PACKAGES:
    while ( my ( $package, $package_version ) = each %{$package_href} ) {
        if ( defined $package_version ) {

            # Check that the version number matches pattern
            if ( $package_version !~ qr/\d+.\d+ | \d+.\d+.\d+/xms ) {
                croak q{The version number does not match defiend pattern for }
                  . q{package: }
                  . $package
                  . $SPACE
                  . $package_version;
            }
            push @packages,
              $package . $package_version_separator . $package_version;
        }
        else {
            push @packages, $package;
        }
    }
    return @packages;
}

sub _create_target_link_paths {

## Function  : Creates paths to bioconda target binaries and links.
##           : Custom solutions for bwakit picard snpeff snpsift manta.
##           : Returns a hash ref consisting of the paths.
## Returns   : %target_link_paths
## Arguments : $bioconda_packages_href => Hash with bioconda packages {REF}
##           : $conda_env_path         => Path to conda environment
##           : $FILEHANDLE             => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten arguments
    my $bioconda_packages_href;
    my $conda_env_path;
    my $FILEHANDLE;

    my $tmpl = {
        bioconda_packages_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$bioconda_packages_href
        },
        conda_env_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$conda_env_path,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use File::Spec::Functions qw{ catfile catdir };
    use MIP::Gnu::Findutils qw{ gnu_find };

    my %target_link_paths;

    my %binaries = (
        bwakit => [
            qw{
              k8                seqtk
              bwa-postalt.js    run-HLA
              typeHLA.sh        fermi2
              fermi2.pl         ropebwt2
              typeHLA-selctg.js typeHLA.js
              }
        ],
        snpeff  => [qw{ snpEff.jar snpEff.config }],
        manta   => [qw{ configManta.py configManta.py.ini }],
        snpsift => [qw{ SnpSift.jar }],
        picard  => [qw{ picard.jar }],
    );

    ## Variables to store the full path in
    my %program_path_aliases = (
        bwakit  => q{BWAKIT_PATH},
        snpeff  => q{SNPEFF_PATH},
        snpsift => q{SNPSIFT_PATH},
        manta   => q{MANTA_PATH},
        picard  => q{PICARD_PATH},
    );

    say {$FILEHANDLE} q{## Find exact path to program and store it for linking};

  PROGRAM:
    foreach my $program ( keys %binaries ) {

        # Check if the program has been set to be installed via shell and
        # thus has been removed from the bioconda_packages hash
        next PROGRAM if ( not $bioconda_packages_href->{$program} );

        ## Capture the full path including the conda patch in a variable
        print {$FILEHANDLE} $program_path_aliases{$program} . q{=} . $BACKTICK;
        my $search_path = catdir( $conda_env_path, q{share},
            $program . q{-} . $bioconda_packages_href->{$program} . q{*} );
        gnu_find(
            {
                search_path   => $search_path,
                test_criteria => q{-type d},
                action        => q{-prune},
                FILEHANDLE    => $FILEHANDLE
            }
        );
        print {$FILEHANDLE} $BACKTICK . $NEWLINE;
        my $program_dir_path = q/${/ . $program_path_aliases{$program} . q/}/;

      BINARY:
        foreach my $binary ( @{ $binaries{$program} } ) {

            ## Construct target path
            my $target_path;
            if ( $program eq q{manta} ) {
                $target_path = catfile( $program_dir_path, q{bin}, $binary );
            }
            else {
                $target_path = catfile( $program_dir_path, $binary );
            }
            ## Construct link_path
            my $link_path = catfile( $conda_env_path, q{bin}, $binary );
            ## Add paths to hash
            $target_link_paths{$target_path} = $link_path;
        }
    }
    return %target_link_paths;
}

1;
