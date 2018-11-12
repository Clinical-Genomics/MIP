package MIP::Recipes::Install::Conda;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;

## CPANM
use Readonly;

## Constants
Readonly my $BACKTICK   => q{`};
Readonly my $COLON      => q{:};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $PIPE       => q{|};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ check_conda_installation get_conda_dir_path install_conda_packages };

}

sub check_conda_installation {

## Function  : Checks that conda is installed sets conda directory path
## Returns   :
## Arguments : $conda_dir_path    => Path to conda environment (default: conda root)
##           : $disable_env_check => Disable environment check
##           : $parameter_href    => The entire parameter hash {REF}
##           : $quiet             => Log only warnings and above
##           : $verbose           => Log debug messages

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_dir_path;
    my $disable_env_check;
    my $parameter_href;
    my $quiet;
    my $verbose;

    my $tmpl = {
        conda_dir_path => {
            required    => 1,
            store       => \$conda_dir_path,
            strict_type => 1,
        },
        disable_env_check => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$disable_env_check,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
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

    use File::Spec::Functions qw{ catdir splitdir };
    use IPC::Cmd qw{ run };
    use MIP::Check::Unix qw{ check_binary_in_path };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };
    use MIP::Package_manager::Conda qw{ conda_check_env_status };

    ## Get logger
    my $log = retrieve_log(
        {
            log_name => q{mip_install::check_conda_installation},
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    ## Scan the PATH for conda
    check_binary_in_path(
        {
            binary       => q{conda},
            log          => $log,
            program_name => q{conda},
        }
    );

    ## Establish path to conda, Exit if not found
    $parameter_href->{conda_dir_path} = get_conda_dir_path(
        {
            conda_dir_path => $conda_dir_path,
            log            => $log,
        }
    );

    ## Check for active conda environment (exit if true)
    conda_check_env_status(
        {
            disable_env_check => $disable_env_check,
            log               => $log
        }
    );

    return;
}

sub install_conda_packages {

## Function  : Install conda packages into a new or existing conda environment.
## Returns   :
## Arguments : $conda_env                  => Name of conda environment
##           : $conda_env_path             => Path to conda environment (default: conda root)
##           : $conda_no_update_dep        => Do not update dependencies
##           : $conda_packages_href        => Hash holding conda packages and their version numbers {REF}
##           : $conda_update               => Update Conda if defined
##           : $FILEHANDLE                 => Filehandle to write to
##           : $snpeff_genome_versions_ref => Array with the genome versions for the snpeff databases {REF}
##           : $quiet                      => Log only warnings and above
##           : $verbose                    => Log debug messages

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_env;
    my $conda_env_path;
    my $conda_packages_href;
    my $conda_update;
    my $FILEHANDLE;
    my $quiet;
    my $snpeff_genome_versions_ref;
    my $verbose;

    ## Defaults
    my $conda_no_update_dep;

    my $tmpl = {
        conda_env => {
            strict_type => 1,
            store       => \$conda_env,
        },
        conda_env_path => {
            defined     => 1,
            required    => 1,
            store       => \$conda_env_path,
            strict_type => 1,
        },
        conda_no_update_dep => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$conda_no_update_dep,
            strict_type => 1,
        },
        conda_packages_href => {
            default     => {},
            required    => 1,
            store       => \$conda_packages_href,
            strict_type => 1,
        },
        conda_update => {
            store => \$conda_update,
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE,
        },
        quiet => {
            allow => [ undef, 0, 1 ],
            store => \$quiet,
        },
        snpeff_genome_versions_ref => {
            default  => [],
            required => 1,
            store    => \$snpeff_genome_versions_ref,
        },
        verbose => {
            allow => [ undef, 0, 1 ],
            store => \$verbose,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments};

    use Array::Utils qw{ intersect };
    use MIP::Gnu::Bash qw{ gnu_unset };
    use MIP::Gnu::Coreutils qw{ gnu_ln };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };
    use MIP::Package_manager::Conda qw{ conda_create conda_install };

    ## Packages to be installed
    my @conda_packages = keys %{$conda_packages_href};

    ## Return if no packages are to be installed
    return if not @conda_packages;

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => q{mip_install::install_conda_packages},
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

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

    ## Create an array for conda packages that are to be installed from provided hash
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
                q{Writing installation instructions for environment: } . $conda_env );
            say {$FILEHANDLE} q{## Creating conda environment: }
              . $conda_env
              . q{ and install packages};
            conda_create(
                {
                    conda_channels_ref => [qw{ bioconda conda-forge }],
                    env_name           => $conda_env,
                    FILEHANDLE         => $FILEHANDLE,
                    packages_ref       => \@packages,
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }
        else {

            $log->warn(
                q{Conda environment: } . $conda_env . $SPACE . q{already exists} );
            $log->warn(q{Will try to install packages into existing environment});
            $log->info(
                q{Writing installation instructions for conda packages to environment: }
                  . $conda_env );
            say {$FILEHANDLE} q{## Installing conda packages into existing environment};
            conda_install(
                {
                    conda_channels_ref => [qw{ bioconda conda-forge }],
                    no_update_dep      => $conda_no_update_dep,
                    FILEHANDLE         => $FILEHANDLE,
                    env_name           => $conda_env,
                    packages_ref       => \@packages,
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
                conda_channels_ref => [qw{ bioconda conda-forge }],
                FILEHANDLE         => $FILEHANDLE,
                packages_ref       => \@packages,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Linking and custom solutions
    my @custom_solutions = qw{ bwakit | gatk | manta | picard | snpeff | snpsift };

    ## Link conda packages
    # Creating target-link paths
    my %target_link_paths = _create_target_link_paths(
        {
            conda_env_path       => $conda_env_path,
            conda_packages_href  => $conda_packages_href,
            custom_solutions_ref => \@custom_solutions,
            FILEHANDLE           => $FILEHANDLE,
        }
    );

    if (%target_link_paths) {
        say {$FILEHANDLE} q{## Creating symbolic links for conda packages};
      TARGET_AND_LINK_PATHS:
        while ( my ( $target_path, $link_path ) = each %target_link_paths ) {
            gnu_ln(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    force       => 1,
                    link_path   => $link_path,
                    symbolic    => 1,
                    target_path => $target_path,
                }
            );
            print {$FILEHANDLE} $NEWLINE;
        }
        print {$FILEHANDLE} $NEWLINE;
    }

    ## Custom solutions for BWA, SnpEff, Manta and GATK
    ## Copying files, downloading necessary databases and make files executable
    finish_conda_package_install(
        {
            conda_env                  => $conda_env,
            conda_env_path             => $conda_env_path,
            conda_packages_href        => $conda_packages_href,
            custom_solutions_ref       => \@custom_solutions,
            FILEHANDLE                 => $FILEHANDLE,
            log                        => $log,
            quiet                      => $quiet,
            snpeff_genome_versions_ref => $snpeff_genome_versions_ref,
            verbose                    => $verbose,
        }
    );

    ## Unset variables

    if ( intersect( @custom_solutions, @conda_packages ) ) {
        say {$FILEHANDLE} q{## Unset variables};
        my %program_path_aliases = (
            bwakit  => q{BWAKIT_PATH},
            manta   => q{MANTA_PATH},
            picard  => q{PICARD_PATH},
            snpeff  => q{SNPEFF_PATH},
            snpsift => q{SNPSIFT_PATH},
        );

      PROGRAM:
        foreach my $program ( keys %program_path_aliases ) {

            # Check if the program has been set to be installed via shell and
            # thus has been removed from the conda_packages hash
            next PROGRAM if ( not $conda_packages_href->{$program} );

            gnu_unset(
                {
                    bash_variable => $program_path_aliases{$program},
                    FILEHANDLE    => $FILEHANDLE,
                }
            );
            print {$FILEHANDLE} $NEWLINE;
        }
        say {$FILEHANDLE} $NEWLINE;
    }

    return;
}

sub finish_conda_package_install {

## Function  : Custom solutions to finish the install of BWA, SnpEff, Manta and GATK
## Returns   :
## Arguments : $conda_env                  => Name of conda env
##           : $conda_env_path             => Path to conda environment
##           : $conda_packages_href        => Hash with conda packages {REF}
##           : $custom_solutions_ref       => Regex with programs that requires some fiddling
##           : $FILEHANDLE                 => Filehandle to write to
##           : $log                        => Log
##           : $quiet                      => Log only warnings and above
##           : $snpeff_genome_versions_ref => Array with the genome versions for the snpeff databases {REF}
##           : $verbose                    => Log debug messages

    my ($arg_href) = @_;

    ## Flatten arguments
    my $conda_env;
    my $conda_env_path;
    my $conda_packages_href;
    my $custom_solutions_ref;
    my $FILEHANDLE;
    my $log;
    my $quiet;
    my $snpeff_genome_versions_ref;
    my $verbose;

    my $tmpl = {
        conda_env => {
            required => 1,
            store    => \$conda_env,
        },
        conda_env_path => {
            defined     => 1,
            required    => 1,
            store       => \$conda_env_path,
            strict_type => 1,
        },
        conda_packages_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$conda_packages_href,
            strict_type => 1,
        },
        custom_solutions_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$custom_solutions_ref,
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
        quiet => {
            allow => [ undef, 0, 1 ],
            store => \$quiet,
        },
        snpeff_genome_versions_ref => {
            default  => [],
            required => 1,
            store    => \$snpeff_genome_versions_ref,
        },
        verbose => {
            allow => [ undef, 0, 1 ],
            store => \$verbose,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Array::Utils qw{ intersect };
    use File::Spec::Functions qw{ catdir catfile };
    use IPC::Cmd qw{ run };
    use MIP::Gnu::Coreutils qw{ gnu_cp gnu_chmod gnu_rm };
    use MIP::Package_manager::Conda qw{ conda_activate conda_deactivate };
    use MIP::Program::Variantcalling::Snpeff qw{ snpeff_download };
    use MIP::Recipes::Install::Gatk qw{ gatk_download };
    use MIP::Recipes::Install::SnpEff qw{ check_mt_codon_table };

    my @conda_packages = keys %{$conda_packages_href};

    ## Return if no custom solutions are required
    return if not intersect( @{$custom_solutions_ref}, @conda_packages );

    ## Only activate conda environment if supplied by user
    if ($conda_env) {

        ## Activate conda environment
        say {$FILEHANDLE} q{## Activate conda environment};
        conda_activate(
            {
                env_name   => $conda_env,
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Custom BWA
    ## Check if bwakit has been removed from conda installation hash
    if ( $conda_packages_href->{bwakit} ) {

        say {$FILEHANDLE} q{## Custom BWA solutions};

        ## Double quote to avoid expansion of shell variable
        my $infile_path  = catdir( q/"${BWAKIT_PATH}"/, q{resource-human-HLA} );
        my $outfile_path = catdir( $conda_env_path,     q{bin} );

        gnu_cp(
            {
                FILEHANDLE   => $FILEHANDLE,
                force        => 1,
                infile_path  => $infile_path,
                outfile_path => $outfile_path,
                recursive    => 1,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Custom SnpEff
    ## Check if snpeff has been set to be installed via shell or excluded from installation
    if ( $conda_packages_href->{snpeff} ) {

        ## Custom snpeff - Download necessary databases
        ## Check and if required add the vertebrate mitochondrial codon table to snpeff config
        say {$FILEHANDLE} q{## Custom SnpEff solutions};

        ## Get the full version, including patch for snpeff
        my $version = _get_full_snpeff_version(
            {
                log            => $log,
                snpeff_version => $conda_packages_href->{snpeff},
            }
        );

      SNPEFF_GENOME_VERSION:
        foreach my $genome_version ( @{$snpeff_genome_versions_ref} ) {

            my $share_dir = catdir( $conda_env_path, q{share}, q{snpeff-} . $version );

            check_mt_codon_table(
                {
                    config_file    => q{snpEff.config},
                    FILEHANDLE     => $FILEHANDLE,
                    genome_version => $genome_version,
                    quiet          => $quiet,
                    share_dir      => $share_dir,
                    verbose        => $verbose,
                }
            );

            my $snpeff_genome_dir = catdir( $share_dir, q{data}, $genome_version );
            next SNPEFF_GENOME_VERSION if ( -d $snpeff_genome_dir );

            ## Write instructions to download snpeff database.
            ## This is done by install script to avoid race conditin when doing first analysis run in MIP
            say {$FILEHANDLE} q{## Downloading snpeff database};
            my $jar_path         = catfile( $conda_env_path, qw{ bin snpEff.jar} );
            my $config_file_path = catfile( $conda_env_path, qw{bin snpEff.config} );

            snpeff_download(
                {
                    config_file_path        => $config_file_path,
                    FILEHANDLE              => $FILEHANDLE,
                    genome_version_database => $genome_version,
                    jar_path                => $jar_path,
                    temp_directory          => 1,
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }
    }

    ## Custom manta
    ## Check if manta has been removed from conda installation hash
    if ( $conda_packages_href->{manta} ) {

        # Make file executable
        say {$FILEHANDLE} q{## Changing mode of configManta.py to executable};
        my $file_path = catfile( $conda_env_path, qw{bin configManta.py} );
        gnu_chmod(
            {
                FILEHANDLE => $FILEHANDLE,
                file_path  => $file_path,
                permission => q{a+x},
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Custom GATK
    ## Check if GATK has been removed from conda installation hash
    if ( $conda_packages_href->{gatk} ) {
        say {$FILEHANDLE} q{## Custom GATK solutions};

        ## Download gatk .tar.bz2
        my $gatk_tar_path = gatk_download(
            {
                FILEHANDLE   => $FILEHANDLE,
                gatk_version => $conda_packages_href->{gatk},
                quiet        => $quiet,
                verbose      => $verbose,
            }
        );

        ## Hard coding here since GATK 4.0 will be open source.
        ## Then this step will be unnecessary
        say {$FILEHANDLE} q{gatk3-register} . $SPACE . $gatk_tar_path . $NEWLINE;

        gnu_rm(
            {
                FILEHANDLE  => $FILEHANDLE,
                force       => 1,
                infile_path => dirname($gatk_tar_path),
                recursive   => 1,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Remove of /tmp/gatk from gatk 3.8 installation
        my $tmpdir = File::Spec->tmpdir();
        gnu_rm(
            {
                FILEHANDLE  => $FILEHANDLE,
                force       => 1,
                infile_path => catdir( $tmpdir, q{gatk} ),
                recursive   => 1,
            }
        );
        say {$FILEHANDLE} $NEWLINE x 2;
    }

    ## Deactivate conda environment if conda_environment exists
    if ($conda_env) {

        say {$FILEHANDLE} q{## Deactivate conda environment};
        conda_deactivate(
            {
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }
    return;
}

sub _create_package_array {

## Function  : Takes a reference to hash of packages and creates an array with
##           : package and version joined with a supplied separator if value is defined.
##           : Also checks that the version number makes sense
## Returns   : "@packages"
## Arguments : $package_href              => Hash with packages {Hash}
##           : $package_version_separator => Scalar separating the package and the version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $package_href;
    my $package_version_separator;

    my $tmpl = {
        package_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$package_href,
            strict_type => 1,
        },
        package_version_separator => {
            defined     => 1,
            required    => 1,
            store       => \$package_version_separator,
            strict_type => 1,
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
            push @packages, $package . $package_version_separator . $package_version;
        }
        else {
            push @packages, $package;
        }
    }
    return @packages;
}

sub _create_target_link_paths {

## Function  : Creates paths to conda target binaries and links.
##           : Custom solutions for bwakit picard snpeff snpsift manta.
##           : Returns a hash ref consisting of the paths.
## Returns   : %target_link_paths
## Arguments : $conda_packages_href  => Hash with conda packages {REF}
##           : $conda_env_path       => Path to conda environment
##           : $custom_solutions_ref => Array with programs that requires som fiddling {REF}
##           : $FILEHANDLE           => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten arguments
    my $conda_packages_href;
    my $conda_env_path;
    my $custom_solutions_ref;
    my $FILEHANDLE;

    my $tmpl = {
        conda_packages_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$conda_packages_href,
            strict_type => 1,
        },
        conda_env_path => {
            defined     => 1,
            required    => 1,
            store       => \$conda_env_path,
            strict_type => 1,
        },
        custom_solutions_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$custom_solutions_ref,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Array::Utils qw{ intersect };
    use File::Spec::Functions qw{ catfile catdir };
    use MIP::Gnu::Coreutils qw{ gnu_tail };
    use MIP::Gnu::Findutils qw{ gnu_find };

    my @conda_packages = keys %{$conda_packages_href};

    ## Skip if no program requires linking
    return if ( not intersect( @{$custom_solutions_ref}, @conda_packages ) );

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
        manta   => [qw{ configManta.py configManta.py.ini }],
        picard  => [qw{ picard.jar }],
        snpeff  => [qw{ snpEff.jar snpEff.config }],
        snpsift => [qw{ SnpSift.jar }],
    );

    ## Variables to store the full path in
    my %program_path_aliases = (
        bwakit  => q{BWAKIT_PATH},
        manta   => q{MANTA_PATH},
        picard  => q{PICARD_PATH},
        snpeff  => q{SNPEFF_PATH},
        snpsift => q{SNPSIFT_PATH},
    );

    say {$FILEHANDLE} q{## Find exact path to program and store it for linking};

  PROGRAM:
    foreach my $program ( keys %binaries ) {

        # Check if the program has been set to be installed via shell and
        # thus has been removed from the conda_packages hash
        next PROGRAM if ( not $conda_packages_href->{$program} );

        ## Capture the full path including the conda patch in a variable

        # Alias for translation
        my $conda_version = $conda_packages_href->{$program};
        ## Exchange for conda internal format when using full conda version
        ## i.e. version=subpatch
        $conda_version =~ tr/=/-/;

        ## Special case for Manta
        if ( $program eq q{manta} ) {

            $conda_version =~ s/py27_//g;
        }

        print {$FILEHANDLE} $program_path_aliases{$program} . q{=} . $BACKTICK;
        my $search_path =
          catdir( $conda_env_path, q{share}, $program . q{-} . $conda_version . q{*} );
        gnu_find(
            {
                action        => q{-prune},
                FILEHANDLE    => $FILEHANDLE,
                search_path   => $search_path,
                test_criteria => q{-type d},
            }
        );
        ## Pipe to next command
        print {$FILEHANDLE} $PIPE . $SPACE;
        ## Only use the latest latest sub patch
        gnu_tail(
            {
                FILEHANDLE => $FILEHANDLE,
                lines      => q{1},
            }
        );
        say {$FILEHANDLE} $BACKTICK;
        ## Double quotes to avoid expansion in shell
        my $program_dir_path = q/"${/ . $program_path_aliases{$program} . q/}"/;

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

sub get_conda_dir_path {

## Function : Finds the conda directory path and returns it
## Returns  : $conda_dir_path
## Arguments: $conda_dir_path => Path to conda dir
##          : $log            => Log

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_dir_path;
    my $log;

    my $tmpl = {
        conda_dir_path => {
            strict_type => 1,
            store       => \$conda_dir_path,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Contstants
    Readonly my $REMOVE_BIN_CONDA_DIR => -2;

    ## Establish path to conda
    if ( $conda_dir_path and not -d $conda_dir_path ) {

        $log->fatal( q{Could not find miniconda directory in}
              . $COLON
              . $SPACE
              . $conda_dir_path );
        exit 1;
    }
    elsif ( not $conda_dir_path ) {
        my $conda_bin_in_path;
        run(
            command => q{which conda},
            buffer  => \$conda_bin_in_path,
        );

        ## If no output from which command
        if ( not $conda_bin_in_path ) {

            $log->fatal( q{Could not find conda installation path.}
                  . $NEWLINE
                  . q{Specify the path to conda using the --conda_dir_path flag} );
            exit 1;
        }
        ## Get the path to the conda base directory
        my @conda_bin_in_paths = splitdir($conda_bin_in_path);

        # Remove the two last elements of the array (bin and conda)
        splice @conda_bin_in_paths, $REMOVE_BIN_CONDA_DIR;
        $conda_dir_path = catdir(@conda_bin_in_paths);
    }

    return $conda_dir_path;
}

sub _get_full_snpeff_version {

## Function  : Get the snpeff version together with patch that are to be installed via Conda
## Returns   : $version
## Arguments : $log             => Log
##           : $snpeff_version  => Main snpeff version number

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $snpeff_version;
    my $log;

    my $tmpl = {
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        snpeff_version => {
            defined     => 1,
            required    => 1,
            store       => \$snpeff_version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use List::MoreUtils qw{ first_index };

    ## Store shell query
    my $command =

      # Get the snpeff version that will be installed
      qq{conda search -c bioconda -c conda-forge --spec snpeff=$snpeff_version}

      # Isolate the version with the latest sub patch
      . q{ | tail -1 };

    ## Capture output in string
    my $snpeff_search_output;

    run(
        command => $command,
        buffer  => \$snpeff_search_output,
    );

    my @snpeff_array = split $SPACE, $snpeff_search_output;

    ## Get the index of the element matching the snpeff version
    my $version_index = first_index { m/$snpeff_version/xms } @snpeff_array;

    ## Concatenate the version and sub patch based on their indexes in the array
    my $version =
      $snpeff_array[$version_index] . q{-} . $snpeff_array[ $version_index + 1 ];

    ## Something went wrong
    if ( not $version ) {

        $log->logcroak(
            q{Could not find the snpeff version and subpatch from string: }
              . $snpeff_search_output );
    }
    return $version;
}

1;
