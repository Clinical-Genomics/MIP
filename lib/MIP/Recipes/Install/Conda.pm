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
use Readonly;

## Constants
Readonly my $SPACE   => q{ };
Readonly my $NEWLINE => qq{\n};
Readonly my $DOT     => q{.};

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.0.0;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ setup_conda_env install_bioconda_packages finish_bioconda_package_install };

}

## Setup Conda env and install default packages

sub setup_conda_env {

## setup_conda_env

## Function  : Creates necessary conda environment and install package(s) from the default channel.
## Returns   :
## Arguments : $conda_packages_href, $conda_env, $conda_env_path, $FILEHANDLE, $conda_update
##           : $conda_packages_href => Hash with conda packages and their version numbers {REF}
##           : $conda_env           => Name of conda environment
##           : $conda_env_path      => Path to conda environment (could bee path to root env)
##           : $FILEHANDLE          => Filehandle to wrire to
##           : $conda_update        => Update Conda if defined

    my ($arg_href) = @_;

    ## Defaults

    ## Flatten argument(s)
    my $conda_packages_href;
    my $conda_env;
    my $conda_env_path;
    my $FILEHANDLE;
    my $conda_update;

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
        }
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Unix qw{ check_binary_in_path };
    use MIP::PacketManager::Conda
      qw{ conda_create conda_update conda_check conda_install };

    ## Scan the PATH for conda
    check_binary_in_path(
        {
            binary => q{conda},
        }
    );

    ## Check for active conda environment (exit if true)
    conda_check();

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
        if ( !-d $conda_env_path ) {
            ## Create conda environment and install packages
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
            say STDERR q{Conda environment: } . $conda_env . q{ already exists};
            say STDERR q{Will try to install packages into existing envronment};
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

## install_bioconda_packages

## Function  : Install conda packages from the bioconda channel into a conda environment.
## Returns   :
## Arguments : $bioconda_packages_href, $bioconda_patches_href, $conda_env, $conda_env_path, $FILEHANDLE
##           : $bioconda_packages_href => Hash holding bioconda packages and their version numbers {REF}
##           : $bioconda_patches_href  => Hash holding the patches for the bioconda packages {REF}
##           : $conda_env              => Name of conda environment
##           : $conda_env_path         => Path to conda environment (could bee path to root env)
##           : $FILEHANDLE             => Filehandle to write to

    my ($arg_href) = @_;

    ## Defaults

    ## Flatten argument(s)
    my $bioconda_packages_href;
    my $bioconda_patches_href;
    my $conda_env;
    my $conda_env_path;
    my $FILEHANDLE;

    my $tmpl = {
        bioconda_packages_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$bioconda_packages_href,
        },
        bioconda_patches_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$bioconda_patches_href,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments};

    use MIP::PacketManager::Conda qw{ conda_install };
    use MIP::Gnu::Coreutils qw{ gnu_ln };

    ## Create an array for bioconda packages that are to be installed from provided hash
    my @packages = _create_package_array(
        {
            package_href              => $bioconda_packages_href,
            package_version_separator => q{=},
        }
    );

    ## Install bioconda packages
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
            bioconda_patches_href  => $bioconda_patches_href,
            conda_env_path         => $conda_env_path,
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

    return;
}

sub finish_bioconda_package_install {

## finish_bioconda_package_install

## Function   : Custom solutions to finish the install of BWA, SnpEff and Manta
## Returns    :
## Argumemnts : $bioconda_packages_href, $bioconda_patches_href, $conda_env_path, $FILEHANDLE
##            : $bioconda_packages_href     => Hash with bioconda packages {REF}
##            : $bioconda_patches_href      => Hash with package patches {REF}
##            : $snpeff_genome_versions_ref => Hash with the genome versins of the snpeff databases {REF}
##            : $conda_env_path             => Path to conda environment
##            : $FILEHANDLE                 => Filehandle to write to
##            : $conda_env                  => Name of conda env

    my ($arg_href) = @_;

    ## Flatten arguments
    my $bioconda_packages_href;
    my $bioconda_patches_href;
    my $snpeff_genome_versions_ref;
    my $conda_env_path;
    my $FILEHANDLE;
    my $conda_env;

    my $tmpl = {
        bioconda_packages_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$bioconda_packages_href
        },
        bioconda_patches_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$bioconda_patches_href
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use File::Spec::Functions qw{ catdir catfile };
    use MIP::Gnu::Coreutils qw{ gnu_cp gnu_chmod };

    ## Custom BWA
    say {$FILEHANDLE} q{## Custom BWA solutions};
    my $infile_path = catdir(
        $conda_env_path,
        q{share},
        q{bwakit-}
          . $bioconda_packages_href->{bwakit}
          . $bioconda_patches_href->{bioconda_bwakit_patch},
        q{resource-human-HLA}
    );
    my $outfile_path = catdir( $conda_env_path, q{bin} );
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

    ## Custom snpeff - Download necessary databases
    ## Check and if required add the vertebrate mitochondrial codon table to snpeff config
    say {$FILEHANDLE} q{## Custom SnpEff solutions};
  SNPEFF_GENOME_VERSIONS:
    foreach my $genome_version ( @{$snpeff_genome_versions_ref} ) {
        my $share_dir = catdir( $conda_env_path, q{share},
                q{snpeff-}
              . $bioconda_packages_href->{snpeff}
              . $bioconda_patches_href->{bioconda_snpeff_patch} );
        _check_mt_codon_table(
            {
                FILEHANDLE     => $FILEHANDLE,
                share_dir      => $share_dir,
                config_file    => q{snpEff.config},
                genome_version => $genome_version,
            }
        );
        if (
            !-d catdir(
                $conda_env_path,
                q{share},
                q{snpeff-}
                  . $bioconda_packages_href->{snpeff}
                  . $bioconda_patches_href->{bioconda_snpeff_patch},
                q{data},
                $genome_version
            )
          )
        {
            ## Write instructions to download snpeff database.
            ## This is done by install script to avoid race conditin when doing first analysis run in MIP
            say {$FILEHANDLE} q{## Downloading snpeff database};
            _snpeff_download(
                {
                    FILEHANDLE     => $FILEHANDLE,
                    genome_version => $genome_version,
                    conda_env      => $conda_env,
                    conda_env_path => $conda_env_path,
                }
            );
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

    return;
}

sub _create_package_array {

## create_package_list

##Function  : Takes a reference to hash of packages and creates an array with
##          : package and version joined with a supplied separator if value is defined.
##          : Also checks that the version number makes sense
##Returns   : "@packages"
##Arguments : $package_href, $package_version_separator
##          : $package_href              => Hash with packages {Hash}
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

## _create_target_link_paths

## Function   : Creates paths to bioconda target binaries and links.
##            : Custom solutions for bwakit picard snpeff snpsift manta.
##            : Returns a hash ref consisting of the paths.
## Returns    : %target_link_paths
## Argumemnts : $bioconda_packages_href, $bioconda_patches_href, $conda_env_path
##            : $bioconda_packages_href => Hash with bioconda packages {REF}
##            : $bioconda_patches_href  => Hash with bioconda package patches {REF}
##            : $conda_env_path         => Path to conda environment

    my ($arg_href) = @_;

    ## Flatten arguments
    my $bioconda_packages_href;
    my $bioconda_patches_href;
    my $conda_env_path;

    my $tmpl = {
        bioconda_packages_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$bioconda_packages_href
        },
        bioconda_patches_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$bioconda_patches_href
        },
        conda_env_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$conda_env_path,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use File::Spec::Functions qw{ catfile };

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

  PROGRAMS:
    foreach my $program ( keys %binaries ) {

      BINARIES:
        foreach my $binary ( @{ $binaries{$program} } ) {
            my $target_path;
            ## Construct target path
            if ( $program eq q{manta} ) {
                $target_path = catfile(
                    $conda_env_path,
                    q{share},
                    $program . q{-}
                      . $bioconda_packages_href->{$program}
                      . $bioconda_patches_href->{ q{bioconda_}
                          . $program
                          . q{_patch} },
                    q{bin},
                    $binary
                );
            }
            else {
                $target_path = catfile(
                    $conda_env_path,
                    q{share},
                    $program . q{-}
                      . $bioconda_packages_href->{$program}
                      . $bioconda_patches_href->{ q{bioconda_}
                          . $program
                          . q{_patch} },
                    $binary
                );
            }
            ## Construct link_path
            my $link_path = catfile( $conda_env_path, q{bin}, $binary );
            ## Add paths to hash
            $target_link_paths{$target_path} = $link_path;
        }
    }
    return %target_link_paths;
}

sub _check_mt_codon_table {

##_check_mt_codon_table

##Function : Check and if required add the vertebrate mitochondrial codon table to snpeff config
##Returns  : ""
##Arguments: $FILEHANDLE, $share_dir, $config_file, $genome_version
##         : $FILEHANDLE     => FILEHANDLE to write to
##         : $share_dir      => The conda env shared directory
##         : $config_file    => The config config_file
##         : $genome_version => snpeff genome version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $share_dir;
    my $config_file;
    my $genome_version;

    my $tmpl = {
        FILEHANDLE => {
            required => 1,
            defined  => 1,
            store    => \$FILEHANDLE
        },
        share_dir => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$share_dir
        },
        config_file => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$config_file
        },
        genome_version => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$genome_version
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use File::Spec::Functions qw{ catfile };
    use IPC::Cmd qw{ can_run run };
    use MIP::Gnu::Coreutils qw{ gnu_mv };

    my $detect_regexp =

      # Execute perl, loop over input and split on whitespace
      q?perl -nae? . $SPACE

      # Check for a MT codon table with matching genome version
      . q?'if($_=~/? . $genome_version . q?.MT.codonTable/)?

      # Print 1 in case of match
      . q?{print 1}' ?;

    my $add_regexp =

      # Execute perl, loop over input and split on whitespace
      q?perl -nae? . $SPACE

      # Search for  genome version .reference match
      . q?'if($_=~/? . $genome_version . q?.reference/) {?

      # print matching element
      . q?print $_;? . $SPACE

      # print MT codon table that is being downloaded
      . q?print "?
      . $genome_version
      . q?.MT.codonTable : Vertebrate_Mitochondrial\n"}?
      . $SPACE

      # if no match: print line
      . q?else {print $_;}' ?;

    my @ret =
      my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf );

    if ( -f catfile( $share_dir, $config_file ) ) {
        @ret = run( command => ( $detect_regexp, $share_dir / $config_file ) );
    }

    if ( !$success ) {    #No MT.codonTable in config
        say {$FILEHANDLE} q{## Adding }
          . $genome_version
          . q{.MT.codonTable : Vertebrate_Mitochondrial to }
          . $share_dir
          . $config_file;

        ## Add MT.codon Table to config
        say {$FILEHANDLE} $add_regexp
          . $SPACE
          . catfile( $share_dir, $config_file ) . q{ > }
          . catfile( $share_dir, $config_file . $DOT . q{tmp} );
        gnu_mv(
            {
                infile_path =>
                  catfile( $share_dir, $config_file . $DOT . q{tmp} ),
                outfile_path => catfile( $share_dir, $config_file ),
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

    }
    else {
        say STDERR q{Found MT.codonTable in}
          . $SPACE
          . catfile( $share_dir, q{snpEff.config} )
          . $DOT
          . $SPACE
          . q{Skipping addition to snpEff config};
    }
    return;
}

sub _snpeff_download {

##_snpeff_download

##Function : Write instructions to download snpeff database. This is done by install script to avoid race condition when doing first analysis run in MIP
##Returns  : ""
##Arguments: $FILEHANDLE, $genome_version, $conda_env, $conda_env_path
##         : $FILEHANDLE      => FILEHANDLE to write to
##         : $genome_version  => snpeff genome version
##         : $conda_env       => name of conda env
##         : $coonda_env_path => path to conda env

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $genome_version;
    my $conda_env;
    my $conda_env_path;

    my $tmpl = {
        FILEHANDLE => {
            required => 1,
            defined  => 1,
            store    => \$FILEHANDLE
        },
        genome_version => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$genome_version
        },
        conda_env => {
            required => 1,
            store    => \$conda_env
        },
        conda_env_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$conda_env_path,
        },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::PacketManager::Conda
      qw { conda_source_activate conda_source_deactivate };
    use MIP::Language::Java qw{ java_core };
    use File::Spec::Functions qw{ catfile };
    use MIP::Unix::Write_to_file qw{ unix_write_to_file };

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

    ## Build base java command
    my $java_jar = catfile( $conda_env_path, qw{ bin snpEff.jar} );
    my @commands = java_core(
        {
            memory_allocation => q{Xmx2g},
            java_jar          => $java_jar,
        }
    );

    ## Build rest of java command
    push @commands, qw{ download -v };
    push @commands, $genome_version, q{-c};
    push @commands, catfile( $conda_env_path, qw{bin snpEff.config} );

    ## Write rest of java commadn to $FILEHANDLE
    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    say {$FILEHANDLE} $NEWLINE;

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

1;
