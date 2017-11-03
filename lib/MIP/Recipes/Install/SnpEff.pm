package MIP::Recipes::Install::SnpEff;

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
    our @EXPORT_OK = qw{ install_snpeff check_mt_codon_table };
}

## Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub install_snpeff {

## Function : Install SnpEff
## Returns  : ""
## Arguments: $program_parameters_href => Hash with SnpEff specific parameters {REF}
##          : $conda_prefix_path       => Conda prefix path
##          : $conda_environment       => Conda environment
##          : $noupdate                => Do not update
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity
##          : $FILEHANDLE              => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $snpeff_parameters_href;
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
            store       => \$snpeff_parameters_href
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
    use MIP::Gnu::Coreutils qw{ gnu_rm gnu_mkdir gnu_ln };
    use MIP::Gnu::Findutils qw{ gnu_find };
    use MIP::Program::Download::Wget qw{ wget };
    use MIP::Program::Compression::Zip qw{ unzip };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };
    use MIP::Program::Variantcalling::SnpEff qw{ snpeff_download };
    use MIP::Script::Utils qw{ create_temp_dir };

    ## Unpack parameters
    my $snpeff_version = $snpeff_parameters_href->{version};
    my @snpeff_genome_versions =
      @{ $snpeff_parameters_href->{snpeff_genome_versions} };

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => q{mip_install::install_snpeff},
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    ## Store original working directory
    my $pwd = cwd();

    say {$FILEHANDLE} q{### Install SnpEff};

    my $snpeff_install_path =
      catfile( $conda_prefix_path, qw{ share snpEff }, $snpeff_version );

    # Check if snpEff.jar exists (assumes that snpSift also exists if true)
    if ( -f catfile( $snpeff_install_path, q{snpEff.jar} ) ) {
        $log->info(
            q{SnpEff is already installed in the specified conda environment.});

        if ($noupdate) {
            $log->info(
                q{Skipping writting installation instructions for SnpEffi});
            say {$FILEHANDLE}
              q{## Skipping writting installation instructions for SnpEff};
            say {$FILEHANDLE} $NEWLINE;

            return;
        }

        $log->warn(q{This will overwrite the current SnpEff installation.});

        ## Removing jar files but keeps config and data in order to avoid unnecessary downloads
        say {$FILEHANDLE} q{## Removing old SnpEff jar files};
        gnu_rm(
            {
                infile_path => catfile( $snpeff_install_path, q{*.jar} ),
                force       => 1,
                FILEHANDLE  => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Remove dead symlinks
        say {$FILEHANDLE} q{## Removing old SnpEff links};
        gnu_find(
            {
                search_path   => catdir( $conda_prefix_path, q{bin} ),
                test_criteria => q{-xtype l},
                action        => q{-delete},
                FILEHANDLE    => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    $log->info(q{Writing instructions for SnpEff installation via SHELL});

    ## Creating temporary install directory
    say {$FILEHANDLE} q{## Create temporary SnpEff install directory};
    my $temp_dir = create_temp_dir( { FILEHANDLE => $FILEHANDLE } );
    say {$FILEHANDLE} $NEWLINE;

    ## Download
    say {$FILEHANDLE} q{## Download SnpEff};
    my $url =
        q{http://sourceforge.net/projects/snpeff/files/snpEff}
      . $UNDERSCORE
      . $snpeff_version
      . $UNDERSCORE
      . q{core.zip/download};
    my $snpeff_zip_path = catfile( $temp_dir,
        q{snpeff} . $UNDERSCORE . $snpeff_version . $UNDERSCORE . q{core.zip} );
    wget(
        {
            url          => $url,
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $quiet,
            verbose      => $verbose,
            outfile_path => $snpeff_zip_path
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    unzip(
        {
            infile_path => $snpeff_zip_path,
            outdir_path => $temp_dir,
            quiet       => $quiet,
            verbose     => $verbose,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Create dir for SnpEff in conda env if it doesn't exist
    if ( not -d $snpeff_install_path ) {

        say {$FILEHANDLE} q{## Create dir for SnpEff in conda env};
        gnu_mkdir(
            {
                indirectory_path => $snpeff_install_path,
                parents          => 1,
                FILEHANDLE       => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Make available from conda environment
    say {$FILEHANDLE} q{## Make available from conda environment};
    gnu_mv(
        {
            infile_path  => catfile( $temp_dir, qw{ snpEff *.jar } ),
            outfile_path => $snpeff_install_path,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    print {$FILEHANDLE} $NEWLINE;
    ## Only move config if it isn't already there
    if ( not -f catfile( $snpeff_install_path, q{snpEff.config} ) ) {
        gnu_mv(
            {
                infile_path => catfile( $temp_dir, qw{ snpEff snpEff.config } ),
                outfile_path => $snpeff_install_path,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        print {$FILEHANDLE} $NEWLINE;
    }
    print {$FILEHANDLE} $NEWLINE;

    ## Define binaries
    my @snpeff_binaries = qw(snpEff.jar SnpSift.jar snpEff.config);

    foreach my $binary (@snpeff_binaries) {
        ## Specifying target and link path
        my $target_path = catfile( $snpeff_install_path, $binary );
        my $link_path = catfile( $conda_prefix_path, q{bin}, $binary );
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

  GENOME_VERSION:
    foreach my $genome_version (@snpeff_genome_versions) {
        ## Check and if required add the vertebrate mitochondrial codon table to SnpEff config
        check_mt_codon_table(
            {
                FILEHANDLE     => $FILEHANDLE,
                share_dir      => $snpeff_install_path,
                config_file    => q{snpEff.config},
                genome_version => $genome_version,
                quiet          => $quiet,
                verbose        => $verbose,
            }
        );

        next GENOME_VERSION
          if ( -d catdir( $snpeff_install_path, q{data}, $genome_version ) );
        ## Write instructions to download SnpEff database.
        ## This is done by install script to avoid race conditin when doing first analysis run in MIP

        say {$FILEHANDLE} q{## Downloading SnpEff database};
        my $jar_path = catfile( $conda_prefix_path, qw{ bin snpEff.jar} );
        my $config_file_path =
          catfile( $conda_prefix_path, qw{bin snpEff.config} );
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

    ## Remove the temporary install directory
    say {$FILEHANDLE} q{## Remove temporary install directory};
    gnu_rm(
        {
            infile_path => $temp_dir,
            recursive   => 1,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE x 2;

    return;
}

sub check_mt_codon_table {

##Function : Check and if required add the vertebrate mitochondrial codon table to snpeff config
##Returns  : ""
##Arguments: $FILEHANDLE     => FILEHANDLE to write to
##         : $share_dir      => The conda env shared directory
##         : $config_file    => The config config_file
##         : $genome_version => snpeff genome version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $share_dir;
    my $config_file;
    my $genome_version;
    my $verbose;
    my $quiet;

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

    use File::Spec::Functions qw{ catfile };
    use IPC::Cmd qw{ can_run run };
    use MIP::Gnu::Coreutils qw{ gnu_mv };

    ## Get logger
    my $log = retrieve_log(
        {
            log_name => q{mip_install::install_snpeff},
            verbose  => $verbose,
            quiet    => $quiet,
        }
    );

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

      # Print matching element
      . q?print $_;? . $SPACE

      # Print MT codon table that is being downloaded
      . q?print "?
      . $genome_version
      . q?.MT.codonTable : Vertebrate_Mitochondrial\n"}?
      . $SPACE

      # If no match: print line
      . q?else {print $_;}' ?;

    ## Test if config file contains the desired MT codon table
    my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf );
    if ( -f catfile( $share_dir, $config_file ) ) {
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run(
            command => $detect_regexp . catfile( $share_dir, $config_file ) );
    }

    #No MT.codonTable in config
    if ( not $success ) {
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

        my $infile_path = catfile( $share_dir, $config_file . $DOT . q{tmp} );
        gnu_mv(
            {
                infile_path  => $infile_path,
                outfile_path => catfile( $share_dir, $config_file ),
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

    }
    else {
        $log->info( q{Found MT.codonTable for}
              . $SPACE
              . $genome_version
              . $SPACE
              . q{in snpEff.config. Skipping addition to snpEff.config} );
    }
    return;
}

1;
