package MIP::Recipes::Install::SnpEff;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## CPAN
use autodie qw{ :all };
use File::Spec::Functions qw{ catfile };
use IPC::Cmd qw{ can_run run };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $DOT $LOG_NAME $NEWLINE $SPACE $UNDERSCORE };
use MIP::Gnu::Coreutils qw{ gnu_ln gnu_mkdir gnu_mv gnu_rm };
use MIP::Gnu::Findutils qw{ gnu_find };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Program::Compression::Zip qw{ unzip };
use MIP::Program::Download::Wget qw{ wget };
use MIP::Program::Variantcalling::Snpeff qw{ snpeff_download };
use MIP::Script::Utils qw{ create_temp_dir };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_snpeff check_mt_codon_table };
}

sub install_snpeff {

## Function : Install SnpEff
## Returns  :
## Arguments: $conda_environment       => Conda environment
##          : $conda_prefix_path       => Conda prefix path
##          : $FILEHANDLE              => Filehandle to write to
##          : $program_parameters_href => Hash with SnpEff specific parameters {REF}
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment;
    my $conda_prefix_path;
    my $FILEHANDLE;
    my $quiet;
    my $snpeff_parameters_href;
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
        program_parameters_href => {
            default     => {},
            required    => 1,
            store       => \$snpeff_parameters_href,
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
    my $snpeff_version         = $snpeff_parameters_href->{version};
    my @snpeff_genome_versions = @{ $snpeff_parameters_href->{snpeff_genome_versions} };

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => $LOG_NAME,
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    ## Store original working directory
    my $pwd = cwd();

    say {$FILEHANDLE} q{### Install SnpEff};
    $log->info(qq{Writing instructions for SnpEff installation via SHELL});

    my $snpeff_install_path =
      catfile( $conda_prefix_path, qw{ share snpEff }, $snpeff_version );

    # Check if snpEff.jar exists (assumes that snpSift also exists if true)
    if ( -f catfile( $snpeff_install_path, q{snpEff.jar} ) ) {
        $log->info(q{SnpEff is already installed in the specified conda environment.});

        $log->warn(q{This will overwrite the current SnpEff installation.});

        ## Removing jar files but keeps config and data in order to avoid unnecessary downloads
        say {$FILEHANDLE} q{## Removing old SnpEff jar files};
        gnu_rm(
            {
                FILEHANDLE  => $FILEHANDLE,
                force       => 1,
                infile_path => catfile( $snpeff_install_path, q{*.jar} ),
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Remove dead symlinks
        say {$FILEHANDLE} q{## Removing old SnpEff links};
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
            FILEHANDLE   => $FILEHANDLE,
            outfile_path => $snpeff_zip_path,
            ,
            quiet   => $quiet,
            url     => $url,
            verbose => $verbose,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    unzip(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $snpeff_zip_path,
            outdir_path => $temp_dir,
            quiet       => $quiet,
            verbose     => $verbose,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Create dir for SnpEff in conda env if it doesn't exist
    if ( not -d $snpeff_install_path ) {

        say {$FILEHANDLE} q{## Create dir for SnpEff in conda env};
        gnu_mkdir(
            {
                FILEHANDLE       => $FILEHANDLE,
                indirectory_path => $snpeff_install_path,
                parents          => 1,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Make available from conda environment
    say {$FILEHANDLE} q{## Make available from conda environment};
    gnu_mv(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => catfile( $temp_dir, qw{ snpEff *.jar } ),
            outfile_path => $snpeff_install_path,
        }
    );
    print {$FILEHANDLE} $NEWLINE;
    ## Only move config if it isn't already there
    if ( not -f catfile( $snpeff_install_path, q{snpEff.config} ) ) {
        gnu_mv(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => catfile( $temp_dir, qw{ snpEff snpEff.config } ),
                outfile_path => $snpeff_install_path,
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
        my $link_path   = catfile( $conda_prefix_path,   q{bin}, $binary );
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

  GENOME_VERSION:
    foreach my $genome_version (@snpeff_genome_versions) {
        ## Check and if required add the vertebrate mitochondrial codon table to SnpEff config
        check_mt_codon_table(
            {
                config_file    => q{snpEff.config},
                FILEHANDLE     => $FILEHANDLE,
                genome_version => $genome_version,
                quiet          => $quiet,
                share_dir      => $snpeff_install_path,
                verbose        => $verbose,
            }
        );

        my $genome_version_directory_path =
          catdir( $snpeff_install_path, q{data}, $genome_version );

        next GENOME_VERSION if ( -d $genome_version_directory_path );
        ## Write instructions to download SnpEff database.
        ## This is done by install script to avoid race conditin when doing first analysis run in MIP

        say {$FILEHANDLE} q{## Downloading SnpEff database};
        my $jar_path         = catfile( $conda_prefix_path, qw{ bin snpEff.jar} );
        my $config_file_path = catfile( $conda_prefix_path, qw{bin snpEff.config} );
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

    ## Remove the temporary install directory
    say {$FILEHANDLE} q{## Remove temporary install directory};
    gnu_rm(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $temp_dir,
            recursive   => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE x 2;

    return;
}

sub check_mt_codon_table {

## Function : Check and if required add the vertebrate mitochondrial codon table to snpeff config
## Returns  :
## Arguments: $config_file    => The config config_file
##          : $FILEHANDLE     => FILEHANDLE to write to
##          : $genome_version => snpeff genome version
##          : $quiet          => Be quiet
##          : $share_dir      => The conda env shared directory
##          : $verbose        => Set verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $config_file;
    my $FILEHANDLE;
    my $genome_version;
    my $quiet;
    my $share_dir;
    my $verbose;

    my $tmpl = {
        config_file => {
            defined     => 1,
            required    => 1,
            store       => \$config_file,
            strict_type => 1,
        },
        FILEHANDLE => {
            defined  => 1,
            required => 1,
            store    => \$FILEHANDLE,
        },
        genome_version => {
            defined     => 1,
            required    => 1,
            store       => \$genome_version,
            strict_type => 1,
        },
        quiet => {
            allow => [ undef, 0, 1 ],
            store => \$quiet,
        },
        share_dir => {
            defined     => 1,
            required    => 1,
            store       => \$share_dir,
            strict_type => 1,
        },
        verbose => {
            allow => [ undef, 0, 1 ],
            store => \$verbose,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Get logger
    my $log = retrieve_log(
        {
            log_name => $LOG_NAME,
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
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
          run( command => $detect_regexp . catfile( $share_dir, $config_file ) );
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
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $infile_path,
                outfile_path => catfile( $share_dir, $config_file ),
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
