package MIP::Recipes::Install::Star_fusion;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile devnull };
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
use MIP::Check::Installation qw{ check_existing_installation };
use MIP::Constants
  qw{ $ASTERISK $AT $DOLLAR_SIGN $DOUBLE_QUOTE $DOT $EMPTY_STR $LOG $NEWLINE $SINGLE_QUOTE $SPACE $UNDERSCORE };
use MIP::Gnu::Coreutils qw{ gnu_chmod gnu_echo gnu_rm };
use MIP::Gnu::Software::Gnu_make qw{ gnu_make };
use MIP::Gnu::Software::Gnu_sed qw{ gnu_sed };
use MIP::Language::Shell qw{ build_shebang };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Program::Compression::Tar qw{ tar };
use MIP::Program::Download::Wget qw{ wget };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_star_fusion setup_star_fusion };
}

sub install_star_fusion {

## Function : Install STAR-Fusion
## Returns  :
## Arguments: $conda_environment       => Conda environment
##          : $conda_prefix_path       => Conda prefix path
##          : $FILEHANDLE              => Filehandle to write to
##          : $program_parameters_href => Hash with star_fusion specific parameters {REF}
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment;
    my $conda_prefix_path;
    my $FILEHANDLE;
    my $quiet;
    my $star_fusion_parameters_href;
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
            store       => \$star_fusion_parameters_href,
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
    my $star_fusion_version = $star_fusion_parameters_href->{version};

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

    say {$FILEHANDLE} q{### Install Star-Fusion};
    $log->info(qq{Writing instructions for Star-Fusion installation via SHELL});

    ## Check if installation exists and remove directory
    my $star_fusion_dir =
      catdir( $conda_prefix_path, q{share}, q{STAR-Fusion-v} . $star_fusion_version );
    check_existing_installation(
        {
            conda_environment      => $conda_environment,
            conda_prefix_path      => $conda_prefix_path,
            FILEHANDLE             => $FILEHANDLE,
            log                    => $log,
            program_directory_path => $star_fusion_dir,
            program_name           => q{Star-Fusion},
        }
    );

    ## Download
    say {$FILEHANDLE} q{## Download Star-Fusion};
    my $url =
        q{https://github.com/STAR-Fusion/STAR-Fusion/releases/download/STAR-Fusion-v}
      . $star_fusion_version
      . q{/STAR-Fusion-v}
      . $star_fusion_version
      . $DOT
      . q{FULL.tar.gz};
    my $star_fusion_download_path =
      catfile( $conda_prefix_path, q{share},
        q{STAR-Fusion-v} . $star_fusion_version . $DOT . q{FULL.tar.gz} );
    wget(
        {
            FILEHANDLE   => $FILEHANDLE,
            outfile_path => $star_fusion_download_path,
            quiet        => $quiet,
            url          => $url,
            verbose      => $verbose,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    tar(
        {
            extract           => 1,
            file_path         => $star_fusion_download_path,
            FILEHANDLE        => $FILEHANDLE,
            outdirectory_path => catdir( $conda_prefix_path, q{share} ),
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Make
    say {$FILEHANDLE} q{## Recompile};
    gnu_make(
        {
            makefile_dir => $star_fusion_dir,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Make STAR-Fusion and its utils available};
    setup_star_fusion(
        {
            conda_prefix_path    => $conda_prefix_path,
            FILEHANDLE           => $FILEHANDLE,
            star_fusion_dir_path => $star_fusion_dir,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Remove the downloaded tar file
    say {$FILEHANDLE} q{## Remove tar file};
    gnu_rm(
        {
            infile_path => $star_fusion_download_path,
            recursive   => 1,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE x 2;

    return;
}

sub setup_star_fusion {

## Function : Put executables from STAR-Fusion and its utils in condas path
## Returns  :
## Arguments: $conda_prefix_path    => Conda prefix path
##          : $FILEHANDLE           => Filehandle to write to
##          : $star_fusion_dir_path => Path to STAR-Fusion directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_prefix_path;
    my $FILEHANDLE;
    my $star_fusion_dir_path;

    my $tmpl = {
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
        star_fusion_dir_path => {
            required    => 1,
            store       => \$star_fusion_dir_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @shebang = build_shebang(
        {
            bash_bin_path => catfile( dirname( dirname( devnull() ) ), qw{ bin bash } ),
        }
    );

    my @util_executables = qw{
      gtf_file_to_feature_seqs.pl
      make_super_locus.pl
      remove_long_intron_readthru_transcripts.pl
      restrict_genome_to_chr_entries.pl
    };

    my %executable = map {
        $_ => catfile( $star_fusion_dir_path, qw {ctat-genome-lib-builder util }, $_ )
    } @util_executables;

    ## The addition of prep_genome_lib.pl to bin via conda seems to be broken in STAR-Fusion versin 1.7.0-1.
    $executable{q{prep_genome_lib.pl}} =
      catfile( $star_fusion_dir_path, qw{ ctat-genome-lib-builder prep_genome_lib.pl } );

    if ( _parent_sub() eq q{install_star_fusion} ) {

        $executable{q{STAR-Fusion}} = catfile( $star_fusion_dir_path, q{STAR-Fusion} );
    }

    say {$FILEHANDLE} q{## Creating mock executables};
  MOCK_EXECUTABLE:
    foreach my $mock_executable ( keys %executable ) {

        my $file_content =
            $executable{$mock_executable}
          . $SPACE
          . $DOUBLE_QUOTE
          . $DOLLAR_SIGN
          . $AT
          . $DOUBLE_QUOTE;
        my $outfile_path = catfile( $conda_prefix_path, q{bin}, $mock_executable );
        gnu_echo(
            {
                FILEHANDLE     => $FILEHANDLE,
                outfile_path   => $outfile_path,
                string_wrapper => $SINGLE_QUOTE,
                strings_ref    => \@shebang,
            }
        );
        print {$FILEHANDLE} $NEWLINE;
        gnu_echo(
            {
                FILEHANDLE             => $FILEHANDLE,
                no_trailing_newline    => 1,
                stdoutfile_path_append => $outfile_path,
                string_wrapper         => $SINGLE_QUOTE,
                strings_ref            => [$file_content],
            }
        );
        print {$FILEHANDLE} $NEWLINE;
        gnu_chmod(
            {
                file_path  => $outfile_path,
                FILEHANDLE => $FILEHANDLE,
                permission => q{a+x},
            }
        );
        print {$FILEHANDLE} $NEWLINE;
    }
    return;
}

sub _parent_sub {

## Function : Returns the name of the parent subroutine
## Returns  : $parent_sub
## Arguments:

    Readonly my $MINUS_ONE => -1;
    Readonly my $THREE     => 3;

    ## Get full path to parent subroutine
    my $parent_sub = ( caller 2 )[$THREE];

    ## Isolate subroutine
    $parent_sub = ( split /::/xms, $parent_sub )[$MINUS_ONE];

    return $parent_sub;
}
1;

