package MIP::Recipes::Install::Mip_scripts;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ fileparse };
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
use MIP::Constants qw{ $DOT $LOG_NAME $NEWLINE $SPACE $UNDERSCORE };
use MIP::Gnu::Coreutils qw{ gnu_chmod gnu_cp gnu_ln gnu_mkdir};
use MIP::Log::MIP_log4perl qw{ retrieve_log };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.13;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_mip_scripts };
}

sub install_mip_scripts {

## Function : Install mip_scripts
## Returns  :
##          : $conda_environment       => Conda environment
##          : $conda_prefix_path       => Conda prefix path
##          : $filehandle              => Filehandle to write to
##          : $program_parameters_href => Hash with mip_scripts specific parameters {REF}
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment;
    my $conda_prefix_path;
    my $filehandle;
    my $mip_scripts_parameters_href;
    my $quiet;
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
        filehandle => {
            defined  => 1,
            required => 1,
            store    => \$filehandle,
        },
        program_parameters_href => {
            default     => {},
            required    => 1,
            store       => \$mip_scripts_parameters_href,
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

    use MIP::Check::Installation qw{ check_mip_executable };

    ## Unpack parameters
    my $mip_scripts_version = $mip_scripts_parameters_href->{version};

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

    ## Define MIP scripts and yaml files
    my @mip_scripts = qw{ mip };

    my %mip_sub_script = (
        utility_scripts => [qw{ calculate_af.pl max_af.pl }],
        t =>
          [qw{ mip_install.test mip_analyse_rd_dna.test mip_core.t mip_analysis.test }],
        templates => [
            qw{ 643594-miptest_pedigree.yaml
              gene_panel.bed
              grch38_mip_rd_dna_config.yaml
              mip_download_rd_dna_config_-1.0-.yaml
              mip_download_rd_rna_config_-1.0-.yaml
              mip_dragen_rd_dna_config.yaml
              mip_install_rd_dna_config_-1.0-.yaml
              mip_install_rd_rna_config_-1.0-.yaml
              mip_log.yaml
              mip_rd_dna_config.yaml
              mip_rd_dna_vcf_rerun_config.yaml
              mip_rd_rna_config.yaml
              program_test_cmds.yaml
              qc_regexp_-v1.24-.yaml
              rank_model_cmms_-v1.28-.ini
              svrank_model_cmms_-v1.8-.ini
              }
        ],
    );

    my @mip_directories = qw{ lib t definitions };

    say {$filehandle} q{### Install MIP};
    $log->info(q{Writing installation instructions for MIP});

    ## Check if mip installation exists and is executable
    # mip is proxy for all mip scripts
    check_mip_executable(
        {
            conda_prefix_path => $conda_prefix_path,
            log               => $log,
        }
    );

    ## Create directories
    say {$filehandle} q{## Create directories};
  DIRECTORY:
    foreach my $directory ( keys %mip_sub_script ) {

        my $indirectory_path = catdir( $conda_prefix_path, q{bin}, $directory );
        gnu_mkdir(
            {
                filehandle       => $filehandle,
                indirectory_path => $indirectory_path,
                parents          => 1,
            }
        );
        print {$filehandle} $NEWLINE;
    }
    print {$filehandle} $NEWLINE;

    ## Copy directory to conda env
    say {$filehandle} q{## Copy directory to conda env};
  DIRECTORY:
    foreach my $directory (@mip_directories) {

        gnu_cp(
            {
                filehandle   => $filehandle,
                force        => 1,
                infile_path  => catdir( $Bin, $directory ),
                outfile_path => catdir( $conda_prefix_path, q{bin} ),
                recursive    => 1,
            }
        );
        print {$filehandle} $NEWLINE;
    }
    print {$filehandle} $NEWLINE;

    ## Copy mip scripts and sub scripts to conda env and make executable
    say {$filehandle}
      q{## Copy mip scripts and subdirectory scripts to conda env and make executable};

  SCRIPT:
    foreach my $script (@mip_scripts) {

        my $script_no_ending = fileparse( $script, qr/\.[^.]*/xms );
        gnu_cp(
            {
                filehandle   => $filehandle,
                infile_path  => catfile( $Bin, $script ),
                outfile_path => catdir( $conda_prefix_path, q{bin}, $script_no_ending ),
            }
        );
        print {$filehandle} $NEWLINE;

        my $file_path = catfile( $conda_prefix_path, q{bin}, $script_no_ending );
        gnu_chmod(
            {
                file_path  => $file_path,
                filehandle => $filehandle,
                permission => q{a+x},
            }
        );
        say {$filehandle} $NEWLINE;
    }

  DIRECTORY:
    foreach my $directory ( keys %mip_sub_script ) {

      SCRIPT:
        foreach my $script ( @{ $mip_sub_script{$directory} } ) {

            gnu_cp(
                {
                    filehandle   => $filehandle,
                    infile_path  => catfile( $Bin, $directory, $script ),
                    outfile_path => catdir( $conda_prefix_path, q{bin}, $directory ),
                }
            );
            print {$filehandle} $NEWLINE;

            my $file_path = catfile( $conda_prefix_path, q{bin}, $directory, $script );
            gnu_chmod(
                {
                    filehandle => $filehandle,
                    file_path  => $file_path,
                    permission => q{a+x},
                }
            );
            say {$filehandle} $NEWLINE;
        }
    }
    print {$filehandle} $NEWLINE;

    return;
}

1;
