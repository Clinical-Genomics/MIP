package MIP::Main::Install;

use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ dirname basename fileparse };
use File::Spec::Functions qw{ catfile catdir devnull };
use Getopt::Long;
use IO::Handle;
use List::Util qw{ any uniq };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use Time::Piece;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Array::Utils qw{ array_minus };
use Readonly;

## MIPs lib/
use MIP::Set::Parameter qw{ set_programs_for_installation };
use MIP::Gnu::Coreutils qw{ gnu_rm };
use MIP::Language::Shell qw{ create_bash_file };
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Package_manager::Conda qw{ conda_activate conda_deactivate };
use MIP::Set::Parameter qw{ set_conda_env_names_and_paths  };

## Recipes
use MIP::Recipes::Install::Bedtools qw{ install_bedtools };
use MIP::Recipes::Install::Blobfish qw{ install_blobfish };
use MIP::Recipes::Install::BootstrapAnn qw{ install_bootstrapann };
use MIP::Recipes::Install::BootstrapAnn qw{ install_bootstrapann };
use MIP::Recipes::Install::Cnvnator qw{ install_cnvnator };
use MIP::Recipes::Install::Conda qw{ check_conda_installation install_conda_packages };
use MIP::Recipes::Install::Expansionhunter qw{ install_expansionhunter };
use MIP::Recipes::Install::Gtf2bed qw{ install_gtf2bed };
use MIP::Recipes::Install::Mip_scripts qw{ install_mip_scripts };
use MIP::Recipes::Install::Picard qw{ install_picard };
use MIP::Recipes::Install::Pip qw{ install_pip_packages };
use MIP::Recipes::Install::Plink2 qw{ install_plink2 };
use MIP::Recipes::Install::Post_installation
  qw{check_program_installations update_config };
use MIP::Recipes::Install::Reference qw{ download_genome_references };
use MIP::Recipes::Install::Rhocall qw{ install_rhocall };
use MIP::Recipes::Install::Sambamba qw{ install_sambamba };
use MIP::Recipes::Install::SnpEff qw{ install_snpeff };
use MIP::Recipes::Install::Star_fusion qw{ install_star_fusion };
use MIP::Recipes::Install::Svdb qw{ install_svdb };
use MIP::Recipes::Install::Tiddit qw{ install_tiddit };
use MIP::Recipes::Install::Vcf2cytosure qw{ install_vcf2cytosure };
use MIP::Recipes::Install::Vep qw{ install_vep };
use MIP::Recipes::Install::Vt qw{ install_vt };

BEGIN {
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = q{1.4.1};

    # Functions and variables that can be optionally exported
    our @EXPORT_OK = qw{ mip_install };

}

## Constants
Readonly my $CLOSED_BRACKET => q{]};
Readonly my $COLON          => q{:};
Readonly my $COMMA          => q{,};
Readonly my $DOT            => q{.};
Readonly my $NEWLINE        => qq{\n};
Readonly my $OPEN_BRACKET   => q{[};
Readonly my $SPACE          => q{ };
Readonly my $TAB            => qq{\t};
Readonly my $UNDERSCORE     => q{_};

sub mip_install {

## Function : Main script for generating MIP installation scripts
## Returns  :
## Arguments: $parameter_href => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments};

    ## Transfer to lexical variables
    my %parameter = %{$parameter_href};

    ## Create default log name
    if ( not $parameter{log_file} ) {

        ## Get local time
        my $date_time       = localtime;
        my $date_time_stamp = $date_time->datetime;

        $parameter{log_file} =
          catfile( q{mip_install} . $UNDERSCORE . $date_time_stamp . $DOT . q{log} );
    }

    ## Initiate logger
    my $log = initiate_logger(
        {
            file_path => $parameter{log_file},
            log_name  => q{mip_install},
        }
    );
    $log->info( q{Writing log messages to} . $COLON . $SPACE . $parameter{log_file} );

    ## Check the conda installation and set the conda path
    check_conda_installation(
        {
            conda_dir_path    => $parameter{conda_dir_path},
            disable_env_check => $parameter{disable_env_check},
            parameter_href    => \%parameter,
            quiet             => $parameter{quiet},
            verbose           => $parameter{verbose},
        }
    );

    ## Set environment names and environment specific conda paths
    set_conda_env_names_and_paths(
        {
            log            => $log,
            parameter_href => \%parameter,
        }
    );

    ##########
    ###MAIN###
    ##########

    ## Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Installation instruction file
    my $file_name_path = catfile( cwd(), q{mip.sh} );

    open $FILEHANDLE, q{>}, $file_name_path
      or $log->logcroak( q{Cannot write to}
          . $SPACE . q{'}
          . $file_name_path . q{'}
          . $SPACE
          . $COLON
          . $OS_ERROR
          . $NEWLINE );

    ## Create bash file for writing install instructions
    create_bash_file(
        {
            FILEHANDLE         => $FILEHANDLE,
            file_name          => $file_name_path,
            invoke_login_shell => $parameter{sbatch_mode},
            log                => $log,
            parameter_href     => \%parameter,
            remove_dir         => catfile( cwd(), $DOT . q{MIP} ),
            sbatch_mode        => $parameter{sbatch_mode},
            set_errexit        => $parameter{bash_set_errexit},
            set_nounset        => $parameter{bash_set_nounset},
        }
    );

    $log->info( q{Writing install instructions to:} . $SPACE . $file_name_path );

    ## Source conda
    if ( not $parameter{sbatch_mode} ) {
        say {$FILEHANDLE} q{## Source conda};
        say {$FILEHANDLE} q{source}
          . $SPACE
          . catfile( $parameter{conda_dir_path}, qw{ etc profile.d conda.sh } )
          . $NEWLINE;
    }

    ## Loop over the selected installations
    foreach my $installation ( @{ $parameter{installations} } ) {
        my $env_name = $parameter{environment_name}{$installation};

        ## If the main MIP installation is to be installed in conda base env
        if ( not $env_name ) {
            $env_name = q{conda base};
        }

        ## Create some space
        $log->info( $OPEN_BRACKET . $installation . $CLOSED_BRACKET );
        $log->info( q{Working on environment: } . $env_name );

        ## Process input parameters to get a correct combination of programs that are to be installed
        set_programs_for_installation(
            {
                installation   => $installation,
                log            => $log,
                parameter_href => \%parameter,
            }
        );

        ## Installing Conda packages
        install_conda_packages(
            {
                conda_env           => $parameter{environment_name}{$installation},
                conda_env_path      => $parameter{$installation}{conda_prefix_path},
                conda_no_update_dep => $parameter{conda_no_update_dep},
                conda_packages_href => $parameter{$installation}{conda},
                conda_update        => $parameter{conda_update},
                FILEHANDLE          => $FILEHANDLE,
                snpeff_genome_versions_ref =>
                  $parameter{$installation}{snpeff_genome_versions},
                quiet   => $parameter{quiet},
                verbose => $parameter{verbose},
            }
        );

        ## Install PIP packages
        install_pip_packages(
            {
                conda_env         => $parameter{environment_name}{$installation},
                FILEHANDLE        => $FILEHANDLE,
                pip_packages_href => $parameter{$installation}{pip},
                quiet             => $parameter{quiet},
            }
        );

        ### Install shell programs
        ## Create dispatch table for shell installation subs
        my %shell_subs = (
            bootstrapann    => \&install_bootstrapann,
            bedtools        => \&install_bedtools,
            blobfish        => \&install_blobfish,
            bootstrapann    => \&install_bootstrapann,
            cnvnator        => \&install_cnvnator,
            expansionhunter => \&install_expansionhunter,
            gtf2bed         => \&install_gtf2bed,
            mip_scripts     => \&install_mip_scripts,
            picard          => \&install_picard,
            plink2          => \&install_plink2,
            rhocall         => \&install_rhocall,
            sambamba        => \&install_sambamba,
            snpeff          => \&install_snpeff,
            star_fusion     => \&install_star_fusion,
            svdb            => \&install_svdb,
            tiddit          => \&install_tiddit,
            vcf2cytosure    => \&install_vcf2cytosure,
            vep             => \&install_vep,
            vt              => \&install_vt,
        );

        ## Launch shell installation subroutines
      SHELL_PROGRAM:
        for my $shell_program ( keys %{ $parameter{$installation}{shell} } ) {

            $shell_subs{$shell_program}->(
                {
                    conda_environment => $parameter{environment_name}{$installation},
                    conda_prefix_path => $parameter{$installation}{conda_prefix_path},
                    FILEHANDLE        => $FILEHANDLE,
                    noupdate          => $parameter{noupdate},
                    program_parameters_href =>
                      $parameter{$installation}{shell}{$shell_program},
                    quiet   => $parameter{quiet},
                    verbose => $parameter{verbose},
                }
            );
        }

        ## Download reference genome if requested, only downloads to MIPs main environment
        if ( ( $parameter{reference_dir} ) and ( $installation eq q{emip} ) ) {

            download_genome_references(
                {
                    conda_environment  => $parameter{environment_name}{emip},
                    conda_prefix_path  => $parameter{emip}{conda_prefix_path},
                    FILEHANDLE         => $FILEHANDLE,
                    pipeline           => $parameter{pipeline},
                    reference_dir_path => $parameter{reference_dir},
                    reference_genome_versions_ref =>
                      $parameter{reference_genome_versions},
                    quiet   => $parameter{quiet},
                    verbose => $parameter{verbose},
                }
            );
        }
    }

    ## Write tests
    foreach my $installation ( @{ $parameter{installations} } ) {
        ## Get the programs that mip has tried to install
        my @programs_to_test = (
            keys %{ $parameter{$installation}{conda} },
            keys %{ $parameter{$installation}{pip} },
            keys %{ $parameter{$installation}{shell} },
        );

        check_program_installations(
            {
                env_name                  => $parameter{environment_name}{$installation},
                FILEHANDLE                => $FILEHANDLE,
                installation              => $installation,
                log                       => $log,
                programs_ref              => \@programs_to_test,
                program_test_command_href => $parameter{program_test_command},
            }
        );
    }

    ## Update/create config
    update_config(
        {
            env_name_href     => $parameter{environment_name},
            FILEHANDLE        => $FILEHANDLE,
            installations_ref => $parameter{installations},
            log               => $log,
            pipeline          => $parameter{pipeline},
            update_config     => $parameter{update_config},
            write_config      => $parameter{write_config},
        }
    );

    $log->info(q{Finished writing installation instructions for MIP});

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});
    return;
}

1;
