package MIP::Recipes::Pipeline::Install_rd_dna;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use List::Util qw{ any };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## MIPs lib/
use MIP::Constants qw{
  $CLOSE_BRACKET
  $COLON
  $DOT
  $LOG
  $NEWLINE
  $OPEN_BRACKET
  $SINGLE_QUOTE
  $SPACE
};
use MIP::Script::Setup_script qw{ setup_install_script };
use MIP::Set::Parameter qw{ set_programs_for_installation };

## Recipes
use MIP::Recipes::Install::Bedtools qw{ install_bedtools };
use MIP::Recipes::Install::Conda qw{ install_conda_packages };
use MIP::Recipes::Install::Cnvnator qw{ install_cnvnator };
use MIP::Recipes::Install::Expansionhunter qw{ install_expansionhunter };
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
use MIP::Recipes::Install::Svdb qw{ install_svdb };
use MIP::Recipes::Install::Tiddit qw{ install_tiddit };
use MIP::Recipes::Install::Vcf2cytosure qw{ install_vcf2cytosure };
use MIP::Recipes::Install::Vep qw{ install_vep };
use MIP::Recipes::Install::Vt qw{ install_vt };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ pipeline_install_rd_dna };
}

sub pipeline_install_rd_dna {

## Function : Download references recipes for rd_dna pipeline
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this download hash {REF}
##          : $quiet                 => Be quiet
##          : $verbose               => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    ## Default(s)
    my $quiet;
    my $verbose;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$quiet,
            strict_type => 1,
        },
        verbose => {
            default     => $arg_href->{active_parameter_href}{verbose},
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger($LOG);

    ## Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Installation instruction file
    my $file_name_path = catfile( cwd(), q{mip.sh} );

    open $FILEHANDLE, q{>}, $file_name_path
      or $log->logcroak( q{Cannot write to}
          . $SPACE
          . $SINGLE_QUOTE
          . $file_name_path
          . $SINGLE_QUOTE
          . $SPACE
          . $COLON
          . $OS_ERROR
          . $NEWLINE );

    ## Create bash file for writing install instructions
    setup_install_script(
        {
            FILEHANDLE            => $FILEHANDLE,
            file_name             => $file_name_path,
            invoke_login_shell    => $active_parameter_href->{sbatch_mode},
            log                   => $log,
            active_parameter_href => $active_parameter_href,
            remove_dir            => catfile( cwd(), $DOT . q{MIP} ),
            sbatch_mode           => $active_parameter_href->{sbatch_mode},
            set_errexit           => $active_parameter_href->{bash_set_errexit},
            set_nounset           => $active_parameter_href->{bash_set_nounset},
        }
    );

    $log->info( q{Writing install instructions to:} . $SPACE . $file_name_path );

    ## Source conda
    if ( not $active_parameter_href->{sbatch_mode} ) {
        say {$FILEHANDLE} q{## Source conda};
        say {$FILEHANDLE} q{source}
          . $SPACE
          . catfile( $active_parameter_href->{conda_path}, qw{ etc profile.d conda.sh } )
          . $NEWLINE;
    }

    ## Make sure that the cnvnator environment is installed last
    if ( any { $_ eq q{ecnvnator} } @{ $active_parameter_href->{installations} } ) {

        @{ $active_parameter_href->{installations} } =
          grep { !m/ecnvnator/xms } @{ $active_parameter_href->{installations} };
        push @{ $active_parameter_href->{installations} }, q{ecnvnator};
    }

    ## Loop over the selected installations
    foreach my $installation ( @{ $active_parameter_href->{installations} } ) {
        my $env_name = $active_parameter_href->{environment_name}{$installation};

        ## Create some space
        $log->info( $OPEN_BRACKET . $installation . $CLOSE_BRACKET );
        $log->info( q{Working on environment: } . $env_name );

        ## Process input parameters to get a correct combination of programs that are to be installed
        set_programs_for_installation(
            {
                installation          => $installation,
                log                   => $log,
                active_parameter_href => $active_parameter_href,
            }
        );

        ## Installing Conda packages
        install_conda_packages(
            {
                conda_env => $active_parameter_href->{environment_name}{$installation},
                conda_env_path =>
                  $active_parameter_href->{$installation}{conda_prefix_path},
                conda_no_update_dep => $active_parameter_href->{conda_no_update_dep},
                conda_packages_href => $active_parameter_href->{$installation}{conda},
                conda_update        => $active_parameter_href->{conda_update},
                FILEHANDLE          => $FILEHANDLE,
                snpeff_genome_versions_ref =>
                  $active_parameter_href->{$installation}{snpeff_genome_versions},
                quiet   => $active_parameter_href->{quiet},
                verbose => $active_parameter_href->{verbose},
            }
        );

        ## Install PIP packages
        install_pip_packages(
            {
                conda_env  => $active_parameter_href->{environment_name}{$installation},
                FILEHANDLE => $FILEHANDLE,
                pip_packages_href => $active_parameter_href->{$installation}{pip},
                quiet             => $active_parameter_href->{quiet},
            }
        );

        ### Install shell programs
        ## Create dispatch table for shell installation subs
        my %shell_subs = (
            bedtools        => \&install_bedtools,
            cnvnator        => \&install_cnvnator,
            expansionhunter => \&install_expansionhunter,
            mip_scripts     => \&install_mip_scripts,
            picard          => \&install_picard,
            plink2          => \&install_plink2,
            rhocall         => \&install_rhocall,
            sambamba        => \&install_sambamba,
            snpeff          => \&install_snpeff,
            svdb            => \&install_svdb,
            tiddit          => \&install_tiddit,
            vcf2cytosure    => \&install_vcf2cytosure,
            vep             => \&install_vep,
            vt              => \&install_vt,
        );

        ## Launch shell installation subroutines
      SHELL_PROGRAM:
        for my $shell_program ( keys %{ $active_parameter_href->{$installation}{shell} } )
        {

            $shell_subs{$shell_program}->(
                {
                    conda_environment =>
                      $active_parameter_href->{environment_name}{$installation},
                    conda_prefix_path =>
                      $active_parameter_href->{$installation}{conda_prefix_path},
                    FILEHANDLE => $FILEHANDLE,
                    noupdate   => $active_parameter_href->{noupdate},
                    program_parameters_href =>
                      \%{ $active_parameter_href->{$installation}{shell}{$shell_program}
                      },
                    quiet   => $active_parameter_href->{quiet},
                    verbose => $active_parameter_href->{verbose},
                }
            );
        }

        ## Download reference genome if requested, only downloads to MIPs main environment
        if (    ( $active_parameter_href->{reference_dir} )
            and ( $installation eq q{emip} ) )
        {

            download_genome_references(
                {
                    conda_environment => $active_parameter_href->{environment_name}{emip},
                    conda_prefix_path =>
                      $active_parameter_href->{emip}{conda_prefix_path},
                    FILEHANDLE         => $FILEHANDLE,
                    pipeline           => $active_parameter_href->{pipeline},
                    reference_dir_path => $active_parameter_href->{reference_dir},
                    reference_genome_versions_ref =>
                      $active_parameter_href->{reference_genome_versions},
                    quiet   => $active_parameter_href->{quiet},
                    verbose => $active_parameter_href->{verbose},
                }
            );
        }
    }

    ## Write tests
    foreach my $installation ( @{ $active_parameter_href->{installations} } ) {
        ## Get the programs that mip has tried to install
        my @programs_to_test = (
            keys %{ $active_parameter_href->{$installation}{conda} },
            keys %{ $active_parameter_href->{$installation}{pip} },
            keys %{ $active_parameter_href->{$installation}{shell} },
        );

        check_program_installations(
            {
                env_name     => $active_parameter_href->{environment_name}{$installation},
                FILEHANDLE   => $FILEHANDLE,
                installation => $installation,
                programs_ref => \@programs_to_test,
                program_test_command_href =>
                  $active_parameter_href->{program_test_command},
            }
        );
    }

    ## Update/create config
    update_config(
        {
            env_name_href     => $active_parameter_href->{environment_name},
            FILEHANDLE        => $FILEHANDLE,
            installations_ref => $active_parameter_href->{installations},
            log               => $log,
            pipeline          => $active_parameter_href->{pipeline},
            update_config     => $active_parameter_href->{update_config},
            write_config      => $active_parameter_href->{write_config},
        }
    );

    $log->info(q{Finished writing installation instructions for MIP});

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});
    return;
}

1;
