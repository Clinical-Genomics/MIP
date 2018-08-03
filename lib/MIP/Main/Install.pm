package MIP::Main::Install;

use Carp;
use charnames qw{ :full :short };
use Cwd;
use Cwd qw{ abs_path };
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
use MIP::Gnu::Coreutils qw{ gnu_rm };
use MIP::Language::Shell qw{ create_bash_file };
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Package_manager::Conda
  qw{ conda_source_activate conda_source_deactivate };
use MIP::Set::Parameter qw{ set_conda_env_names_and_paths  };

## Recipes
use MIP::Recipes::Install::BootstrapAnn qw{ install_bootstrapann };
use MIP::Recipes::Install::Bedtools qw{ install_bedtools };
use MIP::Recipes::Install::Blobfish qw{ install_blobfish };
use MIP::Recipes::Install::BootstrapAnn qw{ install_bootstrapann };
use MIP::Recipes::Install::Cnvnator qw{ install_cnvnator };
use MIP::Recipes::Install::Conda
  qw{ check_conda_installation install_conda_packages };
use MIP::Recipes::Install::Expansionhunter qw{ install_expansionhunter };
use MIP::Recipes::Install::Mip_scripts qw{ install_mip_scripts };
use MIP::Recipes::Install::Picard qw{ install_picard };
use MIP::Recipes::Install::Pip qw{ install_pip_packages };
use MIP::Recipes::Install::Plink2 qw{ install_plink2 };
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
    our $VERSION = q{1.2.7};

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

        $parameter{log_file} = catfile(
            q{mip_install} . $UNDERSCORE . $date_time_stamp . $DOT . q{log} );
    }

    ## Initiate logger
    my $log = initiate_logger(
        {
            file_path => $parameter{log_file},
            log_name  => q{mip_install},
        }
    );
    $log->info(
        q{Writing log messages to} . $COLON . $SPACE . $parameter{log_file} );

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
            FILEHANDLE  => $FILEHANDLE,
            file_name   => $file_name_path,
            log         => $log,
            remove_dir  => catfile( cwd(), $DOT . q{MIP} ),
            set_errexit => $parameter{bash_set_errexit},
            set_nounset => $parameter{bash_set_nounset},
        }
    );

    $log->info(
        q{Writing install instructions to:} . $SPACE . $file_name_path );

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
        get_programs_for_installation(
            {
                log            => $log,
                installation   => $installation,
                parameter_href => \%parameter,
            }
        );

        ## Installing Conda packages
        install_conda_packages(
            {
                conda_env      => $parameter{environment_name}{$installation},
                conda_env_path => $parameter{$installation}{conda_prefix_path},
                conda_packages_href => $parameter{$installation}{conda},
                conda_update        => $parameter{conda_update},
                FILEHANDLE          => $FILEHANDLE,
                snpeff_genome_versions_ref =>
                  $parameter{$installation}{shell}{snpeff}
                  {snpeff_genome_versions},
                quiet   => $parameter{quiet},
                verbose => $parameter{verbose},
            }
        );

        ## Install PIP packages
        install_pip_packages(
            {
                conda_env  => $parameter{environment_name}{$installation},
                FILEHANDLE => $FILEHANDLE,
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
        for my $shell_program (
            @{ $parameter{$installation}{shell_programs_to_install} } )
        {

            $shell_subs{$shell_program}->(
                {
                    conda_environment =>
                      $parameter{environment_name}{$installation},
                    conda_prefix_path =>
                      $parameter{$installation}{conda_prefix_path},
                    FILEHANDLE => $FILEHANDLE,
                    noupdate   => $parameter{noupdate},
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
                    reference_dir_path => $parameter{reference_dir},
                    reference_genome_versions_ref =>
                      $parameter{reference_genome_versions},
                    quiet   => $parameter{quiet},
                    verbose => $parameter{verbose},
                }
            );
        }
    }

    foreach my $installation ( @{ $parameter{installations} } ) {
        ## Add final message to FILEHANDLE
        display_final_message(
            {
                conda_env_name => $parameter{environment_name}{$installation},
                conda_programs_href => $parameter{$installation}{conda},
                FILEHANDLE          => $FILEHANDLE,
                pip_programs_href   => $parameter{$installation}{pip},
                shell_programs_ref =>
                  $parameter{$installation}{shell_programs_to_install},
            }
        );
    }

    $log->info(q{Finished writing installation instructions for MIP});

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    return;
}

#################
###SubRoutines###
#################

sub get_programs_for_installation {

## Function : Proccess the lists of programs that has been selected for or omitted from installation
##          : and update the environment packages
## Returns  :
## Arguments: $installation   => Environment to be installed
##          : $log            => Log
##          : $parameter_href => The entire parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $installation;
    my $log;
    my $parameter_href;

    my $tmpl = {
        installation => {
            defined     => 1,
            required    => 1,
            store       => \$installation,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::Parameter qw{ get_programs_for_shell_installation };
    use MIP::Check::Installation qw{ check_python_compability };

    ## Set install modes to loop over
    my @install_modes = qw{ conda pip shell };

    ## Remove selected programs from installation and gather the rest in an array
    my @programs;
  INSTALL_MODE:
    foreach my $install_mode (@install_modes) {
        delete @{ $parameter_href->{$installation}{$install_mode} }
          { @{ $parameter_href->{skip_programs} } };
        push @programs,
          keys %{ $parameter_href->{$installation}{$install_mode} };
    }
    @programs = uniq @programs;

    ## Exit if a python 2 env has ben specified for a python 3 program
    check_python_compability(
        {
            installation_set_href => $parameter_href->{$installation},
            log                   => $log,
            python3_programs_ref  => $parameter_href->{python3_programs},
            python_version => $parameter_href->{$installation}{conda}{python},
            select_programs_ref => $parameter_href->{select_programs},
        }
    );

    ## Programs that are not installed via conda can have dependencies that
    ## needs to be explicetly installed. Also, depedening on how the analysis
    ## recipes have been written, a module can be dependent on more than one
    ## conda program to function.
    my %dependency = (
        blobfish =>
          [qw{ bioconductor-deseq2 bioconductor-tximport r-optparse r-readr }],
        bootstrapann => [qw{ numpy scipy }],
        chanjo       => [qw{ sambamba }],
        cnvnator     => [qw{ bcftools gcc samtools }],
        peddy        => [qw{ bcftools }],
        picard       => [qw{ java-jdk }],
        star_fusion  => [qw{ star }],
        svdb         => [qw{ bcftools cython htslib numpy picard vcfanno vt }],
        tiddit       => [qw{ cmake numpy scikit-learn }],
        vep          => [qw{ bcftools htslib }],
    );

    if ( @{ $parameter_href->{select_programs} } ) {

        ## Check that only one environment has been specified for installation
        if ( scalar @{ $parameter_href->{installations} } > 1 ) {
            $log->fatal(
q{Please select a single installation environment when using the option select_programs.}
            );
            exit 1;
        }

        ## Add pip since it is required in many cases
        push @{ $parameter_href->{select_programs} }, q{pip};

        ## Add neccessary dependencies
      DEPENDENT:
        foreach my $dependent ( keys %dependency ) {
            if (
                any { $_ eq $dependent }
                @{ $parameter_href->{select_programs} }
              )
            {
                push @{ $parameter_href->{select_programs} },
                  @{ $dependency{$dependent} };
            }
        }

        ## Remove all programs except those selected for installation
        my @programs_to_skip =
          array_minus( @programs, @{ $parameter_href->{select_programs} } );

      INSTALL_MODE:
        foreach my $install_mode (@install_modes) {
            delete @{ $parameter_href->{$installation}{$install_mode} }
              {@programs_to_skip};
        }
    }

    ## Some programs have conflicting dependencies and require seperate environments to function properly
    ## These are excluded from installation unless specified with the select_programs flag

    my @shell_programs_to_install = get_programs_for_shell_installation(
        {
            conda_programs_href => $parameter_href->{$installation}{conda},
            log                 => $log,
            prefer_shell        => $parameter_href->{prefer_shell},
            shell_install_programs_ref => $parameter_href->{shell_install},
            shell_programs_href => $parameter_href->{$installation}{shell},
        }
    );

    ## Removing the conda packages that has been selected to be installed via SHELL
    delete @{ $parameter_href->{$installation}{conda} }
      {@shell_programs_to_install};
    ## Special case for snpsift since it is installed together with SnpEff
    ## if shell installation of SnpEff has been requested.
    if ( any { $_ eq q{snpeff} } @shell_programs_to_install ) {
        delete $parameter_href->{$installation}{conda}{snpsift};
    }
    $parameter_href->{$installation}{shell_programs_to_install} =
      [@shell_programs_to_install];

    return;
}

sub display_final_message {

## Function : Displays a final message to the user at the end of the installation process
## Returns  :
## Arguments: $conda_programs_href => Hash with conda programs {REF}
##          : $conda_env_name         => Name of conda environment
##          : $conda_programs_href    => Hash with conda programs {REF}
##          : pip_programs_href       => Hash with pip programs {REF}
##          : shell_programs_ref      => Array with shell programs {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_env_name;
    my $conda_programs_href;
    my $FILEHANDLE;
    my $pip_programs_href;
    my $shell_programs_ref;

    my $tmpl = {
        conda_env_name => {
            store => \$conda_env_name,
        },
        conda_programs_href => {
            default => {},
            store   => \$conda_programs_href,
        },
        pip_programs_href => {
            default => {},
            store   => \$pip_programs_href,
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE,
        },
        shell_programs_ref => {
            default => [],
            store   => \$shell_programs_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};
    use Array::Utils qw{ unique };

    ## Get the programs that mip has tried to install
    my @programs_to_install = unique(
        keys %{$conda_programs_href},
        keys %{$pip_programs_href},
        @{$shell_programs_ref},
    );

    ## Set conda env variable to Root if no conda environment was specified
    if ( not $conda_env_name ) {
        $conda_env_name = q{Root environment};
    }

    say {$FILEHANDLE}
q{echo -e '\n##############################################################\n'};
    if ( any { $_ eq q{cnvnator} } @programs_to_install ) {
        say {$FILEHANDLE}
q{echo -e "\tMIP's installation script has attempted to install CNVnator"};
        say {$FILEHANDLE} q{echo -e "\tin the specified conda environment: }
          . $conda_env_name . q{\n"};
        say {$FILEHANDLE}
          q{echo -e "\tPlease exit the current session before continuing"};
    }
    else {
        say {$FILEHANDLE}
          q{echo -e "\tMIP has attempted to install the following programs"};
        say {$FILEHANDLE} q{echo -e "\tin the specified conda environment: }
          . $conda_env_name . q{\n"};

        foreach my $program_to_install ( sort @programs_to_install ) {
            say {$FILEHANDLE} q{echo -e "\t"} . $program_to_install;
        }
    }
    say {$FILEHANDLE}
q{echo -e '\n##############################################################\n'};
    return;
}

1;
