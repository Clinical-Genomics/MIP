package MIP::Main::Install;

use Carp;
use charnames qw{ :full :short };
use Cwd;
use Cwd qw{ abs_path };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname basename fileparse };
use File::Spec::Functions qw{ catfile catdir devnull };

#use FindBin qw{ $Bin };
use Getopt::Long;
use IO::Handle;
use List::Util qw{ any };
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

## Recipes
use MIP::Recipes::Install::Bedtools qw{ install_bedtools };
use MIP::Recipes::Install::Cnvnator qw{ install_cnvnator };
use MIP::Recipes::Install::Conda
  qw{ check_conda_installation setup_conda_env install_bioconda_packages };
use MIP::Recipes::Install::Mip_scripts qw{ install_mip_scripts };
use MIP::Recipes::Install::Picard qw{ install_picard };
use MIP::Recipes::Install::Pip qw{ install_pip_packages };
use MIP::Recipes::Install::Plink2 qw{ install_plink2 };
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
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables that can be optionally exported
    our @EXPORT_OK = qw{ mip_install };

}

## Constants
Readonly my $COLON      => q{:};
Readonly my $COMMA      => q{,};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $TAB        => qq{\t};
Readonly my $UNDERSCORE => q{_};

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

    our $VERSION = q{1.2.35};

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

    ## Check the conda installation and get the conda path
    $parameter{conda_prefix_path} = check_conda_installation(
        {
            conda_dir_path    => $parameter{conda_dir_path},
            conda_env         => $parameter{conda_environment},
            disable_env_check => $parameter{disable_env_check},
            quiet             => $parameter{quiet},
            verbose           => $parameter{verbose},
        }
    );

    if ( not $parameter{vep_cache_dir} ) {

        # Cache directory
        $parameter{shell}{vep}{vep_cache_dir} =
          catdir( $parameter{conda_prefix_path},
            q{ensembl-tools-release-} . $parameter{shell}{vep}{version},
            q{cache} );
    }

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

    ## Process input parameters to get a correct combination of programs that are to be installed
    %parameter = get_programs_for_installation(
        {
            log            => $log,
            parameter_href => \%parameter,
        }
    );

    ## Seting up conda environment and installing default packages
    setup_conda_env(
        {
            conda_packages_href => $parameter{conda_packages},
            conda_env           => $parameter{conda_environment},
            conda_env_path      => $parameter{conda_prefix_path},
            FILEHANDLE          => $FILEHANDLE,
            conda_update        => $parameter{conda_update},
            quiet               => $parameter{quiet},
            verbose             => $parameter{verbose},
        }
    );

    ## Installing bioconda packages
    install_bioconda_packages(
        {
            bioconda_packages_href => $parameter{bioconda},
            conda_env              => $parameter{conda_environment},
            conda_env_path         => $parameter{conda_prefix_path},
            FILEHANDLE             => $FILEHANDLE,
            snpeff_genome_versions_ref =>
              $parameter{shell}{snpeff}{snpeff_genome_versions},
            quiet   => $parameter{quiet},
            verbose => $parameter{verbose},
        }
    );

    ## Install PIP packages
    install_pip_packages(
        {
            conda_env         => $parameter{conda_environment},
            FILEHANDLE        => $FILEHANDLE,
            pip_packages_href => $parameter{pip},
            quiet             => $parameter{quiet},
        }
    );

    ### Install shell programs
    ## Create dispatch table for shell installation subs
    my %shell_subs = (
        bedtools     => \&install_bedtools,
        cnvnator     => \&install_cnvnator,
        mip_scripts  => \&install_mip_scripts,
        picard       => \&install_picard,
        plink2       => \&install_plink2,
        rhocall      => \&install_rhocall,
        sambamba     => \&install_sambamba,
        snpeff       => \&install_snpeff,
        svdb         => \&install_svdb,
        tiddit       => \&install_tiddit,
        vep          => \&install_vep,
        vcf2cytosure => \&install_vcf2cytosure,
        vt           => \&install_vt,
    );

    ## Launch shell installation subroutines
  SHELL_PROGRAM:
    for my $shell_program ( @{ $parameter{shell_programs_to_install} } ) {

        $shell_subs{$shell_program}->(
            {
                conda_environment       => $parameter{conda_environment},
                conda_prefix_path       => $parameter{conda_prefix_path},
                FILEHANDLE              => $FILEHANDLE,
                noupdate                => $parameter{noupdate},
                program_parameters_href => $parameter{shell}{$shell_program},
                quiet                   => $parameter{quiet},
                verbose                 => $parameter{verbose},
            }
        );
    }

    ## Download reference genome if requested
    if ( $parameter{reference_dir} ) {

        download_genome_references(
            {
                conda_environment  => $parameter{conda_environment},
                conda_prefix_path  => $parameter{conda_prefix_path},
                FILEHANDLE         => $FILEHANDLE,
                reference_dir_path => $parameter{reference_dir},
                reference_genome_versions_ref =>
                  $parameter{reference_genome_versions},
                quiet   => $parameter{quiet},
                verbose => $parameter{verbose},
            }
        );
    }

    ## Add final message to FILEHANDLE
    display_final_message(
        {
            bioconda_programs_href => $parameter{bioconda},
            conda_env_name         => $parameter{conda_environment},
            conda_programs_href    => $parameter{conda_packages},
            FILEHANDLE             => $FILEHANDLE,
            pip_programs_href      => $parameter{pip},
            shell_programs_ref     => $parameter{shell_programs_to_install},
        }
    );

    $log->info(q{Finished writing installation instructions for MIP});

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    return;
}

#################
###SubRoutines###
#################

sub get_programs_for_installation {

## Function : Procces the lists of programs that has been seleceted for installation and returns them
## Returns  : %{ $parameter_href }
## Arguments: $log            => Log
##          : $parameter_href => The entire parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $log;

    my $tmpl = {
        log => {
            required => 1,
            defined  => 1,
            store    => \$log,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Remove selected programs from installation
    if ( @{ $parameter_href->{skip_program} } ) {
      PROGRAM:
        foreach my $program ( @{ $parameter_href->{skip_program} } ) {
            delete $parameter_href->{shell}{$program};
            delete $parameter_href->{bioconda}{$program};
            delete $parameter_href->{pip}{$program};
        }
    }
    ## Exit if a python 3 env has ben specified for something else than Chanjo or Genmod
    _assure_python_3_compability(
        {
            log                => $log,
            py3_packages_ref   => $parameter_href->{py3_packages},
            python_version     => $parameter_href->{conda_packages}{python},
            select_program_ref => $parameter_href->{select_program}
        }
    );

    ## Programs that are not installed via conda can have dependencies that
    ## needs to be explicetly installed. Also, depedening on how the analysis
    ## recipes have been written, a module can be dependent on more than one
    ## conda program to function.
    if ( any { $_ eq q{chanjo} } @{ $parameter_href->{select_program} } ) {
        push @{ $parameter_href->{select_program} }, q{sambamba};
    }
    if ( any { $_ eq q{cnvnator} } @{ $parameter_href->{select_program} } ) {
        push @{ $parameter_href->{select_program} }, qw{ samtools bcftools };
    }
    if ( any { $_ eq q{peddy} } @{ $parameter_href->{select_program} } ) {
        push @{ $parameter_href->{select_program} }, qw{ bcftools };
    }
    if ( any { $_ eq q{svdb} } @{ $parameter_href->{select_program} } ) {
        push @{ $parameter_href->{select_program} },
          qw{ bcftools htslib picard vcfanno vt };
    }
    if ( any { $_ eq q{vcf2cytosure} } @{ $parameter_href->{select_program} } )
    {
        push @{ $parameter_href->{select_program} }, qw{ libxml2 libxslt };
    }
    if ( any { $_ eq q{vep} } @{ $parameter_href->{select_program} } ) {
        push @{ $parameter_href->{select_program} }, qw{ bcftools htslib };
    }

    ## Remove all programs except those selected for installation
    if ( @{ $parameter_href->{select_program} } ) {
        my @programs = (
            keys %{ $parameter_href->{shell} },
            keys %{ $parameter_href->{bioconda} },
            keys %{ $parameter_href->{pip} },
        );
        my @programs_to_skip =
          array_minus( @programs, @{ $parameter_href->{select_program} } );

      PROGRAM:
        foreach my $program (@programs_to_skip) {
            delete $parameter_href->{shell}{$program};
            delete $parameter_href->{bioconda}{$program};
            delete $parameter_href->{pip}{$program};
        }
    }

    ## Some programs have conflicting dependencies and require seperate environments to function properly
    ## These are excluded from installation unless specified with the select_program flag
    my @conflicting_programs = qw{ cnvnator peddy svdb vep };
  CONFLICTING_PROGRAM:
    foreach my $conflicting_program (@conflicting_programs) {
        if (
            not any { $_ eq $conflicting_program }
            @{ $parameter_href->{select_program} }
          )
        {
            delete $parameter_href->{shell}{$conflicting_program};
            delete $parameter_href->{bioconda}{$conflicting_program};
        }
    }

    ## Assure a gcc version of 4.8 in the case of a cnvnator installation
    if ( any { $_ eq q{cnvnator} } @{ $parameter_href->{select_program} } ) {
        $parameter_href->{conda_packages}{gcc} = q{4.8};
    }

    ## Exclude Chanjo, Genmod and Variant_integrity unless a python 3 env has been specified
    $parameter_href = _assure_python_2_compability(
        {
            log            => $log,
            parameter_href => $parameter_href,
        }
    );

    my @shell_programs_to_install = get_programs_for_shell_installation(
        {
            conda_programs_href        => $parameter_href->{bioconda},
            log                        => $log,
            prefer_shell               => $parameter_href->{prefer_shell},
            shell_install_programs_ref => $parameter_href->{shell_install},
            shell_programs_href        => $parameter_href->{shell},
        }
    );

    ## Removing the bioconda packages that has been selected to be installed via SHELL
    delete @{ $parameter_href->{bioconda} }{@shell_programs_to_install};
    ## Special case for snpsift since it is installed together with SnpEff
    ## if shell installation of SnpEff has been requested.
    if ( any { $_ eq q{snpeff} } @shell_programs_to_install ) {
        delete $parameter_href->{bioconda}{snpsift};
    }
    $parameter_href->{shell_programs_to_install} = [@shell_programs_to_install];

    return %{$parameter_href};
}

sub _assure_python_3_compability {

## Function : Test if specified programs are to be installed in a python 3 environment
## Returns  :
## Arguments: $sub_log            => Log
##          : $py3_packages_ref   => Array with packages that requires python 3 {REF}
##          : $python_version     => The python version that are to be used for the environment
##          : $select_program_ref => Programs selected for installation by the user {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sub_log;
    my $py3_packages_ref;
    my $python_version;
    my $select_program_ref;

    my $tmpl = {
        log => {
            required => 1,
            defined  => 1,
            store    => \$sub_log,
        },
        py3_packages_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$py3_packages_ref,
        },
        python_version => {
            required => 1,
            defined  => 1,
            allow    => qr{
                         ^( [23] )    # Assert that the python major version starts with 2 or 3
                         [.]            # Major version separator
                         ( \d+$        # Assert that the minor version is a digit
                         | \d+ [.] \d+$ ) # Case when minor and patch version has been supplied, allow only digits
                         }xms,
            store => \$python_version,
        },
        select_program_ref => {
            required => 1,
            default  => [],
            store    => \$select_program_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Array::Utils qw{ array_minus };

    ## Check if a python 3 environment has been specified and a python 2 program has been specified for installation
    if (
        $python_version =~ m{
        3 [.] \d+ |    # Python 3 release with minor version eg 3.6
        3 [.] \d+ [.] \d+ # Python 3 release with minor and patch e.g. 3.6.2
        }xms
        and array_minus( @{$select_program_ref}, @{$py3_packages_ref} )
      )
    {
        $sub_log->fatal(
q{A python 3 env has been specified. Please use a python 2 environment for all programs except:}
              . $NEWLINE
              . join $TAB,
            @{$py3_packages_ref}
        );
        exit 1;
    }
    return;
}

sub _assure_python_2_compability {

## Function : Exclude programs that are not compatible with python 2 and exit when a program with python 3 dependency has been selected.
## Returns  : $parameter_href
## Arguments: $log            => Log
##          : $parameter_href => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $parameter_href;

    my $tmpl = {
        log => {
            required => 1,
            defined  => 1,
            store    => \$log,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Array::Utils qw{ intersect };

    if (
        $parameter_href->{conda_packages}{python} =~
        m{ 2 [.] \d+ |    # Python 2 release with minor version eg 2.7.14
           2 [.] \d+ [.] \d+ # Python 3 release with minor and patch e.g. 3.6.2
        }xms
      )
    {
        ## Delete python 3 packages if a python 2 env has been specified
        foreach my $py3_package ( @{ $parameter_href->{py3_packages} } ) {
            delete $parameter_href->{pip}{$py3_package};
        }

        ## Check if a python 3 package has been selected for installation in a python 2 environment
        if (
            intersect(
                @{ $parameter_href->{select_program} },
                @{ $parameter_href->{py3_packages} }
            )
          )
        {
            $log->fatal(
                q{Please specify a python 3 environment for:}
                  . $NEWLINE
                  . join $TAB,
                @{ $parameter_href->{py3_packages} }
            );
            exit 1;
        }
    }

    return $parameter_href;
}

sub get_programs_for_shell_installation {

## Function  : Get the programs that are to be installed via SHELL
## Returns   : @shell_programs
## Arguments : $conda_programs_href        => Hash with conda progrmas {REF}
##           : $log                        => Log
##           : $prefer_shell               => Path to conda environment
##           : $shell_install_programs_ref => Array with programs selected for shell installation {REF}
##           : $shell_programs_href        => Hash with shell programs {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_programs_href;
    my $log;
    my $prefer_shell;
    my $shell_install_programs_ref;
    my $shell_programs_href;

    my $tmpl = {
        conda_programs_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$conda_programs_href,
        },
        log => {
            required => 1,
            defined  => 1,
            store    => \$log,
        },
        prefer_shell => {
            required    => 1,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$prefer_shell
        },
        shell_install_programs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$shell_install_programs_ref,
        },
        shell_programs_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$shell_programs_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Array::Utils qw{ intersect array_minus unique };

    my @shell_programs = keys %{$shell_programs_href};
    my @conda_programs = keys %{$conda_programs_href};

    if ($prefer_shell) {

        # Only get the selected programs otherwise leave the array unaltered
        if ( @{$shell_install_programs_ref} ) {

            # Get the intersect between the two arrays
            @shell_programs =
              intersect( @shell_programs, @{$shell_install_programs_ref} );
        }
    }
    elsif ( @{$shell_install_programs_ref} ) {

        # Assert that the selected program has shell install instructions.
        my @faulty_selects =
          array_minus( @{$shell_install_programs_ref}, @shell_programs );
        if (@faulty_selects) {
            $log->fatal(
                q{No shell installation instructions available for}
                  . $COLON
                  . join $SPACE,
                @faulty_selects
            );
            exit 1;
        }

        # Get elements in @shell_programs that are not part of the conda hash
        my @shell_only_programs =
          array_minus( @shell_programs, @conda_programs );

        # Add the selected program(s) and remove possible duplicates
        @shell_programs =
          unique( @shell_only_programs, @{$shell_install_programs_ref} );
    }
    else {
        # If no shell preferences only add programs lacking conda counterpart
        @shell_programs = array_minus( @shell_programs, @conda_programs );
    }

    return @shell_programs;
}

sub display_final_message {

## Function : Displays a final message to the user at the end of the installation process
## Returns  :
## Arguments: $bioconda_programs_href => Hash with bioconda programs {REF}
##          : $conda_env_name         => Name of conda environment
##          : $conda_programs_href    => Hash with conda programs {REF}
##          : pip_programs_href       => Hash with pip programs {REF}
##          : shell_programs_ref      => Array with shell programs {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bioconda_programs_href;
    my $conda_env_name;
    my $conda_programs_href;
    my $FILEHANDLE;
    my $pip_programs_href;
    my $shell_programs_ref;

    my $tmpl = {
        bioconda_programs_href => {
            default => {},
            store   => \$bioconda_programs_href,
        },
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
        keys %{$bioconda_programs_href},
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
          q{echo -e "\tMIP's installation script has finished\n"};
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
