package MIP::Check::Installation;

use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## Cpanm
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_and_add_dependencies
      check_existing_installation
      check_python_compability
    };
}

## Constants
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };
Readonly my $TAB     => qq{\t};

sub check_existing_installation {

## Function : Checks if the program has already been installed and optionally removes the current installation.
##          : Returns "1" if the program is found and a noupdate flag has been provided
## Returns  : $install_check
## Arguments: $conda_environment      => Conda environment
##          : $conda_prefix_path      => Path to conda environment
##          : $FILEHANDLE             => Filehandle to write to
##          : $log                    => Log to write messages to
##          : $noupdate               => Do not update
##          : $program_directory_path => Path to program directory
##          : $program_name                => Program name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment;
    my $conda_prefix_path;
    my $FILEHANDLE;
    my $log;
    my $noupdate;
    my $program_directory_path;
    my $program_name;

    my $tmpl = {
        program_directory_path => {
            required    => 1,
            store       => \$program_directory_path,
            strict_type => 1,
        },
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
            strict_type => 1,
        },
        conda_prefix_path => {
            defined     => 1,
            required    => 1,
            store       => \$conda_prefix_path,
            strict_type => 1,
        },
        conda_environment => {
            store       => \$conda_environment,
            strict_type => 1,
        },
        noupdate => {
            allow       => [ undef, 0, 1 ],
            store       => \$noupdate,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Modules
    use MIP::Gnu::Coreutils qw{ gnu_rm };
    use MIP::Gnu::Findutils qw{ gnu_find };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };

    ## Set default for conda environment if undef
    if ( not $conda_environment ) {
        $conda_environment = q{root/base};
    }

    ## Check if installation directory exists
    if ( -d $program_directory_path ) {
        $log->info( $program_name
              . $SPACE
              . q{is already installed in conda environment: }
              . $conda_environment );

        if ($noupdate) {
            $log->info( q{Skipping writting installation instructions for }
                  . $program_name );
            say {$FILEHANDLE}
              q{## Skipping writting installation instructions for }
              . $program_name;
            say {$FILEHANDLE} $NEWLINE;

            return 1;
        }

        $log->warn(
            qq{This will overwrite the current $program_name installation});

        say {$FILEHANDLE} qq{## Removing old $program_name directory};
        gnu_rm(
            {
                FILEHANDLE  => $FILEHANDLE,
                force       => 1,
                infile_path => $program_directory_path,
                recursive   => 1,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        say {$FILEHANDLE} qq{## Removing old $program_name links};
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

    $log->info(
        qq{Writing instructions for $program_name installation via SHELL});

    return 0;
}

sub check_python_compability {

## Function : Test if specified programs are to be installed in a python 3 environment
## Returns  :
## Arguments: $installation_set_href => The environment specific installation hash {REF}
##          : $log                   => Log
##          : $python3_programs_ref  => Programs requiring python 3
##          : $python_version        => The python version that are to be used for the environment
##          : $select_programs_ref   => Programs selected for installation by the user {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $installation_set_href;
    my $log;
    my $python3_programs_ref;
    my $python_version;
    my $select_programs_ref;

    my $tmpl = {
        installation_set_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$installation_set_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        python3_programs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$python3_programs_ref,
            strict_type => 1,
        },
        python_version => {
            required => 1,
            store    => \$python_version,
        },
        select_programs_ref => {
            default  => [],
            required => 1,
            store    => \$select_programs_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Array::Utils qw{ intersect unique };

    ## Display warning if python isn't part of the installation
    if ( not $python_version ) {
        $log->warn(
q{Python is not part of the installation. Skipping python compability check.}
        );
        return;
    }

    ## Check format of python version
    if (
        $python_version !~ m{
        ^(?: [23] )      # Assert that the python major version starts with 2 or 3
        [.]              # Major version separator
        (?: \d+$         # Assert that the minor version is a digit
        | \d+ [.] \d+$ ) # Case when minor and patch version has been supplied, allow only digits 
        }xms
      )
    {
        $log->fatal( q{Please specify a python 2 or 3 version, given: }
              . $python_version );
        exit 1;
    }

    ## Cover the case where no program has been actively chosen for installation
    my @programs_to_check;
    if ( not defined $select_programs_ref ) {
        @programs_to_check = unique(
            keys %{ $installation_set_href->{conda} },
            keys %{ $installation_set_href->{shell} },
            keys %{ $installation_set_href->{pip} },
        );
    }
    else {
        @programs_to_check = @{$select_programs_ref};
    }

    my @conflicts = intersect( @programs_to_check, @{$python3_programs_ref} );

    ## Check if a python 2 environment has been specified and a python 3
    ## program has been specified for installation in that environment
    if ( ( $python_version =~ m/^2/xms ) and ( scalar @conflicts > 0 ) ) {
        $log->fatal(
            q{Please use a python 3 environment for:} . $NEWLINE . join $TAB,
            @conflicts );
        exit 1;
    }
    return;
}

sub check_and_add_dependencies {

## Function : Check if shell program dependencies are already part of the installation
## Returns  :
## Arguments: $conda_program_href => Hash with conda programs to be installed {REF}
##          : $dependency_href    => Hash with dependencies {REF}
##          : $log                => Log
##          : $shell_program      => Shell program

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_program_href;
    my $dependency_href;
    my $log;
    my $shell_program;

    my $tmpl = {
        conda_program_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$conda_program_href,
            strict_type => 1,
        },
        dependency_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$dependency_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        shell_program => {
            defined     => 1,
            required    => 1,
            store       => \$shell_program,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  DEPENDENCY:
    foreach my $dependency ( keys %{$dependency_href} ) {

        ## Add dependency to conda installation if missing
        if ( not $conda_program_href->{$dependency} ) {
            $conda_program_href->{$dependency} =
              $dependency_href->{$dependency};
            next DEPENDENCY;
        }

        ## Check if version is specified, do nothing if the same version is already part of the installation
        if ( defined $dependency_href->{$dependency} ) {

            ## Exit if the version of the dependency conflicts with what is already part of the conda installation
            if (
                ( defined $conda_program_href->{$dependency} )
                and ( $dependency_href->{$dependency} ne
                    $conda_program_href->{$dependency} )
              )
            {
                $log->fatal(
qq{$shell_program is dependent on $dependency version: $dependency_href->{$dependency}}
                );
                $log->fatal(
qq{The conda installation specifies version: $conda_program_href->{$dependency} of $shell_program}
                );
                exit 1;
            }
        }
    }
    return;
}

1;
