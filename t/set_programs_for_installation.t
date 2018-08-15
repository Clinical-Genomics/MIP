#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
use FindBin qw{ $Bin };
use Getopt::Long;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2014 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $COMMA   => q{,};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

### User Options
GetOptions(

    # Display help text
    q{h|help} => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },

    # Display version number
    q{v|version} => sub {
        done_testing();
        say {*STDOUT} $NEWLINE
          . basename($PROGRAM_NAME)
          . $SPACE
          . $VERSION
          . $NEWLINE;
        exit;
    },
    q{vb|verbose} => $VERBOSE,
  )
  or (
    done_testing(),
    help(
        {
            USAGE     => $USAGE,
            exit_code => 1,
        }
    )
  );

BEGIN {

    ### Check all internal dependency modules and imports
## Modules with import
    ## Modules with import
    my %perl_module = (
        q{MIP::Log::MIP_log4perl} => [qw{ initiate_logger }],
        q{MIP::Script::Utils}     => [qw{ help }],
    );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

    ## Modules
    my @modules = (q{MIP::Set::Parameter});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Set::Parameter qw{ set_programs_for_installation };
use MIP::Test::Fixtures qw{ test_log };

diag(   q{Test set_programs_for_installation from Set::Parameter.pm v}
      . $MIP::Set::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create temp logger
my $log = test_log();

## Create starting hash
my %parameter = (
    python3_programs => [],
    select_programs  => [],
    shell_install    => [],
    skip_programs    => [],
    test_env         => {
        conda => {
            bio_prog_1   => 1,
            dependency_1 => 1,
            python       => 2.7,
            snpeff       => 1,
            snpsift      => 1,
        },
        pip => {
            py_prog_1 => 1,
            py_prog_2 => 1,
        },
        shell => {
            bio_prog_1 => {
                conda_dependency => {
                    dependency_1 => 1,
                },
                version => 1,
            },
            bio_prog_2 => {
                conda_dependency => {
                    dependency_1 => 1,
                },
                version => 2,
            },
            snpeff => {
                snpeff_genome_versions => [qw{ GRCh37 }],
                version                => 2,
            },
        },
    },
);

## Given a parameter hash with conflicting options
$parameter{select_programs} = [qw{ bio_prog_1 }];
$parameter{skip_programs}   = [qw{ bio_prog_1 bio_prog_2 }];

## When subroutine is executed
trap {
    set_programs_for_installation(
        {
            installation   => q{test_env},
            parameter_href => \%parameter,
            log            => $log,
        }
      )
};

## Then print FATAL log message and exit
like( $trap->stderr, qr/FATAL/xms, q{Fatal log message} );
ok( $trap->exit, q{Exit signal} );

## Given a parameter hash with the select_program option and installation of severeal environements
$parameter{skip_programs} = [];
$parameter{installations} = [qw{ test_env test_env2 }];

## When subroutine is executed
trap {
    set_programs_for_installation(
        {
            installation   => q{test_env},
            parameter_href => \%parameter,
            log            => $log,
        }
      )
};

## Then print FATAL log message and exit
like( $trap->stderr, qr/FATAL/xms, q{Fatal log message} );
ok( $trap->exit, q{Exit signal} );

## Given a parameter hash with a request to skip programs
$parameter{select_programs} = [];
delete $parameter{installations};
$parameter{skip_programs} = [qw{ py_prog_1 bio_prog_1 python }];

## When subroutine is executed
trap {
    set_programs_for_installation(
        {
            installation   => q{test_env},
            parameter_href => \%parameter,
            log            => $log,
        }
      )
};

##Then warn for no python and produce a matching paramter hash
my %installation = (
    python3_programs => [],
    select_programs  => [],
    shell_install    => [],
    skip_programs    => [qw{ py_prog_1 bio_prog_1 python }],
    test_env         => {
        conda => {
            snpeff       => 1,
            snpsift      => 1,
            dependency_1 => 1,
        },
        pip => {
            py_prog_2 => 1,
        },
        shell => {
            bio_prog_2 => {
                conda_dependency => {
                    dependency_1 => 1,
                },
                version => 2,
            },
        },
        snpeff_genome_versions => [qw{ GRCh37 }],
    },
);
like( $trap->stderr, qr/WARN/xms, q{Warn for no python} );
is_deeply( \%parameter, \%installation, q{Solve installation} );

## Given a selective installation
%parameter = (
    installations    => [qw{ test_env }],
    python3_programs => [],
    select_programs  => [qw{ python snpeff bio_prog_2 }],
    shell_install    => [qw{ snpeff }],
    skip_programs    => [],
    test_env         => {
        conda => {
            bio_prog_1   => 1,
            dependency_1 => 1,
            python       => 2.7,
            snpeff       => 1,
            snpsift      => 1,
        },
        pip => {
            py_prog_1 => 1,
            py_prog_2 => 1,
        },
        shell => {
            bio_prog_1 => {
                conda_dependency => {
                    dependency_1 => 1,
                },
                version => 1,
            },
            bio_prog_2 => {
                conda_dependency => {
                    dependency_1 => 1,
                },
                version => 2,
            },
            snpeff => {
                snpeff_genome_versions => [qw{ GRCh37 }],
                version                => 2,
            },
        },
    },
);

## When subroutine is executed
set_programs_for_installation(
    {
        installation   => q{test_env},
        parameter_href => \%parameter,
        log            => $log,
    }
);

## Then solve the installation as such
%installation = (
    python3_programs => [],
    select_programs  => [qw{ python snpeff bio_prog_2 }],
    shell_install    => [qw{ snpeff }],
    skip_programs    => [],
    installations    => [qw{ test_env }],
    test_env         => {
        conda => {
            dependency_1 => 1,
            python       => 2.7,
        },
        pip   => {},
        shell => {
            bio_prog_2 => {
                conda_dependency => {
                    dependency_1 => 1,
                },
                version => 2,
            },
            snpeff => {
                snpeff_genome_versions => [qw{ GRCh37 }],
                version                => 2,
            },
        },
        snpeff_genome_versions => [qw{ GRCh37 }],
    },
);
is_deeply( \%parameter, \%installation, q{Solve installation} );

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

## Function  : Build the USAGE instructions
## Returns   :
## Arguments : $program_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $program_name;

    my $tmpl = {
        program_name => {
            default     => basename($PROGRAM_NAME),
            store       => \$program_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help     Display this help message
    -v/--version  Display version
END_USAGE
}
