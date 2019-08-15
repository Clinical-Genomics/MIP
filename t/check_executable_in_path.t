#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Check::Path}    => [qw{ check_executable_in_path }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Path qw{ check_executable_in_path };

diag(   q{Test check_executable_in_path from Path.pm v}
      . $MIP::Check::Path::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

my %active_parameter = (
    conda_path => catfile( $Bin, qw{ data modules miniconda } ),
    samtools   => 0,
);
my %parameter;

## Given switched off active recipe parameter, when no parameter defined
my $return = check_executable_in_path(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        parameter_href        => \%parameter,
    }
);

## Then return undef
is( $return, undef, q{Skip check when no parameter defined in parameter hash} );

## Given switched off active parameter and defined recipe parameter, when no defined program_executables
%parameter = ( samtools => { type => q{recipe}, }, );

$return = check_executable_in_path(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        parameter_href        => \%parameter,
    }
);

## Then return undef
is( $return, undef, q{Skip check when no program_executables defined in parameter hash} );

## Given switched on active parameter, defined program and program_executables parameter, when program is in path and executable
%active_parameter = (
    conda_path => catfile( $Bin, qw{ data modules miniconda } ),
    samtools   => 1,
);

%parameter = (
    samtools => {
        program_executables => [qw{ samtools }],
        type                => q{recipe},
    },
);

trap {
    check_executable_in_path(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
        }
    )
};

## Then INFO message should broadcast
like( $trap->stderr, qr/INFO/xms, q{Found bin and executable: Throw INFO log message} );

## Given switched on active parameter, defined recipe and program_executables parameter, when program is in path and executable
%active_parameter = (
    conda_path => catfile( $Bin, qw{ data modules miniconda } ),
    samtools   => 1,
);

%parameter = (
    samtools => {
        program_executables => [qw{ no_binary }],
        type                => q{recipe},
    },
);

trap {
    check_executable_in_path(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
        }
    )
};

## Then FATAL message should broadcast
like( $trap->stderr, qr/FATAL/xms, q{No bin and executable - Throw FATAL log message} );

## Given switched on active parameter, defined recipe and program_executables
## parameter, when program is in env path and executable for
## module source environment command
%active_parameter = (
    conda_path => catfile( $Bin, qw{ data modules miniconda } ),
    load_env   => {
        test_env_1 => {
            rankvariant => undef,
            method      => q{conda},
        },
    },
    rankvariant => 1,
);

%parameter = (
    rankvariant => {
        program_executables => [qw{ genmod }],
        type                => q{recipe},
    },
);

trap {
    check_executable_in_path(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
        }
    )
};

## Then INFO message should broadcast
like( $trap->stderr, qr/INFO/xms,
    q{Found bin and executable load env cmd: Throw INFO log message} );

## Given switched on active parameter, defined recipe and program_executables
## parameter, when program is in env path and not executable for
## load env command
%active_parameter = (
    conda_path => catfile( $Bin, qw{ data modules miniconda } ),
    load_env   => {
        test_env_1 => {
            rankvariant => undef,
            method      => q{conda},
        },
    },
    rankvariant => 1,
);

%parameter = (
    rankvariant => {
        program_executables => [qw{ no_binary }],
        type                => q{recipe},
    },
);

trap {
    check_executable_in_path(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
        }
    )
};

## Then FATAL message should broadcast
like( $trap->stderr, qr/FATAL/xms,
    q{No bin and executable in load env cmd - Throw FATAL log message} );

done_testing();
