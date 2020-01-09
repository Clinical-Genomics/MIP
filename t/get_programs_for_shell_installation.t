#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Get::Parameter} => [qw{ get_programs_for_shell_installation }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::Parameter qw{ get_programs_for_shell_installation };

diag(   q{Test get_programs_for_shell_installation from Parameter.pm v}
      . $MIP::Get::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

my @returned_shell_programs;
my $expected_shell_programs_ref;

## Given an install set where two programs overlapp betwen bioconda and shell.
## No general shell preference is given,
## but one of the programs have been selected for shell installation.
my %parameter = (
    env_1 => {
        bioconda => {
            test_1 => 123,
            test_2 => 123,
            test_3 => 123,
        },
        shell => {
            test_2 => 234,
            test_3 => 234,
            test_4 => 234,

        },
    },
    prefer_shell  => 0,
    shell_install => [qw{ test_3 }],
);

## When the subroutine is executed
@returned_shell_programs = get_programs_for_shell_installation(
    {
        conda_programs_href        => $parameter{env_1}{bioconda},
        log                        => $log,
        prefer_shell               => $parameter{prefer_shell},
        shell_install_programs_ref => $parameter{shell_install},
        shell_programs_href        => $parameter{env_1}{shell},
    }
);
@returned_shell_programs = sort @returned_shell_programs;

## Then
$expected_shell_programs_ref = [qw{ test_3 test_4 }];
is_deeply( \@returned_shell_programs, $expected_shell_programs_ref,
    q{Testing shell_install} );

## Given a general preference for shell installation
$parameter{prefer_shell}  = 1;
@parameter{shell_install} = [];

## When the subroutine is executed
@returned_shell_programs = get_programs_for_shell_installation(
    {
        conda_programs_href        => $parameter{env_1}{bioconda},
        log                        => $log,
        prefer_shell               => $parameter{prefer_shell},
        shell_install_programs_ref => $parameter{shell_install},
        shell_programs_href        => $parameter{env_1}{shell},
    }
);
@returned_shell_programs = sort @returned_shell_programs;

## Then
$expected_shell_programs_ref = [qw{ test_2 test_3 test_4 }];
is_deeply( \@returned_shell_programs, $expected_shell_programs_ref, q{Prefer shell} );

## Given a general preference for shell installation
@parameter{shell_install} = [qw{ test_3 }];

## When the subroutine is executed
@returned_shell_programs = get_programs_for_shell_installation(
    {
        conda_programs_href        => $parameter{env_1}{bioconda},
        log                        => $log,
        prefer_shell               => $parameter{prefer_shell},
        shell_install_programs_ref => $parameter{shell_install},
        shell_programs_href        => $parameter{env_1}{shell},
    }
);
@returned_shell_programs = sort @returned_shell_programs;

## Then
$expected_shell_programs_ref = [qw{ test_3 }];
is_deeply( \@returned_shell_programs, $expected_shell_programs_ref,
    q{Prefer shell and selected shell program } );

## Given no special requests
$parameter{prefer_shell}  = 0;
@parameter{shell_install} = [];

## When the subroutine is executed
@returned_shell_programs = get_programs_for_shell_installation(
    {
        conda_programs_href        => $parameter{env_1}{bioconda},
        log                        => $log,
        prefer_shell               => $parameter{prefer_shell},
        shell_install_programs_ref => $parameter{shell_install},
        shell_programs_href        => $parameter{env_1}{shell},
    }
);
@returned_shell_programs = sort @returned_shell_programs;

## Then
$expected_shell_programs_ref = [qw{ test_4 }];
is_deeply( \@returned_shell_programs, $expected_shell_programs_ref,
    q{No special requests} );

done_testing();
