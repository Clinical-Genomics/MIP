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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.03;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Install::Post_installation} =>
          [qw{ build_perl_program_check_command }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Install::Post_installation qw{ build_perl_program_check_command };

diag(   q{Test build_perl_program_check_command from Check.pm v}
      . $MIP::Recipes::Install::Post_installation::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a test program with both test methods
my @programs             = qw{ program };
my %program_test_command = (
    program => {
        path      => q{program},
        execution => q{program --version},
    },
);

## Then return perl oneliner with test code for execution
my @perl_commands = build_perl_program_check_command(
    {
        programs_ref              => \@programs,
        program_test_command_href => \%program_test_command,
    }
);
my $perl_command = join $SPACE, @perl_commands;

my $expected_return =
q{perl -e ' use IPC::Cmd qw{ can_run run }; use Test::More; ok(run(command => q{program --version}, timeout => 60), q{Can execute: program}); done_testing; '};
is( $perl_command, $expected_return, q{Perl oneliner execution test built} );

## Given a program with only a path test
delete $program_test_command{program}{execution};

## Then return perl oneliner with test code for path
@perl_commands = build_perl_program_check_command(
    {
        programs_ref              => \@programs,
        program_test_command_href => \%program_test_command,
    }
);
$perl_command = join $SPACE, @perl_commands;

# Join array to string
$expected_return =
q{perl -e ' use IPC::Cmd qw{ can_run run }; use Test::More; ok(can_run( q{program} ), q{Program in path: program}); done_testing; '};
is( $perl_command, $expected_return, q{Perl oneliner path test built} );

## Given a program that lacks testing methods
%program_test_command = ();

## Then return empty
@perl_commands = build_perl_program_check_command(
    {
        programs_ref              => \@programs,
        program_test_command_href => \%program_test_command,
    }
);
$expected_return = [];

# Test
is_deeply( \@perl_commands, $expected_return, q{Returns on missing methods} );

done_testing;
