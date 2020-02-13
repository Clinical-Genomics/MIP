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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
          [qw{ check_program_installations }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Install::Post_installation qw{ check_program_installations };

diag(   q{Test check_program_installations from Post_installation.pm v}
      . $MIP::Recipes::Install::Post_installation::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log => test_log( { no_screen => 0, } );

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open my $filehandle, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## Given no programs to test
my $env_name     = q{an_env_name};
my $installation = q{emip};
my @programs;
my %program_test_command = (
    fastqc => {
        execution => q{fastqc --version},
        path      => q{fastqc},
    }
);

trap {
    check_program_installations(
        {
            env_name                  => $env_name,
            filehandle                => $filehandle,
            installation              => $installation,
            programs_ref              => \@programs,
            program_test_command_href => \%program_test_command,
        }
    )
};

## Then return before testing
like(
    $trap->stderr,
    qr/No\s+tests\s+available\s+for\s+programs/xms,
    q{Return if no programs}
);

## Given programs to test
@programs = qw{ fastqc };

my @is_ok = trap {
    check_program_installations(
        {
            env_name                  => $env_name,
            filehandle                => $filehandle,
            installation              => $installation,
            programs_ref              => \@programs,
            program_test_command_href => \%program_test_command,
        }
    )
};

## Then write tests for programs
ok( $is_ok[0], q{Wrote test for programs} );

## Close the filehandle
close $filehandle;

done_testing();
