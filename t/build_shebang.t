#! /usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile devnull  };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use List::MoreUtils qw{ any };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

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
        q{MIP::Language::Shell} => [qw{ build_shebang }],
        q{MIP::Test::Fixtures}  => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Language::Shell qw{ build_shebang };
use MIP::Test::Writefile qw{ test_write_to_file };

diag(   q{Test build_shebang from Shell.pm v}
      . $MIP::Language::Shell::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $separator = q{\n};

## Base arguments
my $batch_shebang = q{#!} . $SPACE;

my $bash_bin_path =
  catfile( dirname( dirname( devnull() ) ), qw(usr bin env bash) );

## Specific arguments
my %argument = (
    bash_bin_path => {
        input           => $bash_bin_path,
        expected_output => $batch_shebang . $bash_bin_path . q{ --login},
    },
    invoke_login_shell => {
        input           => 1,
        expected_output => $batch_shebang . $bash_bin_path . q{ --login},
    },
);

my @commands = build_shebang(
    {
        bash_bin_path      => $argument{bash_bin_path}{input},
        invoke_login_shell => $argument{invoke_login_shell}{input},
    }
);

## Testing return of commands
foreach my $key ( keys %argument ) {

    # Unpack expected output
    my $expected_output = $argument{$key}{expected_output};

    ## Then the returned commands should match the expected output
    ok( ( any { $_ eq $expected_output } @commands ), q{Argument: } . $key );
}

## Testing write to file

## Given fake arguments
my @args = (
    bash_bin_path => $bash_bin_path,
    filehandle    => undef,
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&build_shebang;

my @function_base_commands = ( $batch_shebang . $bash_bin_path );

## Then write to file should return without errors
test_write_to_file(
    {
        args_ref             => \@args,
        base_commands_ref    => \@function_base_commands,
        module_function_cref => $module_function_cref,
        separator            => $separator,
    }
);

done_testing();
