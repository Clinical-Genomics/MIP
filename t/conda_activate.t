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
use autodie qw{ :all };
use Modern::Perl qw{ 2018 };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SEMICOLON $SPACE };
use MIP::Test::Commands qw{ test_function };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Program::Conda} => [qw{ conda_activate }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Conda qw{ conda_activate };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test conda_activate from Conda.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @source_commands = ( q{source}, catfile(qw{path to conda etc profile.d conda.sh}), $SEMICOLON );
my @function_base_commands = ( @source_commands, qw{ conda activate } );

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base and/or %specific to enable testing of each individual argument
my %required_argument = (
    conda_init_path => {
        input           => catfile(qw{path to conda etc profile.d conda.sh}),
        expected_output => q{source}
          . $SPACE
          . catfile(qw{path to conda etc profile.d conda.sh})
          . $SPACE
          . $SEMICOLON,
    },
    env_name => {
        input           => q{test_env},
        expected_output => q{test_env},
    },
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&conda_activate;

## Test both base and function specific arguments
my @arguments = ( \%base_argument );

foreach my $argument_href (@arguments) {
    my @commands = test_function(
        {
            argument_href              => $argument_href,
            do_test_base_command       => 1,
            function_base_commands_ref => \@function_base_commands,
            module_function_cref       => $module_function_cref,
            required_argument_href     => \%required_argument,
        }
    );
}

done_testing();
