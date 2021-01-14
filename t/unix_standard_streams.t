#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname  };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Unix::Standard_streams} => [qw{ unix_standard_streams }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test unix_standard_streams from Standard_streams.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ < stdinfile.test };

## Can be duplicated with %base and/or %specific to enable testing of each individual argument
my %required_argument = (
    stdinfile_path => {
        input           => q{stdinfile.test},
        expected_output => q{<} . $SPACE . q{stdinfile.test},
    },
);

my %specific_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    stdinfile_path => {
        input           => q{stdinfile.test},
        expected_output => q{<} . $SPACE . q{stdinfile.test},
    },
    stdoutfile_path => {
        input           => q{stdoutfile.test},
        expected_output => q{1>} . $SPACE . q{stdoutfile.test},
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2>} . $SPACE . q{stderrfile.test},
    },
    stderrfile_path_append => {
        input           => q{stderrfile.test},
        expected_output => q{2>>} . $SPACE . q{stderrfile.test},
    },
    stdoutfile_path_append => {
        input           => q{stdoutfile.test},
        expected_output => q{1>>} . $SPACE . q{stdoutfile.test},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&unix_standard_streams;

## Test both arguments
my @arguments = ( \%specific_argument );

foreach my $argument_href (@arguments) {

    my @commands = test_function(
        {
            argument_href              => $argument_href,
            do_test_base_command       => 0,
            function_base_commands_ref => \@function_base_commands,
            module_function_cref       => $module_function_cref,
            required_argument_href     => \%required_argument,
        }
    );
}

done_testing();
