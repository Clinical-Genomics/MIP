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
use autodie qw{ :all };
use Modern::Perl qw{ 2018 };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Commands qw{ test_function };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Gnu::Software::Gnu_less} => [qw{ gnu_less }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Gnu::Software::Gnu_less qw{ gnu_less };

diag(   q{Test gnu_less from Gnu_less.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ less };

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    stderrfile_path_append => {
        input           => q{stderrfile.test},
        expected_output => q{2>> stderrfile.test},
    },
);

## Specific arguments
my %specific_argument = (
    infile_path => {
        input           => q{infile.test},
        expected_output => q{infile.test},
    },
    outfile_path => {
        input           => q{outfile.test},
        expected_output => q{> outfile.test},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gnu_less;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

foreach my $argument_href (@arguments) {

    test_function(
        {
            argument_href              => $argument_href,
            function_base_commands_ref => \@function_base_commands,
            module_function_cref       => $module_function_cref,
        }
    );
}

done_testing();
