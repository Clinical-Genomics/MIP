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
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $EQUALS $SPACE };
use MIP::Test::Commands qw{ test_function };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Program::Gnu::Software::Gnu_grep} => [qw{ gnu_grep }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Gnu::Software::Gnu_grep qw{ gnu_grep };

diag(   q{Test gnu_grep from Gnu_grep.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

use MIP::Program::Gnu::Software::Gnu_grep qw(gnu_grep);
use MIP::Test::Commands qw(test_function);

## Base arguments
my @function_base_commands = qw{ grep };

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
    stdoutfile_path => {
        input           => q{stdoutfile.test},
        expected_output => q{1> stdoutfile.test},
    },
);

## Specific arguments
my %specific_argument = (
    count => {
        input           => 1,
        expected_output => q{--count},
    },
    filter_file_path => {
        input           => q{test_file},
        expected_output => q{--file=test_file},
    },
    infile_path => {
        input           => q{infile.test},
        expected_output => q{infile.test},
    },
    invert_match => {
        input           => 1,
        expected_output => q{--invert-match},
    },
    pattern => {
        input           => q{^chr},
        expected_output => q{--regexp} . $EQUALS . q{^chr},
    },
    word_regexp => {
        input           => 1,
        expected_output => q{--word-regexp},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gnu_grep;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

foreach my $argument_href (@arguments) {

    my @commands = test_function(
        {
            argument_href              => $argument_href,
            do_test_base_command       => 1,
            function_base_commands_ref => \@function_base_commands,
            module_function_cref       => $module_function_cref,
        }
    );
}

done_testing();
