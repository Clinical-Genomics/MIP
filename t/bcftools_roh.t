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
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Commands qw{ test_function };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Program::Bcftools} => [qw{ bcftools_roh }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Bcftools qw{ bcftools_roh };

diag(   q{Test bcftools_roh from Bcftools}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ bcftools };

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
        input           => q{stdoutfile_path},
        expected_output => q{1> stdoutfile_path},
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    infile_path => {
        input           => q{infile.test},
        expected_output => q{infile.test},
    },
);

my %specific_argument = (
    af_file_path => {
        input           => catfile(qw{ path to af_file }),
        expected_output => q{--AF-file} . $SPACE . catfile(qw{ path to af_file }),
    },
    af_tag => {
        input           => q{GNOMAD_AF},
        expected_output => q{--AF-tag GNOMAD_AF},
    },
    outfile_path => {
        input           => q{outfile.txt},
        expected_output => q{-o outfile.txt},
    },
    samples_ref => {
        inputs_ref      => [qw{ idref1 idref2 idref3 }],
        expected_output => q{--samples idref1,idref2,idref3},
    },
    skip_indels => {
        input           => 1,
        expected_output => q{--skip-indels},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&bcftools_roh;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {
    my @commands = test_function(
        {
            argument_href              => $argument_href,
            required_argument_href     => \%required_argument,
            module_function_cref       => $module_function_cref,
            function_base_commands_ref => \@function_base_commands,
            do_test_base_command       => 1,
        }
    );
}

done_testing();
