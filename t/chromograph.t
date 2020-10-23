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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.05;

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
        q{MIP::Program::Chromograph} => [qw{ chromograph }],
        q{MIP::Test::Fixtures}       => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Chromograph qw{ chromograph };

diag(   q{Test chromograph from Chromograph.pm v}
      . $MIP::Program::Chromograph::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $STEP => 10_000;

## Base arguments
my @function_base_commands = qw{ chromograph };

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

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    outdir_path => {
        input           => catdir(qw{ path to out_dir }),
        expected_output => q{--outd} . $SPACE . catdir(qw{ path to out_dir }),
    },
);

my %specific_argument = (
    coverage_file_path => {
        input           => catfile(qw{ path to wig }),
        expected_output => q{--coverage} . $SPACE . catfile(qw{ path to wig }),
    },
    euploid => {
        input           => 1,
        expected_output => q{--euploid},
    },
    ideo_file_path => {
        input           => catfile(qw{ a file.bed }),
        expected_output => q{--ideo} . $SPACE . catfile(qw{ a file.bed }),
    },
    outdir_path => {
        input           => catdir(qw{ path to out_dir }),
        expected_output => q{--outd} . $SPACE . catdir(qw{ path to out_dir }),
    },
    step => {
        input           => $STEP,
        expected_output => q{--step} . $SPACE . $STEP,
    },
    upd_regions_file_path => {
        input           => catfile(qw{ path to bed }),
        expected_output => q{--regions} . $SPACE . catfile(qw{ path to bed }),
    },
    upd_sites_file_path => {
        input           => catfile(qw{ path to bed }),
        expected_output => q{--sites} . $SPACE . catfile(qw{ path to bed }),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&chromograph;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
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
