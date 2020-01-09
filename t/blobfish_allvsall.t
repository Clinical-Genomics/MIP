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
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Commands qw{ test_function };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Program::Blobfish} => [qw{ blobfish_allvsall }],
        q{MIP::Test::Fixtures}    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Blobfish qw{ blobfish_allvsall };

diag(   q{Test blobfish_allvsall from Blobfish.pm v}
      . $MIP::Program::Blobfish::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = (q{BlobFish.py --allvsall});

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
    conditions_ref => {
        inputs_ref      => [qw{ affected unaffected }],
        expected_output => q{--conditions affected unaffected},
    },
    indir_paths_ref => {
        inputs_ref => [ catdir(qw{ path to sample1 }), catdir(qw{ path to sample2 }) ],
        expected_output => q{--paths }
          . catdir(qw{ path to sample1 })
          . $SPACE
          . catdir(qw{ path to sample2 }),
    },
    outdir_path => {
        input           => catdir(qw{ out dir }),
        expected_output => q{--dir } . catdir(qw{ out dir }),
    },
    tx2gene_file_path => {
        input           => catfile(qw{ path to file }),
        expected_output => q{--tx } . catfile(qw{ path to file }),
    },
);

my %specific_argument = (
    conditions_ref => {
        inputs_ref      => [qw{ affected unaffected }],
        expected_output => q{--conditions affected unaffected},
    },
    indir_paths_ref => {
        inputs_ref => [ catdir(qw{ path to sample1 }), catdir(qw{ path to sample2 }) ],
        expected_output => q{--paths }
          . catdir(qw{ path to sample1 })
          . $SPACE
          . catdir(qw{ path to sample2 }),
    },
    outdir_path => {
        input           => catdir(qw{ out dir }),
        expected_output => q{--dir } . catdir(qw{ out dir }),
    },
    tx2gene_file_path => {
        input           => catfile(qw{ path to file }),
        expected_output => q{--tx } . catfile(qw{ path to file }),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&blobfish_allvsall;

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
