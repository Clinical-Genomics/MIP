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
    my %perl_module = (
        q{MIP::Program::Glnexus} => [qw{ glnexus_merge }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Glnexus qw{ glnexus_merge };

diag(   q{Test glnexus_merge from Glnexus.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $MEMORY  => 10;
Readonly my $THREADS => 10;

## Base arguments

my @function_base_commands = qw{ glnexus_cli };

my %base_argument = (
    config => {
        input           => q{DeepVariant},
        expected_output => q{--config } . q{DeepVariant},
    },
    infile_paths_ref => {
        inputs_ref      => [ catfile(qw{ dir infile1.vcf }) ],
        expected_output => catfile(qw{ dir infile1.vcf }),
    },
    memory => {
        input           => $MEMORY,
        expected_output => q{--mem-gbytes } . $MEMORY,
    },
    threads => {
        input           => $THREADS,
        expected_output => q{--threads } . $THREADS,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    config => {
        input           => q{DeepVariant},
        expected_output => q{--config } . q{DeepVariant},
    },
    dir => {
        input           => catfile(qw{ dir glnexus }),
        expected_output => q{--dir } . catfile(qw{ dir glnexus }),
    },
    infile_paths_ref => {
        inputs_ref      => [ catfile(qw{ dir infile1.vcf }) ],
        expected_output => catfile(qw{ dir infile1.vcf }),
    },
    threads => {
        input           => $THREADS,
        expected_output => q{--threads } . $THREADS,
    },
);

my %specific_argument = (
    config => {
        input           => q{DeepVariant},
        expected_output => q{--config } . q{DeepVariant},
    },
    dir => {
        input           => catfile(qw{ dir glnexus }),
        expected_output => q{--dir } . catfile(qw{ dir glnexus }),
    },
    infile_paths_ref => {
        inputs_ref      => [ catfile(qw{ dir infile1.vcf }) ],
        expected_output => catfile(qw{ dir infile1.vcf }),
    },
    memory => {
        input           => $MEMORY,
        expected_output => q{--mem-gbytes } . $MEMORY,
    },
    threads => {
        input           => $THREADS,
        expected_output => q{--threads } . $THREADS,
    },
);
## Coderef - enables generalized use of generate call
my $module_function_cref = \&glnexus_merge;

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
