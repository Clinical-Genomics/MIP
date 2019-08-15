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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA               => q{,};
Readonly my $CPU                 => 16;
Readonly my $MAX_INTRON_DISTANCE => 10_000;
Readonly my $SPACE               => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Variantcalling::Trinity} => [qw{ trinity_genome_guided }],
        q{MIP::Test::Fixtures}                   => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Variantcalling::Trinity qw{ trinity_genome_guided };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test trinity_genome_guided from Trinity.pm v}
      . $MIP::Program::Variantcalling::Trinity::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ Trinity };

my %base_argument = (
    FILEHANDLE => {
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
    infile_path => {
        input           => q{infile_path},
        expected_output => q{--genome_guided_bam infile_path},
    },
);

my %specific_argument = (
    infile_path => {
        input           => q{infile_path},
        expected_output => q{--genome_guided_bam infile_path},
    },
    max_intron_distance => {
        input           => $MAX_INTRON_DISTANCE,
        expected_output => q{--genome_guided_max_intron} . $SPACE . $MAX_INTRON_DISTANCE,
    },
    number_cpu => {
        input           => $CPU,
        expected_output => q{--CPU} . $SPACE . $CPU,
    },
    max_memory => {
        input           => $CPU,
        expected_output => q{--max_memory} . $SPACE . $CPU . q{G},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&trinity_genome_guided;

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
