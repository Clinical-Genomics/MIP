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
use MIP::Constants qw{ $COMMA $NEWLINE $SPACE };
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

## Constants
Readonly my $N_PROCESSORS => 4;

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Qc::Peddy} => [qw{ peddy }],
        q{MIP::Test::Fixtures}     => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Qc::Peddy qw{ peddy };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test peddy from Peddy.pm v}
      . $MIP::Program::Qc::Peddy::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ peddy };

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

# Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    case_file_path => {
        input           => catfile(qw{ outcase_file_directory case_id .fam }),
        expected_output => catfile(qw{ outcase_file_directory case_id .fam }),
    },
    infile_path => {
        input           => catfile(qw{ temp_directory infile_prefix .vcf.gz }),
        expected_output => catfile(qw{ temp_directory infile_prefix .vcf.gz }),
    },
    outfile_prefix_path => {
        input           => catfile(qw{ outcase_directory case_id }),
        expected_output => q{--prefix}
          . $SPACE
          . catfile(qw{ outcase_directory case_id }),
    },
);

my %specific_argument = (
    genome_site => {
        input           => q{hg38},
        expected_output => q{--sites} . $SPACE . q{hg38},
    },
    plot => {
        input           => 1,
        expected_output => q{--plot},
    },
    processor_number => {
        input           => $N_PROCESSORS,
        expected_output => q{--procs} . $SPACE . $N_PROCESSORS,
    },
);

# Coderef - enables generalized use of generate call
my $module_function_cref = \&peddy;

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
