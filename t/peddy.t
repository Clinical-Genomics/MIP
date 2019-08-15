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
use MIP::Test::Commands qw{ test_function };
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
Readonly my $COMMA        => q{,};
Readonly my $NEWLINE      => qq{\n};
Readonly my $N_PROCESSORS => 4;
Readonly my $SPACE        => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Variantcalling::Peddy} => [qw{ peddy }],
        q{MIP::Test::Fixtures}                 => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Variantcalling::Peddy qw{ peddy };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test peddy from Peddy.pm v}
      . $MIP::Program::Variantcalling::Peddy::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ python -m peddy };

my %base_argument = (
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
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
    case_file_path => {
        input           => catfile(qw{ outcase_file_directory case_id .fam }),
        expected_output => catfile(qw{ outcase_file_directory case_id .fam }),
    },
);

my %specific_argument = (
    processor_number => {
        input           => $N_PROCESSORS,
        expected_output => q{--procs} . $SPACE . $N_PROCESSORS,
    },
    plot => {
        input           => 1,
        expected_output => q{--plot},
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
            required_argument_href     => \%required_argument,
            module_function_cref       => $module_function_cref,
            function_base_commands_ref => \@function_base_commands,
            do_test_base_command       => 1,
        }
    );
}

done_testing();
