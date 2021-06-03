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
        q{MIP::Program::Plink} => [qw{ plink_sex_check }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Plink qw{ plink_sex_check };

diag(   q{Test plink_sex_check from MIP::Program::Plink}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ plink2 };

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    binary_fileset_prefix => {
        input           => catfile(qw{ temp_directory $case_id _data }),
        expected_output => q{--bfile}
          . $SPACE
          . catfile(qw{ temp_directory $case_id _data }),
    },
    extract_file => {
        input           => catfile(qw{ temp_directory case_id _data.prune.in }),
        expected_output => q{--extract}
          . $SPACE
          . catfile(qw{ temp_directory case_id _data.prune.in }),
    },
    outfile_prefix => {
        input           => catfile(qw{ temp_directory $case_id _data }),
        expected_output => q{--out}
          . $SPACE
          . catfile(qw{ temp_directory $case_id _data }),
    },
    read_freqfile_path => {
        input           => catfile(qw{ temp_directory case_id _data.frqx }),
        expected_output => q{--read-freq}
          . $SPACE
          . catfile(qw{ temp_directory case_id _data.frqx }),
    },
);

my %specific_argument = (
    binary_fileset_prefix => {
        input           => catfile(qw{ temp_directory $case_id _data }),
        expected_output => q{--bfile}
          . $SPACE
          . catfile(qw{ temp_directory $case_id _data }),
    },
    extract_file => {
        input           => catfile(qw{ temp_directory case_id _data.prune.in }),
        expected_output => q{--extract}
          . $SPACE
          . catfile(qw{ temp_directory case_id _data.prune.in }),
    },
    outfile_prefix => {
        input           => catfile(qw{ temp_directory $case_id _data }),
        expected_output => q{--out}
          . $SPACE
          . catfile(qw{ temp_directory $case_id _data }),
    },
    read_freqfile_path => {
        input           => catfile(qw{ temp_directory case_id _data.frqx }),
        expected_output => q{--read-freq}
          . $SPACE
          . catfile(qw{ temp_directory case_id _data.frqx }),
    },
    sex_check_min_f => {
        input           => q{0.2 0.75},
        expected_output => q{0.2 0.75},
    },
);

# Coderef - enables generalized use of generate call
my $module_function_cref = \&plink_sex_check;

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
