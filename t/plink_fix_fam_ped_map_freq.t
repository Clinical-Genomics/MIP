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
our $VERSION = 1.01;

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
        q{MIP::Program::Plink} => [qw{ plink_fix_fam_ped_map_freq }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Plink qw{ plink_fix_fam_ped_map_freq };

diag(   q{Test plink_fix_fam_ped_map_freq from Plink.pm v}
      . $MIP::Program::Plink::VERSION
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
    fam_file_path => {
        input           => catfile(qw{ outcase_file_directory case_id .fam }),
        expected_output => q{--fam}
          . $SPACE
          . catfile(qw{ outcase_file_directory case_id .fam }),
    },
    outfile_prefix => {
        input           => catfile(qw{ temp_directory $case_id _data }),
        expected_output => q{--out}
          . $SPACE
          . catfile(qw{ temp_directory $case_id _data }),
    },
);

my %specific_argument = (
    allow_no_sex => {
        input           => 1,
        expected_output => q{--allow-no-sex},
    },
    binary_fileset_prefix => {
        input           => catfile(qw{ temp_directory $case_id _data }),
        expected_output => q{--bfile}
          . $SPACE
          . catfile(qw{ temp_directory $case_id _data }),
    },
    fam_file_path => {
        input           => catfile(qw{ outcase_file_directory case_id .fam }),
        expected_output => q{--fam}
          . $SPACE
          . catfile(qw{ outcase_file_directory case_id .fam }),
    },
    freqx => {
        input           => 1,
        expected_output => q{--freqx},
    },
    make_just_fam => {
        input           => 1,
        expected_output => q{--make-just-fam},
    },
    outfile_prefix => {
        input           => catfile(qw{ temp_directory $case_id _data }),
        expected_output => q{--out}
          . $SPACE
          . catfile(qw{ temp_directory $case_id _data }),
    },
    recode => {
        input           => 1,
        expected_output => q{--recode},
    },
);

# Coderef - enables generalized use of generate call
my $module_function_cref = \&plink_fix_fam_ped_map_freq;

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
