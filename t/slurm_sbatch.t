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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
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
        q{MIP::Workloadmanager::Slurm} => [qw{ slurm_sbatch }],
        q{MIP::Test::Fixtures}         => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Workloadmanager::Slurm qw{ slurm_sbatch };

diag(   q{Test slurm_sbatch from Slurm.pm v}
      . $MIP::Workloadmanager::Slurm::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ sbatch };

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
        input           => q{outfile.test},
        expected_output => q{1> outfile.test},
    },
);

## Can be duplicated with %base and/or %specific to enable testing of each individual argument
my %required_argument = (
    dependency_type => {
        input           => q{afterok},
        expected_output => q{--dependency=afterok},
    },
    infile_path => {
        input           => q{infile.test},
        expected_output => q{infile.test},
    },
    job_ids_string => {
        input           => q{:123:124},
        expected_output => q{:123:124},
    },
);

## Specific arguments
my %specific_argument = (
    dependency_type => {
        input           => q{afterok},
        expected_output => q{--dependency=afterok:123:124},
    },
    infile_path => {
        input           => q{infile.test},
        expected_output => q{infile.test},
    },
    job_ids_string => {
        input           => q{:123:124},
        expected_output => q{--dependency=afterok:123:124},
    },
    reservation_name => {
        input           => q{development},
        expected_output => q{--reservation=development},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&slurm_sbatch;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {

    my @commands = test_function(
        {
            argument_href              => $argument_href,
            function_base_commands_ref => \@function_base_commands,
            module_function_cref       => $module_function_cref,
            required_argument_href     => \%required_argument,
        }
    );
}

done_testing();
