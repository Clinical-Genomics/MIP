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
our $VERSION = 1.00;

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
        q{MIP::Program::Mip}   => [qw{ mip_vercollect }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Mip qw{ mip_vercollect };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test mip_vercollect from Mip.pm v}
      . $MIP::Program::Mip::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ mip vercollect };

my %base_argument = (
    stdoutfile_path => {
        input           => q{stdoutfile.test},
        expected_output => q{1> stdoutfile.test},
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    stderrfile_path_append => {
        input           => q{stderrfile.test},
        expected_output => q{2>> stderrfile.test},
    },
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    infile_path => {
        input           => catfile(qw{ outdata_dir case_id case_id_qc_sample_info.yaml }),
        expected_output => q{--infile}
          . $SPACE
          . catfile(qw{ outdata_dir case_id case_id_qc_sample_info.yaml }),
    },
    outfile_path => {
        input           => catfile(qw{ outcase_directory case_id _qc_metrics.yaml }),
        expected_output => q{--outfile}
          . $SPACE
          . catfile(qw{ outcase_directory case_id _qc_metrics.yaml }),
    },
);

my %specific_argument = (
    log_file_path => {
        input           => catfile(qw{ outcase_directory case_id _vercollect.log }),
        expected_output => q{--log_file}
          . $SPACE
          . catfile(qw{ outcase_directory case_id _vercollect.log }),
    },
    infile_path => {
        input           => catfile(qw{ outdata_dir case_id case_id_qc_sample_info.yaml }),
        expected_output => q{--infile}
          . $SPACE
          . catfile(qw{ outdata_dir case_id case_id_qc_sample_info.yaml }),
    },
    outfile_path => {
        input           => catfile(qw{ outcase_directory case_id _qc_metrics.yaml }),
        expected_output => q{--outfile}
          . $SPACE
          . catfile(qw{ outcase_directory case_id _qc_metrics.yaml }),
    },
);

# Coderef - enables generalized use of generate call
my $module_function_cref = \&mip_vercollect;

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
