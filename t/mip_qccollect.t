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
    my %perl_module = ( q{MIP::Program::Mip} => [qw{ mip_qccollect }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Mip qw{ mip_qccollect };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test mip_qccollect from Mip.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ mip qccollect };

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
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    infile_path => {
        input           => catfile(qw{ outdata_dir case_id case_id_qc_sample_info.yaml }),
        expected_output => q{--sample_info_file}
          . $SPACE
          . catfile(qw{ outdata_dir case_id case_id_qc_sample_info.yaml }),
    },
    outfile_path => {
        input           => catfile(qw{ outcase_directory case_id _qc_metrics.yaml }),
        expected_output => q{--outfile}
          . $SPACE
          . catfile(qw{ outcase_directory case_id _qc_metrics.yaml }),
    },
    regexp_file_path => {
        input           => q{qc_regexp_-v1.13-.yaml},
        expected_output => q{--regexp_file qc_regexp_-v1.13-.yaml},
    },
);

my %specific_argument = (
    eval_metric_file => {
        input           => catfile(q{metrics.yaml}),
        expected_output => q{--eval_metric_file} . $SPACE . catfile(q{metrics.yaml}),
    },
    infile_path => {
        input           => catfile(qw{ outdata_dir case_id case_id_qc_sample_info.yaml }),
        expected_output => q{--sample_info_file}
          . $SPACE
          . catfile(qw{ outdata_dir case_id case_id_qc_sample_info.yaml }),
    },
    log_file_path => {
        input           => catfile(qw{ outcase_directory case_id _qccollect.log }),
        expected_output => q{--log_file}
          . $SPACE
          . catfile(qw{ outcase_directory case_id _qccollect.log }),
    },
    outfile_path => {
        input           => catfile(qw{ outcase_directory case_id _qc_metrics.yaml }),
        expected_output => q{--outfile}
          . $SPACE
          . catfile(qw{ outcase_directory case_id _qc_metrics.yaml }),
    },
    regexp_file_path => {
        input           => q{qc_regexp_-v1.13-.yaml},
        expected_output => q{--regexp_file qc_regexp_-v1.13-.yaml},
    },
    skip_evaluation => {
        input           => 1,
        expected_output => q{--skip_evaluation},
    },
    store_metrics_outfile => {
        input => catfile(qw{ outcase_directory case_id case_id_metrics_deliverables.yaml }),
        expected_output => q{--store_metrics_outfile}
          . $SPACE
          . catfile(qw{ outcase_directory case_id case_id_metrics_deliverables.yaml }),
    },
);

# Coderef - enables generalized use of generate call
my $module_function_cref = \&mip_qccollect;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {
    test_function(
        {
            argument_href              => $argument_href,
            function_base_commands_ref => \@function_base_commands,
            do_test_base_command       => 1,
            module_function_cref       => $module_function_cref,
            required_argument_href     => \%required_argument,
        }
    );
}

done_testing();
