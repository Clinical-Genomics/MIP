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
        q{MIP::Program::BootstrapAnn} => [qw{ bootstrapann }],
        q{MIP::Test::Fixtures}        => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::BootstrapAnn qw{ bootstrapann };

diag(   q{Test bootstrapann from BootstrapAnn.pm v}
      . $MIP::Program::BootstrapAnn::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my @function_base_commands = qw{ BootstrapAnn.py };

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

my %required_argument = (
    ase_file_path => {
        input           => catfile(qw{ path to ase }),
        expected_output => q{--ase} . $SPACE . catfile(qw{ path to ase }),
    },
    vcf_infile_path => {
        input           => catfile(qw{ path to input_vcf }),
        expected_output => q{--vcf} . $SPACE . catfile(qw{ path to input_vcf }),
    },
);

my %specific_argument = (
    ase_file_path => {
        input           => catfile(qw{ path to ase }),
        expected_output => q{--ase} . $SPACE . catfile(qw{ path to ase }),
    },
    vcf_infile_path => {
        input           => catfile(qw{ path to input_vcf }),
        expected_output => q{--vcf} . $SPACE . catfile(qw{ path to input_vcf }),
    },
);

my $module_function_cref = \&bootstrapann;

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
