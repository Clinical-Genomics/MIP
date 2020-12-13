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


## Constants
Readonly my $MAX_SV_SIZE => 50_000_000;
Readonly my $MIN_SV_SIZE => 15;

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Delly} => [qw{ delly_merge }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Delly qw{ delly_merge };

diag(   q{Test delly_merge from Delly.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ delly merge };

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

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    infile_paths_ref => {
        inputs_ref      => [qw{ s1.bcf s2.bcf s3.bcf }],
        expected_output => q{s1.bcf s2.bcf s3.bcf},
    },
    sv_type => {
        input           => q{DEL},
        expected_output => q{--type DEL},
    },
);

my %specific_argument = (
    max_size => {
        input           => $MAX_SV_SIZE,
        expected_output => q{--maxsize} . $SPACE . $MAX_SV_SIZE,
    },
    min_size => {
        input           => $MIN_SV_SIZE,
        expected_output => q{--minsize} . $SPACE . $MIN_SV_SIZE,
    },
    outfile_path => {
        input           => catfile(qw{ outfile_path_prefix SV_type.txt }),
        expected_output => q{--outfile}
          . $SPACE
          . catfile(qw{ outfile_path_prefix SV_type.txt }),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&delly_merge;

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
