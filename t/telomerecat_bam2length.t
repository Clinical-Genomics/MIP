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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Telomerecat} => [qw{ telomerecat_bam2length }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Telomerecat qw{ telomerecat_bam2length };

diag(   q{Test telomerecat_bam2length from Telomerecat.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ telomerecat bam2length };

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
        inputs_ref      => [qw{ alignment1.bam alignment2.bam }],
        expected_output => q{alignment1.bam alignment2.bam},
    },
    outfile_path => {
        input           => q{telomere_length.csv},
        expected_output => q{--output telomere_lengths.csv},
    },
);

my %specific_argument = (
    infile_paths_ref => {
        inputs_ref      => [qw{ alignment1.bam alignment2.bam }],
        expected_output => q{alignment1.bam alignment2.bam},
    },
    outfile_path => {
        input           => q{telomere_lengths.csv},
        expected_output => q{--output telomere_lengths.csv},
    },
    processes => {
        input           => 2,
        expected_output => q{-p 2},
    },
    temp_directory => {
        input           => catdir(qw{ temp path }),
        expected_output => q{--temp_dir} . $SPACE . catdir(qw{ temp path }),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&telomerecat_bam2length;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {

    test_function(
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
