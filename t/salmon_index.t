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
    my %perl_module = ( q{MIP::Program::Salmon} => [qw{ salmon_index }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Salmon qw{ salmon_index };

diag(   q{Test salmon_index from Salmon.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

Readonly my $CORES => 12;

my @function_base_commands = qw{ salmon index };

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

my %required_argument = (
    fasta_path => {
        input           => catfile(qw{ dir transcripts.fasta }),
        expected_output => q{--transcripts} . $SPACE . catfile(qw{ dir transcripts.fasta }),
    },
    outfile_path => {
        input           => catfile(qw{ dir index_dir }),
        expected_output => q{--index} . $SPACE . catfile(qw{ dir index_dir }),
    },
);

my %specific_argument = (
    decoy_path => {
        input           => catfile(qw{ path to decoy.txt }),
        expected_output => q{--decoy} . $SPACE . catfile(qw{ path to decoy.txt}),
    },
    fasta_path => {
        input           => catfile(qw{ dir transcripts.fasta }),
        expected_output => q{--transcripts} . $SPACE . catfile(qw{ dir transcripts.fasta }),
    },
    gencode => {
        input           => 1,
        expected_output => q{--gencode},
    },
    outfile_path => {
        input           => catfile(qw{ dir index_dir }),
        expected_output => q{--index} . $SPACE . catfile(qw{ dir index_dir }),
    },
    temp_directory => {
        input           => catdir(qw{ path to temp }),
        expected_output => q{--tmpdir} . $SPACE . catdir(qw{ path to temp }),
    },
    threads => {
        input           => $CORES,
        expected_output => q{--threads} . $SPACE . $CORES,
    },
);

my $module_function_cref = \&salmon_index;

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
