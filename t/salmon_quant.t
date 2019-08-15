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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.03;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Variantcalling::Salmon} => [qw{ salmon_quant }],
        q{MIP::Test::Fixtures}                  => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Variantcalling::Salmon qw{ salmon_quant };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test salmon_quant from Salmon.pm v}
      . $MIP::Program::Variantcalling::Salmon::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $READ_FILES_COMMAND => q{pigz -dc};

my @function_base_commands = qw{ salmon quant };

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

my %required_argument = (
    index_path => {
        input           => catfile(qw{ dir index_dir }),
        expected_output => q{--index} . $SPACE . catfile(qw{ dir index_dir }),
    },
    outdir_path => {
        input           => catfile(qw{ dir output_dir }),
        expected_output => q{--output} . $SPACE . catdir(qw{ dir output_dir }),
    },
    read_1_fastq_paths_ref => {
        inputs_ref =>
          [ catfile(qw{ dir hi_r1.fastq.gz }), catfile(qw{ dir hello_r1.fastq.gz }) ],
        expected_output => q{-1}
          . $SPACE . q{<(}
          . $READ_FILES_COMMAND
          . $SPACE
          . catfile(qw{ dir hi_r1.fastq.gz })
          . $SPACE
          . catfile(qw{ dir hello_r1.fastq.gz })
          . $SPACE . q{)}
    },
);

my %specific_argument = (
    gc_bias => {
        input           => 1,
        expected_output => q{--gcBias},
    },
    index_path => {
        input           => catfile(qw{ dir index_dir }),
        expected_output => q{--index} . $SPACE . catfile(qw{ dir index_dir }),
    },
    outdir_path => {
        input           => catfile(qw{ dir output_dir }),
        expected_output => q{--output} . $SPACE . catdir(qw{ dir output_dir }),
    },
    read_1_fastq_paths_ref => {
        inputs_ref =>
          [ catfile(qw{ dir hi_r1.fastq.gz }), catfile(qw{ dir hello_r1.fastq.gz }) ],
        expected_output => q{-1}
          . $SPACE . q{<(}
          . $READ_FILES_COMMAND
          . $SPACE
          . catfile(qw{ dir hi_r1.fastq.gz })
          . $SPACE
          . catfile(qw{ dir hello_r1.fastq.gz })
          . $SPACE . q{)}
    },
    read_2_fastq_paths_ref => {
        inputs_ref =>
          [ catfile(qw{ dir hi_r2.fastq.gz }), catfile(qw{ dir hello_r2.fastq.gz }) ],
        expected_output => q{-2}
          . $SPACE . q{<(}
          . $READ_FILES_COMMAND
          . $SPACE
          . catfile(qw{ dir hi_r2.fastq.gz })
          . $SPACE
          . catfile(qw{ dir hello_r2.fastq.gz })
          . $SPACE . q{)}
    },
    validate_mappings => {
        input           => 1,
        expected_output => q{--validateMappings},
    },
);

my $module_function_cref = \&salmon_quant;

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
