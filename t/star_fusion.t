#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname  };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Program::Variantcalling::Star_fusion} => [qw{ star_fusion }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Variantcalling::Star_fusion qw{ star_fusion };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test star_fusion from Star_fusion.pm v}
      . $MIP::Program::Variantcalling::Star_fusion::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $READ_LENGTH   => 150;
Readonly my $THREAD_NUMBER => 16;

my @function_base_commands = qw{ STAR-Fusion };

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
    genome_lib_dir_path => {
        input           => catfile(qw{ dir genome_lib_dir_path }),
        expected_output => q{--genome_lib_dir}
          . $SPACE
          . catfile(qw{ dir genome_lib_dir_path }),
    },
    sjdb_path => {
        input           => catfile(qw{ dir junctions.tab }),
        expected_output => q{-J} . $SPACE . catfile(qw{ dir junctions.tab }),
    },
    output_directory_path => {
        input           => catfile(qw{ dir }),
        expected_output => q{--output_dir} . $SPACE . catfile(qw{ dir }),
    },
);

my %specific_argument = (
    genome_lib_dir_path => {
        input           => catfile(qw{ dir genome_lib_dir_path }),
        expected_output => q{--genome_lib_dir}
          . $SPACE
          . catfile(qw{ dir genome_lib_dir_path }),
    },
    samples_file_path => {
        input           => catfile(qw{ dir samples_file.txt }),
        expected_output => q{--samples_file}
          . $SPACE
          . catfile(qw{ dir samples_file.txt }),
    },
    sjdb_path => {
        input           => catfile(qw{ dir junctions.tab }),
        expected_output => q{-J} . $SPACE . catfile(qw{ dir junctions.tab }),
    },
    output_directory_path => {
        input           => catfile(qw{ dir }),
        expected_output => q{--output_dir} . $SPACE . catfile(qw{ dir }),
    },
);

my $module_function_cref = \&star_fusion;

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
