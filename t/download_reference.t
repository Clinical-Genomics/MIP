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
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Program::Download::Download_reference} => [qw{ download_reference }],
        q{MIP::Test::Fixtures}                        => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Download::Download_reference qw{ download_reference };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test download_reference from Download_reference.pm v}
      . $MIP::Program::Download::Download_reference::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = q{mip download};

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
    config_file_path => {
        input           => catfile(qw{ a test file }),
        expected_output => q{--config_file} . $SPACE . catfile(qw{ a test file }),
    },
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    pipeline => {
        input           => q{rna},
        expected_output => q{rna},
    },
    reference_genome_versions_ref => {
        inputs_ref => [qw{ test_v1 test_v2 }],
        expected_output =>
          q{--reference_genome_versions test_v1 --reference_genome_versions test_v2},
    },
    reference_dir_path => {
        input           => catdir(qw{ a test path }),
        expected_output => q{--reference_dir} . $SPACE . catdir(qw{ a test path }),
    },
);

my %specific_argument = (
    config_file_path => {
        input           => catfile(qw{ a test file }),
        expected_output => q{--config_file} . $SPACE . catfile(qw{ a test file }),
    },
    pipeline => {
        input           => q{rna},
        expected_output => q{rna},
    },
    reference_genome_versions_ref => {
        inputs_ref => [qw{ test_v1 test_v2 }],
        expected_output =>
          q{--reference_genome_versions test_v1 --reference_genome_versions test_v2},
    },
    reference_dir_path => {
        input           => catdir(qw{ a test path }),
        expected_output => q{--reference_dir} . $SPACE . catdir(qw{ a test path }),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&download_reference;

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

######################
####SubRoutines#######
######################

sub build_usage {

## Function  : Build the USAGE instructions
## Returns   :
## Arguments : $program_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $program_name;

    my $tmpl = {
        program_name => {
            default     => basename($PROGRAM_NAME),
            strict_type => 1,
            store       => \$program_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
