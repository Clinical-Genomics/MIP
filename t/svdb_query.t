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
        q{MIP::Program::Variantcalling::Svdb} => [qw{ svdb_query }],
        q{MIP::Test::Fixtures}                => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Variantcalling::Svdb qw{ svdb_query };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test svdb_query from Svdb.pm v}
      . $MIP::Program::Variantcalling::Svdb::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $BND_DISTANCE  => 10_000;
Readonly my $EVENT_OVERLAP => 0.6;

## Base arguments
my @function_base_commands = qw{ svdb --query };

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
    dbfile_path => {
        input           => catfile(qw{ a test databasefile }),
        expected_output => q{--db} . $SPACE . catfile(qw{ a test databasefile }),
    },
    infile_path => {
        input           => catfile(qw{ a test infile }),
        expected_output => q{--query_vcf} . $SPACE . catfile(qw{ a test infile }),
    },
);

my %specific_argument = (
    bnd_distance => {
        input           => $BND_DISTANCE,
        expected_output => q{--bnd_distance} . $SPACE . $BND_DISTANCE,
    },
    dbfile_path => {
        input           => catfile(qw{ a test databasefile }),
        expected_output => q{--db} . $SPACE . catfile(qw{ a test databasefile }),
    },
    in_allele_count_tag => {
        input           => q{OCC},
        expected_output => q{--in_occ} . $SPACE . q{OCC},
    },
    in_frequency_tag => {
        input           => q{FRQ},
        expected_output => q{--in_frq} . $SPACE . q{FRQ},
    },
    infile_path => {
        input           => catfile(qw{ a test infile }),
        expected_output => q{--query_vcf} . $SPACE . catfile(qw{ a test infile }),
    },
    out_allele_count_tag => {
        input           => q{OCC_out},
        expected_output => q{--out_occ} . $SPACE . q{OCC_out},
    },
    out_frequency_tag => {
        input           => q{FRQ_out},
        expected_output => q{--out_frq} . $SPACE . q{FRQ_out},
    },
    overlap => {
        input           => $EVENT_OVERLAP,
        expected_output => q{--overlap} . $SPACE . $EVENT_OVERLAP,
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&svdb_query;

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
