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
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
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
        q{MIP::Dependency_tree} => [qw{ get_dependency_tree }],
        q{MIP::Test::Fixtures}  => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Dependency_tree qw{ get_dependency_tree_chain };

diag(   q{Test get_dependency_tree_chain from Dependency_tree.pm v}
      . $MIP::Dependency_tree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %parameter;
my %dependency_tree = (
    CHAIN_ALL => [
        q{program_0},
        q{program_1},
        {
            CHAIN_0 => {
                PARALLEL => [
                    q{parallel_program_0},
                    q{parallel_program_1},
                    {
                        PARALLEL_PROGRAM_2 =>
                          [ q{parallel_program_2}, q{parallel_program_3}, ],
                    },
                ],
            },
        },
        q{program_2},
        {
            CHAIN_1 => {
                PARALLEL => [
                    q{parallel_program_4},
                    q{parallel_program_5},
                    {
                        CHAIN_MAIN => [ q{parallel_program_6}, q{parallel_program_7}, ],
                    },
                ],
            },
        },
        q{program_3},
    ],
);

get_dependency_tree_chain(
    {
        dependency_tree_href => \%dependency_tree,
        parameter_href       => \%parameter,
    }
);

my %expected_chain = (
    program_0          => { chain => q{ALL}, },
    program_1          => { chain => q{ALL}, },
    program_2          => { chain => q{ALL}, },
    program_3          => { chain => q{ALL}, },
    parallel_program_0 => { chain => q{PARALLEL_PROGRAM_0}, },
    parallel_program_1 => { chain => q{PARALLEL_PROGRAM_1}, },
    parallel_program_2 => { chain => q{PARALLEL_PROGRAM_2}, },
    parallel_program_3 => { chain => q{PARALLEL_PROGRAM_2}, },
    parallel_program_4 => { chain => q{PARALLEL_PROGRAM_4}, },
    parallel_program_7 => { chain => q{MAIN}, },
);

foreach my $program ( keys %expected_chain ) {

    is(
        $parameter{$program}{chain},
        $expected_chain{$program}{chain},
        q{Chain - } . $program
    );
}
done_testing();
