#!/usr/bin/env perl

use 5.018;
use Array::Utils qw{ array_diff };
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use Getopt::Long;
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.1';

## Constants
Readonly my $COMMA   => q{,};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

### User Options
GetOptions(

    # Display help text
    q{h|help} => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },

    # Display version number
    q{v|version} => sub {
        done_testing();
        say {*STDOUT} $NEWLINE
          . basename($PROGRAM_NAME)
          . $SPACE
          . $VERSION
          . $NEWLINE;
        exit;
    },
    q{vb|verbose} => $VERBOSE,
  )
  or (
    done_testing(),
    help(
        {
            USAGE     => $USAGE,
            exit_code => 1,
        }
    )
  );

BEGIN {

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Script::Utils} => [qw{ help }], );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Get::Analysis});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Get::Analysis qw{ get_dependency_tree };

diag(   q{Test get_dependency_tree from Analysis.pm v}
      . $MIP::Get::Analysis::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given an initiation map
my %dependency_tree = (
    ALL => [
        q{program_12},
        {
            CHAIN_MAIN => [ qw{ program_0 }, ],
        },
        {
            CHAIN_0 => [ qw{ program_1 }, ],
        },
        q{program_2},
        {
            CHAIN_MAIN => [
                qw{ program_3 program_4 },
                {
                    CHAIN_1 => [qw{ program_5 }],
                },
                {
                    CHAIN_2 => [
                        {
                            PARALLEL => [
                                qw{ parallel_program_0 parallel_program_1 },
                                {
                                    parallel_program_2 => [
                                        qw{ parallel_program_2 parallel_program_3 },
                                    ],
                                },
                            ],
                        },
                        qw{ program_6 program_7 },
                    ],
                },
            ],
        },
        {
            CHAIN_MAIN => [
                {
                    PARALLEL => [
                        qw{ parallel_program_4 parallel_program_5 },
                        {
                            parallel_program_6 =>
                              [ qw{ parallel_program_6 parallel_program_7 }, ],
                        },
                    ],
                },
                qw{ program_8 program_9 },
                {
                    CHAIN_3 => [ qw{ program_10 }, ],
                },
            ],
        },
        q{program_11},
    ],
);

## Expected results from traversing the dependency tree with different starting points
my @expected_programs_12 = qw{ program_12 program_2 program_11 };
my @expected_programs_0 =
  qw{ program_0 program_1 program_2 program_3 program_4 program_5 parallel_program_0 parallel_program_1 parallel_program_2 parallel_program_3 program_6 program_7 parallel_program_4 parallel_program_5 parallel_program_6 parallel_program_7 program_8 program_9 program_10 program_11 };
my @expected_programs_1 = qw{ program_1 program_2 program_11 };
my @expected_parallel_programs_0 =
  qw{ parallel_program_0 program_6 program_7 program_11 };
my @expected_parallel_programs_2 =
  qw{ parallel_program_2 parallel_program_3 program_6 program_7 program_11 };
my @expected_parallel_programs_3 =
  qw{ parallel_program_3 program_6 program_7 program_11 };
my @expected_program_2 = qw{ program_2 program_11 };
my @expected_parallel_programs_5 =
  qw{ parallel_program_5 program_8 program_9 program_10 program_11 };
my @expected_program_11 = qw{ program_11 };

## Define tests
my %test = (
    program_12         => \@expected_programs_12,
    program_0          => \@expected_programs_0,
    program_1          => \@expected_programs_1,
    parallel_program_0 => \@expected_parallel_programs_0,
    parallel_program_2 => \@expected_parallel_programs_2,
    parallel_program_3 => \@expected_parallel_programs_3,
    program_2          => \@expected_program_2,
    parallel_program_5 => \@expected_parallel_programs_5,
    program_11         => \@expected_program_11,
);

## Run tests
while ( my ( $test_key, $expected_programs_ref ) = each %test ) {

    my @start_with_programs;
    my $is_program_found   = 0;
    my $is_chain_found     = 0;
    my $start_with_program = $test_key;

    get_dependency_tree(
        {
            dependency_tree_href    => \%dependency_tree,
            is_program_found_ref    => \$is_program_found,
            is_chain_found_ref      => \$is_chain_found,
            program                 => $start_with_program,
            start_with_programs_ref => \@start_with_programs,
        }
    );

## Then match the expected result
    is_deeply(
        \@start_with_programs,
        \@{$expected_programs_ref},
        q{Start with } . $test_key
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
            store       => \$program_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help     Display this help message
    -v/--version  Display version
END_USAGE
}
