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
our $VERSION = '1.0.0';

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

my %dependency_tree = (
    initiation => [
        q{program_0},
        q{program_1},
        {
            CHAIN_0 => {
                PARALLEL => [
                    q{parallel_program_0},
                    q{parallel_program_1},
                    {
                        parallel_program_2 =>
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
                        parallel_program_6 =>
                          [ q{parallel_program_6}, q{parallel_program_7}, ],
                    },
                ],
            },
        },
        q{program_3},
    ],
);
## Expected results from traversing the dependency tree
my @expected_programs_0 =
  qw{ program_0 program_1 parallel_program_0 parallel_program_1 parallel_program_2 parallel_program_3 program_2 parallel_program_4 parallel_program_5 parallel_program_6 parallel_program_7 program_3 };
my @expected_parallel_programs_0 = qw{ parallel_program_0 program_2 program_3 };
my @expected_parallel_programs_2 =
  qw{ parallel_program_2 parallel_program_3 program_2 program_3 };
my @expected_parallel_programs_3 = qw{ parallel_program_3 program_2 program_3 };
my @expected_program_2 =
  qw{ program_2 parallel_program_4 parallel_program_5 parallel_program_6 parallel_program_7 program_3 program_3};
my @expected_parallel_programs_5 = qw{ parallel_program_5 program_3 };
my @expected_program_3           = qw{ program_3 };

## Tests to run
my %test = (
    program_0          => \@expected_programs_0,
    parallel_program_0 => \@expected_parallel_programs_0,
    parallel_program_2 => \@expected_parallel_programs_2,
    parallel_program_3 => \@expected_parallel_programs_3,
    program_2          => \@expected_program_2,
    parallel_program_5 => \@expected_parallel_programs_5,
    program_3          => \@expected_program_3,
);

## Test
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

    is( array_diff( @start_with_programs, @{$expected_programs_ref} ),
        0, q{Start with } . $test_key );
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
