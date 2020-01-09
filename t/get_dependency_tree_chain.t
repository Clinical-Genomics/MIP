#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use Getopt::Long;
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
        say {*STDOUT} $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $VERSION . $NEWLINE;
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

use MIP::Get::Analysis qw{ get_dependency_tree_chain };

diag(   q{Test get_dependency_tree_chain from Analysis.pm v}
      . $MIP::Get::Analysis::VERSION
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
