#!/usr/bin/env perl

#### Copyright 2017 Henrik Stranneheim

use Modern::Perl qw(2014);
use warnings qw(FATAL utf8);
use autodie;
use 5.018;    #Require at least perl 5.18
use utf8;
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use English qw(-no_match_vars);
use Params::Check qw(check allow last_error);

use FindBin qw($Bin);    #Find directory of script
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catfile catdir devnull);
use Getopt::Long;
use Test::More;

use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use Script::Utils qw(help);

our $USAGE = build_usage( {} );

my $VERBOSE = 0;
our $VERSION = '1.0.0';

###User Options
GetOptions(
    'h|help' => sub {
        done_testing();
        print {*STDOUT} $USAGE, "\n";
        exit;
    },    #Display help text
    'v|version' => sub {
        done_testing();
        print {*STDOUT} "\n" . basename($PROGRAM_NAME) . q{  } . $VERSION,
          "\n\n";
        exit;
    },    #Display version number
    'vb|verbose' => $VERBOSE,
  )
  or (
    done_testing(),
    Script::Utils::help(
        {
            USAGE     => $USAGE,
            exit_code => 1,
        }
    )
  );

BEGIN {

### Check all internal dependency modules and imports
    ## Modules with import
    my %perl_module = (
        'Script::Utils'       => [qw{help}],
        'MIP::Check::Cluster' => [qw{check_max_core_number}],
    );

    while ( my ( $module, $module_import ) = each %perl_module ) {

        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load } . $module;
    }

    ## Modules
    my @modules = (qw{MIP::Cluster});

    for my $module (@modules) {

        require_ok($module) or BAIL_OUT q{Cannot load } . $module;
    }
}

use MIP::Cluster qw(get_core_number);

diag(
    "Test get_core_number $MIP::Cluster::VERSION, Perl $^V, $EXECUTABLE_NAME" );

# Core number to test
Readonly my $MODULE_CORE_NUMBER   => 1;
Readonly my $MODIFIER_CORE_NUMBER => 2;
Readonly my $MAX_CORES_PER_NODE   => 2;

my $returned_module_core_number = get_core_number(
    {
        module_core_number   => $MODULE_CORE_NUMBER,
        modifier_core_number => $MODIFIER_CORE_NUMBER,
        max_cores_per_node   => $MAX_CORES_PER_NODE,
    }
);

my $returned_modified_core_number = get_core_number(
    {
        module_core_number   => undef,
        modifier_core_number => $MODIFIER_CORE_NUMBER,
        max_cores_per_node   => $MAX_CORES_PER_NODE,
    }
);

## Test
is( $returned_module_core_number, 1, q{Got module core number} );

is( $returned_modified_core_number, 2, q{Got modifier core number} );

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

##build_usage

##Function : Build the USAGE instructions
##Returns  : ""
##Arguments: $program_name
##         : $program_name => Name of the script

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

    check( $tmpl, $arg_href, 1 ) or croak qw(Could not parse arguments!);

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
