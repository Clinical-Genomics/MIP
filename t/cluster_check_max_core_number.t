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
    my %perl_module;

    $perl_module{'Script::Utils'} = [qw(help)];

    while ( my ( $module, $module_import ) = each %perl_module ) {

        use_ok( $module, @{$module_import} )
          or BAIL_OUT 'Cannot load ' . $module;
    }

    ## Modules
    my @modules = ('MIP::Check::Cluster');

    for my $module (@modules) {

        require_ok($module) or BAIL_OUT 'Cannot load ' . $module;
    }
}

use MIP::Check::Cluster qw(check_max_core_number);

diag(
"Test check_max_core_number $MIP::Check::Cluster::VERSION, Perl $^V, $EXECUTABLE_NAME"
);

# Core number to test
Readonly my $LOWER_THAN_MAX_CORES_PER_NODE    => 1;
Readonly my $EQUALS_MAX_CORES_PER_NODE        => 2;
Readonly my $GREATHER_THAN_MAX_CORES_PER_NODE => 3;

my @test_core_numbers = (
    $LOWER_THAN_MAX_CORES_PER_NODE,
    $EQUALS_MAX_CORES_PER_NODE, $GREATHER_THAN_MAX_CORES_PER_NODE,
);

# Possibly adjusted core numbers according to max core numbers
my @returned_core_numbers;

foreach my $core_number (@test_core_numbers) {

    push @returned_core_numbers,
      check_max_core_number(
        {
            max_cores_per_node    => 2,
            core_number_requested => $core_number,
        }
      );
}

## Test
is(
    $returned_core_numbers[0],
    $LOWER_THAN_MAX_CORES_PER_NODE,
    'Core number requested is lower than max core number'
);

is( $returned_core_numbers[1],
    $EQUALS_MAX_CORES_PER_NODE,
    'Core number requested equals max core number' );

is( $returned_core_numbers[2],
    $EQUALS_MAX_CORES_PER_NODE,
    'Core number requested was greather than max core number' );

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
