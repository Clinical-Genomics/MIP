#!/usr/bin/env perl

#### Copyright 2017 Henrik Stranneheim

use Modern::Perl qw{ 2018 };
use warnings qw(FATAL utf8);
use autodie;
use 5.026;    # Require at least perl 5.18
use utf8;
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use English qw(-no_match_vars);
use Params::Check qw(check allow last_error);

use FindBin qw($Bin);    # Find directory of script
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catfile catdir devnull);
use Getopt::Long;
use Test::More;

use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Script::Utils qw(help);

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
        print {*STDOUT} "\n" . basename($PROGRAM_NAME) . q{  } . $VERSION, "\n\n";
        exit;
    },    #Display version number
    'vb|verbose' => $VERBOSE,
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
    my %perl_module = ( 'MIP::Script::Utils' => [qw{help}], );

    while ( my ( $module, $module_import ) = each %perl_module ) {

        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load } . $module;
    }

    ## Modules
    my @modules = (qw{MIP::Processmanagement::Processes});

    for my $module (@modules) {

        require_ok($module) or BAIL_OUT q{Cannot load } . $module;
    }
}

use MIP::Processmanagement::Processes qw(print_wait);

diag(
"Test print_wait $MIP::Processmanagement::Processes::VERSION, Perl $^V, $EXECUTABLE_NAME"
);

my $filehandle = IO::Handle->new();

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open $filehandle, '>', \$file_content
  or croak 'Cannot write to ' . $file_content . ': ' . $OS_ERROR;

# Process number to test
Readonly my $MAX_PROCESS_NUMBER => 2;

my $process_batch_count = 1;

my @process_counters = ( 1, 2 );

my @returned_process_batch_count;

foreach my $process_counter (@process_counters) {

    push @returned_process_batch_count,
      print_wait(
        {
            process_counter       => $process_counter,
            max_process_number    => $MAX_PROCESS_NUMBER,
            process_batches_count => $process_batch_count,
            filehandle            => $filehandle,
        }
      );
}

close $filehandle;

## Test
is( $returned_process_batch_count[0], 1, q{Did not print wait} );

is( $returned_process_batch_count[1], 2, q{Update process_batches_count for next loop} );

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
