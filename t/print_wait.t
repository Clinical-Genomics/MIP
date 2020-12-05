#!/usr/bin/env perl

#! /usr/bin/env perl

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
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

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
        q{MIP::Processmanagement::Processes} => [qw{ print_wait }],
        q{MIP::Test::Fixtures}               => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Processes qw{ print_wait };

diag(   q{Test print_wait from Processes.pm v}
      . $MIP::Processmanagement::Processes::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $file_content;

## Store file content in memory by using referenced variable
open my $filehandle, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## Given a maximum number of processes
Readonly my $MAX_PROCESS_NUMBER => 2;

## Given a batch counter to scale each batch
my $process_batch_count = 1;

## Given a list of processes
my @process_counters = ( 1, 2 );

my @returned_process_batch_count;

## When iterating over each process
PROCESS_COUNTER:
foreach my $process_counter (@process_counters) {

    push @returned_process_batch_count,
      print_wait(
        {
            filehandle            => $filehandle,
            max_process_number    => $MAX_PROCESS_NUMBER,
            process_batches_count => $process_batch_count,
            process_counter       => $process_counter,
        }
      );
}

close $filehandle;

## Then return $process_batch_count since no print "wait" was needed
is( $returned_process_batch_count[0], $process_batch_count, q{Did not print wait} );

## Then return an incremented $process_batch_count since "wait" for printed
is( $returned_process_batch_count[1], 2, q{Update process_batches_count for next loop} );

done_testing();
