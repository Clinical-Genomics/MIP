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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Environment::Child_process} => [qw{ child_process }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Child_process qw{ child_process };

diag(   q{Test child_process from Child_process.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given proper system call
my $cmd          = q{ls};
my $process_type = q{ipc_cmd_run};

my %process_return = child_process(
    {
        commands_ref => [$cmd],
        process_type => $process_type,
    }
);

## Then capture the output
ok( @{ $process_return{stdouts_ref} }, q{Captured output from system call} );

## Then no errors should be generated
is( @{ $process_return{stderrs_ref} }, 0, q{No error from correct system call} );

## Given system call, when writing to STDERR
$cmd = q{perl -e 'print STDERR q{Ops}'};

%process_return = child_process(
    {
        commands_ref => [$cmd],
        process_type => $process_type,
    }
);

## Then capture the errors
ok( @{ $process_return{stderrs_ref} }, q{Captured error from system call} );

## Then no output should be generated
is( @{ $process_return{stdouts_ref} }, 0, q{No output from incorrect system call} );

done_testing();
