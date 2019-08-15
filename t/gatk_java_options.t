#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname  };
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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {
    use MIP::Test::Fixtures qw{ test_import };
    ### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Test::Fixtures} => [qw{ test_standard_cli }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Base::Gatk qw{ gatk_java_options };

diag(   q{Test gatk_java_options from Base::Gatk.pm v}
      . $MIP::Program::Base::Gatk::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my @function_base_commands = qw{ program };

my $commands_ref = [qw{ gatk test }];

## Given undefined java options
my $java_use_large_pages;
my $memory_allocation;

## When the subroutine is executed
gatk_java_options(
    {
        commands_ref         => $commands_ref,
        java_use_large_pages => $java_use_large_pages,
        memory_allocation    => $memory_allocation,
    }
);

## Return non modified command array
is_deeply( $commands_ref, [qw{ gatk test }], q{No java options} );

## Given large pages and memory allocation command
$java_use_large_pages = 1;
$memory_allocation    = q{Xmx4G};

## When the subroutine is executed
gatk_java_options(
    {
        commands_ref         => $commands_ref,
        java_use_large_pages => $java_use_large_pages,
        memory_allocation    => $memory_allocation,
    }
);

## Return a modified command array
my $expected_commands_ref =
  [ qw{ gatk test}, q{--java-options "-XX:-UseLargePages -Xmx4G"} ];
is_deeply( $commands_ref, $expected_commands_ref, q{Java options} );

## Given that the program will be executed via xargs
my $xargs_mode = 1;
$commands_ref = [qw{ gatk test }];

## When the subroutine is executed
gatk_java_options(
    {
        commands_ref         => $commands_ref,
        java_use_large_pages => $java_use_large_pages,
        memory_allocation    => $memory_allocation,
        xargs_mode           => $xargs_mode,
    }
);

## Return a modified command array with escaped quotation marks
$expected_commands_ref =
  [ qw{ gatk test}, q{--java-options \"-XX:-UseLargePages -Xmx4G\"} ];
is_deeply( $commands_ref, $expected_commands_ref, q{Xargs_mode} );

done_testing();
