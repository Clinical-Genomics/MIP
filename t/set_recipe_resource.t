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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Active_parameter} => [qw{ set_recipe_resource }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ set_recipe_resource };

diag(   q{Test set_recipe_resource from Active_parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $BWA_MEM_CORE_NR => 16;
Readonly my $BWA_MEM_TIME    => 5;
Readonly my $MANTA_CORE_NR   => 36;
Readonly my $SET_BWA_CORE_NR => 36;
Readonly my $SET_MANTA_TIME  => 5;

my %active_parameter = (
    recipe_core_number => {
        bwa_mem => $BWA_MEM_CORE_NR,
        fastqc  => 2,
        manta   => $MANTA_CORE_NR,
    },
    recipe_time => {
        bwa_mem => $BWA_MEM_TIME,
        fastqc  => 1,
        manta   => 1,
    },
    set_recipe_core_number => {
        bwa_mem => $SET_BWA_CORE_NR,
        fastqc  => 1,
    },
    set_recipe_time => { manta => $SET_MANTA_TIME, },
);

## Given
set_recipe_resource( { active_parameter_href => \%active_parameter, } );

my %expected_recipe_parameter = (
    recipe_core_number => {
        bwa_mem => $SET_BWA_CORE_NR,
        fastqc  => 1,
        manta   => $MANTA_CORE_NR,
    },
    recipe_time => {
        bwa_mem => $BWA_MEM_TIME,
        fastqc  => 1,
        manta   => $SET_MANTA_TIME,
    },
);

## Then the resource recipe_core_number for the recipes should have been set
is_deeply(
    \%{ $active_parameter{recipe_core_number} },
    \%{ $expected_recipe_parameter{recipe_core_number} },
    q{Set recipes core number parameters}
);

## Then the resource recipe_time for the recipes should have been set
is_deeply(
    \%{ $active_parameter{recipe_time} },
    \%{ $expected_recipe_parameter{recipe_time} },
    q{Set recipes time parameters}
);

done_testing();
