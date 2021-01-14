#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
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
use MIP::Constants qw{ $COMMA $NEWLINE $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::List}           => [qw{ check_element_exist_hash_of_array }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::List qw{ check_element_exist_hash_of_array };

diag(   q{Test check_element_exist_hash_of_array from List.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a hash of arrays
my %file_info = (
    exome_target_bed => [qw{ .interval_list .pad100.interval_list }],

    # BWA human genome reference file endings
    bwa_build_reference => [qw{ .amb .ann .bwt .pac .sa }],

    # Human genome meta files
    human_genome_reference_file_endings => [qw{ .dict .fai }],
);

## When a key that is not part of hash is used
my $return = check_element_exist_hash_of_array(
    {
        element  => q{.another_file_ending},
        hash_ref => \%file_info,
        key      => q{not_a_key},
    }
);

## Then return zero
is( $return, 0, q{Key is not part of hash} );

## Given an existing key and a missing element

## When not part of the array hash
my $not_exist = check_element_exist_hash_of_array(
    {
        element  => q{.another_file_ending},
        hash_ref => \%file_info,
        key      => q{bwa_build_reference},
    }
);

## Then return undef
is( $not_exist, undef, q{Element not part of hash of arrays} );

## When part of array
my $exist = check_element_exist_hash_of_array(
    {
        element  => q{.amb},
        hash_ref => \%file_info,
        key      => q{bwa_build_reference},
    }
);

## Then return true
is( $exist, 1, q{Element is part of hash of arrays} );

done_testing();
