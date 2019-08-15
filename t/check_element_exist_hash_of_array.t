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
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA   => q{,};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Check::Hash}    => [qw{ check_element_exist_hash_of_array }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Hash qw{ check_element_exist_hash_of_array };

diag(   q{Test check_element_exist_hash_of_array from Hash.pm v}
      . $MIP::Check::Hash::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Create the hash of arrays
my %file_info = (
    exome_target_bed => [qw{ .interval_list .pad100.interval_list }],

    # BWA human genome reference file endings
    bwa_build_reference => [qw{ .amb .ann .bwt .pac .sa }],

    # Human genome meta files
    human_genome_reference_file_endings => [qw{ .dict .fai }],
);

# Test an existing key and an element not part of the array hash
my $not_exist = check_element_exist_hash_of_array(
    {
        element  => q{.another_file_ending},
        hash_ref => \%file_info,
        key      => q{bwa_build_reference},
    }
);

is( $not_exist, 1, q{Element not part of hash of arrays} );

my $exist = check_element_exist_hash_of_array(
    {
        element  => q{.amb},
        hash_ref => \%file_info,
        key      => q{bwa_build_reference},
    }
);

is( $exist, undef, q{Element is part of hash of arrays} );

done_testing();
