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


## Constants
Readonly my $EXPECTED_COVERAGE => 30;

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Active_parameter} => [qw{ get_user_supplied_pedigree_parameter }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ get_user_supplied_pedigree_parameter };

diag(   q{Test get_user_supplied_pedigree_parameter from Active_parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given user supplied parameters
my %active_parameter = (
    analysis_type     => { sample_id => q{wgs}, },
    sample_ids        => [qw{ sample_1 sample_2 }],
    expected_coverage => $EXPECTED_COVERAGE,
);

my %is_user_supplied = get_user_supplied_pedigree_parameter(
    {
        active_parameter_href => \%active_parameter,
    }
);

## Then return 1 for user supplied parameters
is( $is_user_supplied{analysis_type}, 1, q{Got set HASH parameter} );

is( $is_user_supplied{sample_ids}, 1, q{Got set ARRAY parameter} );

is( $is_user_supplied{expected_coverage}, 1, q{Got set SCALAR parameter} );

## and 0 for parameter which was not supplied by user
is( $is_user_supplied{exome_target_bed}, 0, q{No user defined input for parameter} );

done_testing();
