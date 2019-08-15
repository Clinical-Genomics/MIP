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
Readonly my $COMMA             => q{,};
Readonly my $EXPECTED_COVERAGE => 30;
Readonly my $SPACE             => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Get::Parameter} => [qw{ get_user_supplied_info }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::Parameter qw{ get_user_supplied_info };

diag(   q{Test get_user_supplied_info from Parameter.pm v}
      . $MIP::Get::Parameter::VERSION
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

my %user_supply_switch = get_user_supplied_info(
    {
        active_parameter_href => \%active_parameter,
    }
);

## Then return 1 for user supplied parameters
is( $user_supply_switch{analysis_type}, 1, q{Got set HASH parameter} );

is( $user_supply_switch{sample_ids}, 1, q{Got set ARRAY parameter} );

is( $user_supply_switch{expected_coverage}, 1, q{Got set SCALAR parameter} );

## and 0 for parameter which was not supplied by user
is( $user_supply_switch{exome_target_bed}, 0, q{No user defined input for parameter} );

done_testing();
