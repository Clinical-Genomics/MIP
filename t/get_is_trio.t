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
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Pedigree}       => [qw{ get_is_trio }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Pedigree qw{ get_is_trio };

diag(   q{Test get_is_trio from Pedigree.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 1, } );

## Given only single sample
my %active_parameter = ( sample_ids => [qw{ sample_1 }], );
my %sample_info;

my $is_trio = get_is_trio(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
);

## Then do not detect trio
is( $is_trio, undef, q{Single sample - did not detect trio} );

## Given sample info when 3 children are present and no parents
%active_parameter = ( sample_ids => [qw{ child_1 child_2 child_3 }], );
%sample_info      = (
    sample => {
        child_1 => {
            father => 0,
            mother => 0,
        },
        child_2 => {
            father => 0,
            mother => 0,
        },
        child_3 => {
            father => 0,
            mother => 0,
        },
    },
);

$is_trio = get_is_trio(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
);

## Then do not detect trio
is( $is_trio, undef, q{Three children - did not detect trio} );

## Given more samples than a trio, when one child has parents, but not in analysis
%sample_info = (
    sample => {
        child_1 => {
            father => 1,
            mother => 2,
        },
        child_2 => {
            father => 0,
            mother => 0,
        },
        child_3 => {
            father => 0,
            mother => 0,
        },
    },
);

$is_trio = get_is_trio(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
);

## Then do not detect trio
is( $is_trio, undef,
    q{Three children where one is trio but not in analysis - did not detect trio} );

## Given more samples than a trio, when one child has parents in analysis
%active_parameter = ( sample_ids => [qw{ child_1 child_2 child_3 father_1 mother_1 }], );

%sample_info = (
    sample => {
        child_1 => {
            father => q{father_1},
            mother => q{mother_2},
        },
        child_2 => {
            father => 0,
            mother => 0,
        },
        child_3 => {
            father => 0,
            mother => 0,
        },
    },
);

$is_trio = get_is_trio(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
);

## Then do not detect trio due to to many samples
is( $is_trio, undef,
    q{Three children where one is trio in analysis - did not detect trio} );

## Given a trio, when correct number of samples in analysis
%active_parameter = ( sample_ids => [qw{ child_1 father_1 mother_1 }], );

%sample_info = (
    sample => {
        child_1 => {
            father => q{father_1},
            mother => q{mother_1},
        },
        father_1 => {
            father => 0,
            mother => 0,
        },
        mother_1 => {
            father => 0,
            mother => 0,
        },
    },
);

$is_trio = get_is_trio(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
);

## Then do detect trio
is( $is_trio, 1, q{Trio - did detect trio} );

done_testing();
