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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Pedigree}       => [qw{ detect_founders }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Pedigree qw{ detect_founders };

diag(   q{Test detect_founders from Pedigree.pm v}
      . $MIP::Pedigree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given only single sample
my %active_parameter = ( sample_ids => [qw{ sample_1 }], );
my %sample_info;

my $founders_count = detect_founders(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
);

## Then do not detect trio
is( $founders_count, 0, q{Single sample - did not detect founders} );

## Given more samples than a trio, when one child has a single parent in analysis
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

$founders_count = detect_founders(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
);

## Then do find one parent due to to many samples
is( $founders_count, 1,
    q{Three children where one is trio with one parent in analysis - did detect founders}
);

## Given a correct trio, when correct number of samples in analysis
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

$founders_count = detect_founders(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
);

## Then do detect trio
is( $founders_count, 2, q{Correct trio - did detect all founders} );

done_testing();
