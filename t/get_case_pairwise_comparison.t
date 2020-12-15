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
        q{MIP::Qccollect}      => [qw{ get_case_pairwise_comparison }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qccollect qw{ get_case_pairwise_comparison };

diag(   q{Test get_case_pairwise_comparison from Qccollect.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $IS_SELF                     => 1;
Readonly my $PARENT_CHILD_RELATED_METRIC => 0.76;
Readonly my $UNRELATED_METRIC            => 0.62;
Readonly my $PARENTS_RELATED_METRIC      => 0.8;
Readonly my $SIBLING_RELATED_METRIC      => 0.72;
Readonly my $SIBLING_INDEX               => 11;

## Given
my @relationship_values = (
    $IS_SELF,                     $UNRELATED_METRIC,
    $PARENT_CHILD_RELATED_METRIC, $UNRELATED_METRIC,
    $IS_SELF,                     $PARENT_CHILD_RELATED_METRIC,
    $PARENT_CHILD_RELATED_METRIC, $PARENT_CHILD_RELATED_METRIC,
    $IS_SELF,
);
my @samples = qw{ ADM1059A2 ADM1059A3 ADM1059A1 };

my %case = get_case_pairwise_comparison(
    {
        relationship_values_ref => \@relationship_values,
        sample_orders_ref       => \@samples,
    }
);

my %expected_case = (
    ADM1059A2 => {
        ADM1059A1 => [$PARENT_CHILD_RELATED_METRIC],
        ADM1059A2 => [$IS_SELF],
        ADM1059A3 => [$UNRELATED_METRIC],
    },
    ADM1059A1 => {
        ADM1059A1 => [$IS_SELF],
        ADM1059A2 => [$PARENT_CHILD_RELATED_METRIC],
        ADM1059A3 => [$PARENT_CHILD_RELATED_METRIC],
    },
    ADM1059A3 => {
        ADM1059A1 => [$PARENT_CHILD_RELATED_METRIC],
        ADM1059A2 => [$UNRELATED_METRIC],
        ADM1059A3 => [$IS_SELF],
    },
);

## Then a hash with all per sample pairwise comparison should be created
is_deeply( \%case, \%expected_case, q{Got pairwaise comparison for case} );

done_testing();
