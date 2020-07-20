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
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

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
        q{MIP::Qccollect}      => [qw{ relation_check }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qccollect qw{ relation_check };

diag(   q{Test relation_check from Qccollect.pm v}
      . $MIP::Qccollect::VERSION
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
    $PARENT_CHILD_RELATED_METRIC, $PARENT_CHILD_RELATED_METRIC,
    $UNRELATED_METRIC,            $IS_SELF,
    $PARENT_CHILD_RELATED_METRIC, $PARENT_CHILD_RELATED_METRIC,
    $PARENT_CHILD_RELATED_METRIC, $PARENT_CHILD_RELATED_METRIC,
    $IS_SELF,                     $SIBLING_RELATED_METRIC,
    $PARENT_CHILD_RELATED_METRIC, $PARENT_CHILD_RELATED_METRIC,
    $SIBLING_RELATED_METRIC,      $IS_SELF,
);
my @samples = qw{ ADM1059A2 ADM1059A3 ADM1059A1 sample_4 };
my %qc_data = (
    recipe => {
        relation_check => { sample_relation_check => [@relationship_values], },
        pedigree_check => { sample_order          => [@samples], },
    },
);
my %sample_info = test_mip_hashes(
    {
        mip_hash_name => q{qc_sample_info},
        recipe_name   => q{a_recipe},
    }
);

$sample_info{sample}{sample_4} = $sample_info{sample}{ADM1059A1};
$sample_info{sample}{sample_4}{sample_id} = q{sample_4};

relation_check(
    {
        qc_data_href            => \%qc_data,
        relationship_values_ref => \@relationship_values,
        sample_info_href        => \%sample_info,
        sample_orders_ref       => \@{ $qc_data{recipe}{pedigree_check}{sample_order} },
    }
);

SAMPLE:
foreach my $sample_id (@samples) {

## Then all relations should pass
    is( $qc_data{sample}{$sample_id}{relation_check},
        q{PASS}, q{Checked relations for pedigree} );
}

## Given duplicate samples
$relationship_values[1] = $IS_SELF;

relation_check(
    {
        qc_data_href            => \%qc_data,
        relationship_values_ref => \@relationship_values,
        sample_info_href        => \%sample_info,
        sample_orders_ref       => \@{ $qc_data{recipe}{pedigree_check}{sample_order} },
    }
);

## Then fail with duplicates message
is(
    $qc_data{sample}{ADM1059A2}{relation_check},
    q{FAIL: Duplicated sample?},
    q{Found duplicate sample}
);

## Given inbreed parents
$relationship_values[1] = $PARENTS_RELATED_METRIC;

relation_check(
    {
        qc_data_href            => \%qc_data,
        relationship_values_ref => \@relationship_values,
        sample_info_href        => \%sample_info,
        sample_orders_ref       => \@{ $qc_data{recipe}{pedigree_check}{sample_order} },
    }
);

## Then fail with parents related message
is(
    $qc_data{sample}{ADM1059A2}{relation_check},
    q{FAIL: Parents related?},
    q{Found related parents}
);

## Given unrelated siblings
# Reset to PASS value
$relationship_values[1] = $UNRELATED_METRIC;

# Add unrelated sibling
$relationship_values[$SIBLING_INDEX] = $UNRELATED_METRIC;

relation_check(
    {
        qc_data_href            => \%qc_data,
        relationship_values_ref => \@relationship_values,
        sample_info_href        => \%sample_info,
        sample_orders_ref       => \@{ $qc_data{recipe}{pedigree_check}{sample_order} },
    }
);

## Then fail with sample not related message
is(
    $qc_data{sample}{ADM1059A1}{relation_check},
    q{FAIL:ADM1059A1 not related to sample_4},
    q{Found unrelated siblings}
);

## When no relationship_values and no sample_order
delete $qc_data{recipe}{relation_check}{sample_relation_check};
delete $qc_data{recipe}{pedigree_check}{sample_order};

my $return = relation_check(
    {
        qc_data_href => \%qc_data,
        relationship_values_ref =>
          \@{ $qc_data{recipe}{relation_check}{sample_relation_check} },
        sample_info_href  => \%sample_info,
        sample_orders_ref => \@{ $qc_data{recipe}{pedigree_check}{sample_order} },
    }
);

## Then skip relation check and return undef
is( $return, undef, q{Skip if no sample order and no relationship values} );

done_testing();
