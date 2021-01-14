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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Pedigree}       => [qw{ set_pedigree_capture_kit_info }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Pedigree qw{ set_pedigree_capture_kit_info };

diag(   q{Test set_pedigree_capture_kit_info from Pedigree.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %active_parameter;

my %parameter = (
    supported_capture_kit => {
        default => {
            q{agilent_sureselect.v4} =>
              q{genome_reference_source_version_agilent_sureselect_targets_-v4-.bed},
            q{agilent_sureselect.v5} =>
              q{genome_reference_source_version_agilent_sureselect_targets_-v5-.bed},
        },
    },
);

my %pedigree = (
    case    => q{case_1},
    samples => [
        {
            analysis_type => q{wes},
            capture_kit   => q{agilent_sureselect.v5},
            father        => 0,
            mother        => 0,
            phenotype     => q{affected},
            sample_id     => q{sample_1},
            sex           => q{female},
        },
        {
            analysis_type => q{wgs},
            capture_kit   => q{agilent_sureselect.v4},
            father        => 0,
            mother        => 0,
            phenotype     => q{unaffected},
            sample_id     => q{sample_2},
            sex           => q{male},
        },
        {
            analysis_type => q{wts},
            capture_kit   => q{agilent_sureselect.v5},
            father        => 0,
            mother        => 0,
            phenotype     => q{unknown},
            sample_id     => q{sample_3},
            sex           => q{other},
        },
        {
            analysis_type => q{wgs},
            father        => q{sample_1},
            mother        => q{sample_2},
            phenotype     => q{unknown},
            sample_id     => q{sample_4},
            sex           => q{unknown},
        },
    ],
);

my %sample_info = (
    sample => {
        sample_1 => {
            analysis_type     => q{wes},
            expected_coverage => 30,
            capture_kit       => q{agilent_sureselect.v5},
        },
        sample_2 => {
            analysis_type     => q{wes},
            expected_coverage => 30,
            capture_kit       => q{agilent_sureselect.v4},
        },
        sample_3 => {
            analysis_type     => q{wes},
            expected_coverage => 30,
            capture_kit       => q{agilent_sureselect.v5},
        },
    },
);

my %user_supply_switch = ( exome_target_bed => 0, );
set_pedigree_capture_kit_info(
    {
        active_parameter_href => \%active_parameter,
        is_user_supplied_href => \%user_supply_switch,
        parameter_href        => \%parameter,
        pedigree_href         => \%pedigree,
        sample_info_href      => \%sample_info,
    }
);

my $capture_kit_string = $active_parameter{exome_target_bed}
  {q{genome_reference_source_version_agilent_sureselect_targets_-v5-.bed}};
is( $capture_kit_string, q{sample_1,sample_3}, q{Set sample_ids for capture kit 1} );

$capture_kit_string = $active_parameter{exome_target_bed}
  {q{genome_reference_source_version_agilent_sureselect_targets_-v4-.bed}};
is( $capture_kit_string, q{sample_2}, q{Set sample_ids for capture kit 2} );

## Test unknown capture kit
%active_parameter = ();
%sample_info      = (
    sample => {
        sample_1 => {
            analysis_type     => q{wes},
            expected_coverage => 30,
            capture_kit       => q{unknown_capture_kit},
        },
    },
);

set_pedigree_capture_kit_info(
    {
        active_parameter_href => \%active_parameter,
        is_user_supplied_href => \%user_supply_switch,
        parameter_href        => \%parameter,
        pedigree_href         => \%pedigree,
        sample_info_href      => \%sample_info,
    }
);
$capture_kit_string = $active_parameter{exome_target_bed}{q{unknown_capture_kit}};
is( $capture_kit_string, q{sample_1}, q(Unknown capture kit) );

## Test no capture kit
%active_parameter = ();
%sample_info      = ();

set_pedigree_capture_kit_info(
    {
        active_parameter_href => \%active_parameter,
        is_user_supplied_href => \%user_supply_switch,
        parameter_href        => \%parameter,
        pedigree_href         => \%pedigree,
        sample_info_href      => \%sample_info,
    }
);

is( $active_parameter{exome_target_bed},
    undef, q(No capture kit from cmd, config or pedigree) );

done_testing();
