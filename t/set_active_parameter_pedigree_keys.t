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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

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
        q{MIP::Active_parameter} => [qw{ get_user_supplied_pedigree_parameter }],
        q{MIP::Pedigree}         => [qw{ set_active_parameter_pedigree_keys }],
        q{MIP::Test::Fixtures}   => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ get_user_supplied_pedigree_parameter };
use MIP::Pedigree qw{ set_active_parameter_pedigree_keys };

diag(   q{Test set_active_parameter_pedigree_keys from Pedigree.pm v}
      . $MIP::Pedigree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %active_parameter = ( expected_coverage => { sample_1 => 30, }, );

my %pedigree = (
    case    => q{case_1},
    samples => [
        {
            analysis_type => q{wes},
            father        => 0,
            mother        => 0,
            phenotype     => q{affected},
            sample_id     => q{sample_1},
            sex           => q{female},
            time_point    => 0,
        },
        {
            analysis_type => q{wgs},
            father        => 0,
            mother        => 0,
            phenotype     => q{unaffected},
            sample_id     => q{sample_2},
            sex           => q{male},
            time_point    => 1,
        },
        {
            analysis_type => q{wts},
            father        => 0,
            mother        => 0,
            phenotype     => q{unknown},
            sample_id     => q{sample_3},
            sex           => q{other},
            time_point    => 0,
        },
        {
            analysis_type => q{wgs},
            father        => q{sample_1},
            mother        => q{sample_2},
            phenotype     => q{unknown},
            sample_id     => q{sample_4},
            sex           => q{unknown},
            time_point    => 0,
        },
    ],
);

my %sample_info = (
    sample => {
        sample_1 => {
            analysis_type     => q{wes},
            expected_coverage => 30,
        },
        sample_2 => {
            time_point => 1,
        },
    },
);

my %is_user_supplied = get_user_supplied_pedigree_parameter(
    {
        active_parameter_href => \%active_parameter,
    }
);

set_active_parameter_pedigree_keys(
    {
        active_parameter_href => \%active_parameter,
        pedigree_href         => \%pedigree,
        sample_info_href      => \%sample_info,
        is_user_supplied_href => \%is_user_supplied,
    }
);
my $set_analysis_type     = $active_parameter{analysis_type}{sample_1};
my $set_expected_coverage = $sample_info{sample}{expected_coverage}{sample_1};
my $set_time_point        = $sample_info{sample}{sample_2}{time_point};

is( $set_analysis_type,     q{wes}, q{Set analysis type} );
is( $set_expected_coverage, undef,  q(Did not set expected coverage) );
is( $set_time_point, $active_parameter{time_point}{sample_2}, q{Set time point} );

done_testing();
