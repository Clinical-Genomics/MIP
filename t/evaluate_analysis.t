#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

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
        q{MIP::Qccollect}      => [qw{ evaluate_analysis }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qccollect qw{ evaluate_analysis };

diag(   q{Test evaluate_analysis from Qccollect.pm v}
      . $MIP::Qccollect::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

Readonly my $FRACTION_OF_COMMON_VARIANTS => 0.3;

my $log = test_log( {} );

## Given
my $metric_lt = q{fraction_of_common_variants};
my $recipe_lt = q{variant_integrity_ar_father};
my %qc_data   = (
    recipe => {
        $recipe_lt => {
            $metric_lt => $FRACTION_OF_COMMON_VARIANTS
        },
    },
);

my %sample_info = test_mip_hashes(
    {
        mip_hash_name => q{qc_sample_info},
    }
);

my $eval_metric_file =
  catfile( dirname($Bin), qw{ t data references qc_eval_metric_-v1.0-.yaml} );

## Then set the relevant evaluation metrics for the analysis
my $is_ok = evaluate_analysis(
    {
        eval_metric_file => $eval_metric_file,
        qc_data_href     => \%qc_data,
        sample_info_href => \%sample_info,
        skip_evaluation  => 0,
    }
);
ok( $is_ok, q{Evaluated analysis} );

## Given request to skip analysis
trap {
    evaluate_analysis(
        {
            eval_metric_file => $eval_metric_file,
            qc_data_href     => \%qc_data,
            sample_info_href => \%sample_info,
            skip_evaluation  => 1,
        }
    )
};
## Then skip evaluation
like( $trap->stderr, qr/Skipping\sQC\sevaluation/xms, q{Skip QC evaluation} );

## Given an eval metric file that doesn't match the analysis
my $eval_metric_file_no_wgs =
  catfile( dirname($Bin), qw{ t data references qc_eval_metric_no_wgs_-v1.0-.yaml} );
$sample_info{sample}{ADM1059A1}{expected_coverage} = undef;
$sample_info{sample}{ADM1059A2}{expected_coverage} = undef;
$sample_info{sample}{ADM1059A3}{expected_coverage} = undef;

trap {
    evaluate_analysis(
        {
            eval_metric_file => $eval_metric_file_no_wgs,
            qc_data_href     => \%qc_data,
            sample_info_href => \%sample_info,
            skip_evaluation  => 0,
        }
    )
};
## Then skip evaluation
like(
    $trap->stderr,
    qr/Supplied\seval_metric_file/xms,
    q{No relevant evaluation metrics}
);

done_testing();
