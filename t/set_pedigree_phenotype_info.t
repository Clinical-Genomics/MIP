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
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Set::Pedigree}  => [qw{ set_pedigree_phenotype_info }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::Pedigree qw{ set_pedigree_phenotype_info };

diag(   q{Test set_pedigree_phenotype_info from Pedigree.pm v}
      . $MIP::Set::Pedigree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a pedigree with phenotype info
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
        },
        {
            analysis_type => q{wgs},
            father        => 0,
            mother        => 0,
            phenotype     => q{unaffected},
            sample_id     => q{sample_2},
            sex           => q{male},
        },
        {
            analysis_type => q{wts},
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
my %parameter;
set_pedigree_phenotype_info(
    {
        pedigree_href  => \%pedigree,
        parameter_href => \%parameter,
    }
);

my @got_unaffected_samples = @{ $parameter{cache}{unaffected} };
my @got_affected_samples   = @{ $parameter{cache}{affected} };
my @got_unknown_samples    = @{ $parameter{cache}{unknown} };

## Then pedigre members phenotypes should have been set in cache
is( scalar @got_unaffected_samples, 1, q{Got all samples with unaffected phenotype} );

is( scalar @got_affected_samples, 1, q{Got all samples with affected phenotype} );

is( scalar @got_unknown_samples, 2, q{Got all samples with unknown phenotype} );

## As well as phenoype info in Plink format
my $affected_plink_phenotype =
  $parameter{cache}{sample_1}{plink_phenotype};
my $unaffected_plink_phenotype =
  $parameter{cache}{sample_2}{plink_phenotype};
my $unknown_plink_phenotype =
  $parameter{cache}{sample_3}{plink_phenotype};

is( $affected_plink_phenotype, 2, q{Reformated to plink affected phenotype} );

is( $unaffected_plink_phenotype, 1, q{Reformated to plink unaffected phenotype} );

is( $unknown_plink_phenotype, 0, q{Reformated to plink unknown phenotype} );

done_testing();
