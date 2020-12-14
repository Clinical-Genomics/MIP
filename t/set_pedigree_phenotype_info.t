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
        q{MIP::Pedigree}       => [qw{ set_pedigree_phenotype_info }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Pedigree qw{ set_pedigree_phenotype_info };

diag(   q{Test set_pedigree_phenotype_info from Pedigree.pm}
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
        parameter_href => \%parameter,
        pedigree_href  => \%pedigree,
    }
);

my @got_affected_samples   = @{ $parameter{cache}{affected} };
my @got_unaffected_samples = @{ $parameter{cache}{unaffected} };
my @got_unknown_samples    = @{ $parameter{cache}{unknown} };

## Then pedigre members phenotypes should have been set in cache
is( scalar @got_affected_samples, 1, q{Got all samples with affected phenotype} );

is( scalar @got_unaffected_samples, 1, q{Got all samples with unaffected phenotype} );

is( scalar @got_unknown_samples, 2, q{Got all samples with unknown phenotype} );

done_testing();
