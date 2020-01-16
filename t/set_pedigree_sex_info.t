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
        q{MIP::Pedigree}       => [qw{ set_pedigree_sex_info }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Pedigree qw{ set_pedigree_sex_info };

diag(   q{Test set_pedigree_sex_info from Pedigree.pm v}
      . $MIP::Pedigree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a pedigree
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
set_pedigree_sex_info(
    {
        parameter_href => \%parameter,
        pedigree_href  => \%pedigree,
    }
);

my @got_female_samples  = @{ $parameter{cache}{female} };
my @got_male_samples    = @{ $parameter{cache}{male} };
my @got_other_samples   = @{ $parameter{cache}{other} };
my @got_unknown_samples = @{ $parameter{cache}{unknown} };

## Then pedigree members sex should have been set
is( scalar @got_female_samples, 1, q{Got all samples with female sex} );

is( scalar @got_male_samples, 1, q{Got all samples with male sex} );

is( scalar @got_other_samples, 1, q{Got all samples with other sex} );

is( scalar @got_unknown_samples, 1, q{Got all samples with unknown sex} );

## As well as memebers sex in plink format
my $female_plink_sex  = $parameter{cache}{sample_1}{plink_sex};
my $male_plink_sex    = $parameter{cache}{sample_2}{plink_sex};
my $other_plink_sex   = $parameter{cache}{sample_3}{plink_sex};
my $unknown_plink_sex = $parameter{cache}{sample_4}{plink_sex};

is( $female_plink_sex, 2, q{Reformated to plink female sex} );

is( $male_plink_sex, 1, q{Reformated to plink male sex} );

is( $other_plink_sex, q{other}, q{Reformated to plink other sex} );

is( $unknown_plink_sex, q{other}, q{Reformated unknown to plink other sex} );

done_testing();
