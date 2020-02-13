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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

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
        q{MIP::Qccollect}      => [qw{ plink_gender_check }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qccollect qw{ plink_gender_check };

diag(   q{Test plink_gender_check from Qccollect.pm v}
      . $MIP::Qccollect::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given pedigree, qc_data gender and sample_info sex, when represented by letters
my %expected_qc_data = (
    recipe => {
        plink_sexcheck =>
          { plink_gender_check => [qw{ ADM1059A1:PASS ADM1059A2:PASS ADM1059A3:PASS }], },
    },
);
my %qc_data = (
    recipe => {
        plink_sexcheck =>
          { sample_sexcheck => [qw{ ADM1059A1:2 ADM1059A2:1 ADM1059A3:0 }], },
    },
);

my %sample_info = test_mip_hashes(
    {
        mip_hash_name => q{qc_sample_info},
        recipe_name   => q{a_recipe},
    }
);

plink_gender_check(
    {
        qc_data_href     => \%qc_data,
        sample_info_href => \%sample_info,
    }
);

## Then PASS should be added to qc_data for sample
is_deeply(
    \@{ $qc_data{recipe}{plink_sexcheck}{plink_gender_check} },
    \@{ $expected_qc_data{recipe}{plink_sexcheck}{plink_gender_check} },
    q{Passed gender as numbers}
);

## Given pedigree and sample_info sex, when represented by numbers
# Modify to plink format for gender
$sample_info{sample}{ADM1059A1}{sex} = q{female};
$sample_info{sample}{ADM1059A2}{sex} = q{male};
$sample_info{sample}{ADM1059A3}{sex} = q{unknown};

# Remove prior plink gender check
delete $qc_data{recipe}{plink_sexcheck}{plink_gender_check};

plink_gender_check(
    {
        qc_data_href     => \%qc_data,
        sample_info_href => \%sample_info,
    }
);

## Then PASS should be added to qc_data for sample
is_deeply(
    \@{ $qc_data{recipe}{plink_sexcheck}{plink_gender_check} },
    \@{ $expected_qc_data{recipe}{plink_sexcheck}{plink_gender_check} },
    q{Passed as letters}
);

## Given an incorrect gender
my $failed_sample_id = q{ADM1059A1};
$qc_data{recipe}{plink_sexcheck}{sample_sexcheck} =
  [qw{ ADM1059A1:1 ADM1059A2:1 ADM1059A3:0 }];

# Remove prior plink gender check
delete $qc_data{recipe}{plink_sexcheck}{plink_gender_check};
$expected_qc_data{recipe}{plink_sexcheck}{plink_gender_check} =
  [qw{ ADM1059A1:FAIL ADM1059A2:PASS ADM1059A3:PASS }];

plink_gender_check(
    {
        qc_data_href     => \%qc_data,
        sample_info_href => \%sample_info,
    }
);

is_deeply(
    \@{ $qc_data{recipe}{plink_sexcheck}{plink_gender_check} },
    \@{ $expected_qc_data{recipe}{plink_sexcheck}{plink_gender_check} },
    q{Failed gender }
);

## Given an incorrect metrics when trailing garbage
$qc_data{recipe}{plink_sexcheck}{sample_sexcheck} =
  [qw{ ADM1059A1:1:garbage ADM1059A2:1 ADM1059A3:0 }];

trap {
    plink_gender_check(
        {
            qc_data_href     => \%qc_data,
            sample_info_href => \%sample_info,
        }
    )
};
## Then exit and throw FATAL log message
like( $trap->stderr, qr/Unexpected\s+trailing/xms, q{Throw warning log message} );

done_testing();
